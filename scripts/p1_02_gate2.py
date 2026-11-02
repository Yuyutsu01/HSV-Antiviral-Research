#!/usr/bin/env python
"""
Phase 1, step 2 - fetch reference proteins and run Gate 2.

Gate 2 asks whether the reference residue named in each HerpesDRG label matches
the residue actually present at that position in the NCBI reference protein.
HerpesDRG states that it maps all mutations to the NCBI reference strain, so
this is confirmation rather than repair - but it is the check that catches an
off-by-one or a wrong reference strain before any analysis depends on it.

Threshold, declared in the protocol before the data was seen: >95% agreement.

Outputs:
  data/raw/reference/<virus>_<gene>.faa      reference protein, cached
  data/processed/gate2_coordinate_check.tsv  per-label verdict
  reports/phase1/gate2_summary.json
"""

import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF_DIR = ROOT / "data" / "raw" / "reference"
OUT = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase1"
LABELS = OUT / "hsv_mutations_tiered.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
DELAY = 0.12 if API_KEY else 0.36

GATE_2_THRESHOLD = 0.95

REFSEQ = {"HSV1": "NC_001806", "HSV2": "NC_001798"}

# Gene synonyms as they appear in [gene=] and [protein=] qualifiers. Matching is
# deliberately explicit: "DNA polymerase processivity factor" is UL42 and must
# never be captured by the UL30 pattern.
GENES = {
    "UL23": {
        "gene": {"ul23", "tk"},
        "protein_re": r"thymidine kinase",
    },
    "UL30": {
        "gene": {"ul30", "pol"},
        "protein_re": r"dna polymerase(?!.*processivity)(?!.*accessory)",
    },
    # UL5 and UL52 are two subunits of the same complex, and UL8 is the third.
    # Names like "DNA helicase/primase complex-associated protein" (UL8) and the
    # bare "helicase-primase subunit" are ambiguous, so they are left unmatched
    # rather than guessed. UL9 is the origin-binding helicase and must not be
    # captured here either.
    "UL5": {
        "gene": {"ul5"},
        "protein_re": r"helicase[-/ ]primase helicase subunit|^dna helicase$",
    },
    "UL52": {
        "gene": {"ul52"},
        "protein_re": r"helicase[-/ ]primase primase subunit|^dna primase$",
    },
}


def eutil(endpoint, params, retries=5):
    p = {"tool": "hsv-resistance-protocol"}
    if API_KEY:
        p["api_key"] = API_KEY
    p.update(params)
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(p)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:                       # noqa: BLE001 - retry all
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{endpoint} failed: {last}")


def parse_cds_fasta(text):
    """Yield (header, sequence) from an NCBI fasta_cds_* payload."""
    header, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                yield header, "".join(buf)
            header, buf = line, []
        elif line.strip():
            buf.append(line.strip())
    if header:
        yield header, "".join(buf)


def qualifiers(header):
    return {m.group(1).lower(): m.group(2)
            for m in re.finditer(r"\[([A-Za-z_]+)=([^\]]*)\]", header)}


def match_gene(q):
    gene = q.get("gene", "").strip().lower()
    protein = q.get("protein", "").strip().lower()
    for target, spec in GENES.items():
        if gene in spec["gene"]:
            return target
        if gene:
            # An explicit but different gene name is authoritative; do not fall
            # through to a loose protein-name match that could steal a paralog.
            continue
        if protein and re.search(spec["protein_re"], protein):
            return target
    # Second pass for records whose gene qualifier is absent entirely.
    if not gene and protein:
        for target, spec in GENES.items():
            if re.search(spec["protein_re"], protein):
                return target
        # Many submissions give the product as the bare locus name --
        # [protein=UL23] with no gene qualifier at all -- which the descriptive
        # patterns above do not match even though it is completely
        # unambiguous. Skipping these cost 244 accessions, most of the records
        # that appeared to carry no target CDS, and hit HSV-2 hardest.
        #
        # Exact match only. A substring test would let UL5 swallow UL50 to
        # UL56, which are six unrelated genes that between them appear in
        # nearly every record.
        for target, spec in GENES.items():
            if protein in spec["gene"]:
                return target
    return None


def load_reference_proteins():
    REF_DIR.mkdir(parents=True, exist_ok=True)
    refs = {}
    for virus, acc in REFSEQ.items():
        cache = REF_DIR / f"{virus}_cds_aa.fasta"
        if cache.exists() and cache.stat().st_size > 0:
            text = cache.read_text(encoding="utf-8")
        else:
            text = eutil("efetch", {"db": "nuccore", "id": acc,
                                    "rettype": "fasta_cds_aa",
                                    "retmode": "text"})
            cache.write_text(text, encoding="utf-8")
            time.sleep(DELAY)
        for header, seq in parse_cds_fasta(text):
            g = match_gene(qualifiers(header))
            if g and (virus, g) not in refs:
                refs[(virus, g)] = seq.rstrip("*")
                (REF_DIR / f"{virus}_{g}.faa").write_text(
                    f">{virus}_{g} {REFSEQ[virus]}\n{seq.rstrip('*')}\n",
                    encoding="utf-8")
    return refs


def main():
    REP.mkdir(parents=True, exist_ok=True)
    refs = load_reference_proteins()

    print("reference proteins resolved:")
    for (virus, gene), seq in sorted(refs.items()):
        print(f"  {virus:<6} {gene:<6} {len(seq):>5} aa")
    missing = [(v, g) for v in REFSEQ for g in GENES if (v, g) not in refs]
    if missing:
        print(f"  NOT FOUND: {missing}")

    rows = list(csv.DictReader(open(LABELS, encoding="utf-8"), delimiter="\t"))
    subs = [r for r in rows if r["variant_type"] == "substitution"]

    results, stats = [], {}
    for r in subs:
        virus, gene = r["virus"], r["gene"]
        pos = int(r["position"])
        ref_aa = r["ref_aa"]
        seq = refs.get((virus, gene))
        if seq is None:
            verdict, observed = "no_reference", ""
        elif pos < 1 or pos > len(seq):
            verdict, observed = "out_of_range", ""
        else:
            observed = seq[pos - 1]
            verdict = "match" if observed == ref_aa else "mismatch"
        results.append({**r, "ref_protein_aa": observed, "verdict": verdict})
        k = f"{virus}:{gene}"
        stats.setdefault(k, {"match": 0, "mismatch": 0,
                             "out_of_range": 0, "no_reference": 0})
        stats[k][verdict] += 1

    checkable = [x for x in results
                 if x["verdict"] in ("match", "mismatch", "out_of_range")]
    n_match = sum(1 for x in checkable if x["verdict"] == "match")
    rate = n_match / len(checkable) if checkable else 0.0

    with open(OUT / "gate2_coordinate_check.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()),
                           delimiter="\t")
        w.writeheader()
        w.writerows(results)

    summary = {
        "threshold": GATE_2_THRESHOLD,
        "n_substitution_labels": len(subs),
        "n_checkable": len(checkable),
        "n_match": n_match,
        "agreement": round(rate, 4),
        "passed": rate >= GATE_2_THRESHOLD,
        "by_gene": stats,
        "reference_lengths": {f"{v}:{g}": len(s) for (v, g), s in refs.items()},
    }
    (REP / "gate2_summary.json").write_text(json.dumps(summary, indent=2,
                                                       sort_keys=True),
                                            encoding="utf-8")

    print("\nper gene:")
    for k, c in sorted(stats.items()):
        tot = sum(c.values())
        ok = c["match"] / tot if tot else 0
        print(f"  {k:<12} match {c['match']:>3}/{tot:<4} ({ok:6.1%})  "
              f"mismatch {c['mismatch']:>3}  out_of_range {c['out_of_range']:>3}")

    print(f"\nGATE 2  agreement {rate:.1%} (threshold {GATE_2_THRESHOLD:.0%}): "
          f"{'PASS' if summary['passed'] else 'FAIL'}")

    assert summary["passed"], (
        f"GATE 2 FAILED: {rate:.1%} agreement. Coordinates must be repaired "
        f"before any analysis proceeds. See {OUT / 'gate2_coordinate_check.tsv'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
