#!/usr/bin/env python
"""
Phase 1, step 4 - per-position allele frequencies in reference coordinates.

Each isolate protein is aligned pairwise to the NCBI reference protein rather
than placed in a multiple-sequence alignment. Within a herpesvirus species the
sequences are ~98-99% identical, so a reference-anchored pairwise alignment
recovers the same columns an MSA would while mapping every observed residue
directly onto a reference coordinate. That removes the usual source of silent
coordinate drift, which matters here because the resistance labels are
themselves indexed against the reference.

A codon-level multiple alignment is still required for the selection analysis
in Phase 2; this step produces the amino-acid frequency table only.

Frequencies are stratified by a treatment-exposure proxy. Isolates whose record
mentions resistance or an antiviral are enriched for treatment failure, so they
are reported separately and excluded from the headline naive frequencies.

Outputs:
  data/processed/allele_frequencies.tsv
  reports/phase1/frequency_summary.json
"""

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parent.parent
GENE_DIR = ROOT / "data" / "processed" / "genes"
REF_DIR = ROOT / "data" / "raw" / "reference"
INTERIM = ROOT / "data" / "interim"
OUT = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase1"

# Words in a record that indicate the isolate was collected in a resistance or
# antiviral context. Matched against definition, note and isolation_source.
TREATMENT_TERMS = re.compile(
    r"resist|acyclovir|aciclovir|ganciclovir|foscarnet|cidofovir|penciclovir|"
    r"valacyclovir|famciclovir|brivudine|pritelivir|amenamevir|refractory|"
    r"antiviral|treatment.failure|acv[-_ ]?r\b|drug",
    re.I,
)


def wilson(k, n, z=1.96):
    """Wilson score interval - correct at the very low frequencies seen here,
    where the normal approximation gives negative lower bounds."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def read_fasta(path):
    name, buf, out = None, [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def make_aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    # Terminal gaps are free: many submissions are partial CDS and should not
    # be penalised for starting late or ending early.
    a.target_end_gap_score = 0.0
    a.query_end_gap_score = 0.0
    return a


def load_exposure():
    """accession -> 'treatment_associated' | 'naive_proxy'"""
    path = INTERIM / "isolate_metadata.tsv"
    if not path.exists():
        return {}
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        blob = " ".join(r.get(f, "") for f in
                        ("definition", "note", "isolation_source", "strain",
                         "isolate"))
        out[r["accession"]] = ("treatment_associated"
                               if TREATMENT_TERMS.search(blob)
                               else "naive_proxy")
    return out


def load_study_keys():
    """accession -> PubMed id or BioProject, for leave-one-study-out."""
    path = INTERIM / "isolate_metadata.tsv"
    if not path.exists():
        return {}
    return {r["accession"]: r.get("study_key", "")
            for r in csv.DictReader(open(path, encoding="utf-8"),
                                    delimiter="\t")}


def main():
    REP.mkdir(parents=True, exist_ok=True)
    aligner = make_aligner()
    exposure = load_exposure()
    study = load_study_keys()

    rows = []
    calls = []          # long format: one row per isolate per non-reference residue
    coverage_rows = []  # per isolate per gene, so a denominator can be recomputed
    summary = {}

    for faa in sorted(GENE_DIR.glob("*_*.faa")):
        virus, gene = faa.stem.split("_", 1)
        ref_path = REF_DIR / f"{virus}_{gene}.faa"
        if not ref_path.exists():
            continue
        ref = read_fasta(ref_path)[0][1]
        seqs = read_fasta(faa)

        # counts[pos][aa] = {stratum: n}; coverage[pos][stratum] = n
        counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        coverage = defaultdict(lambda: defaultdict(int))
        n_by_stratum = defaultdict(int)

        for acc, aa in seqs:
            stratum = exposure.get(acc, "naive_proxy")
            n_by_stratum[stratum] += 1
            aln = aligner.align(ref, aa)[0]
            covered = 0
            for (rs, re_), (qs, qe) in zip(*aln.aligned):
                for off in range(re_ - rs):
                    pos = rs + off + 1          # 1-based reference coordinate
                    obs = aa[qs + off]
                    if obs in "XB Z*":
                        continue
                    coverage[pos][stratum] += 1
                    counts[pos][obs][stratum] += 1
                    covered += 1
                    if obs != ref[pos - 1]:
                        calls.append({
                            "accession": acc, "virus": virus, "gene": gene,
                            "position": pos, "ref_aa": ref[pos - 1],
                            "alt_aa": obs, "stratum": stratum,
                            "study_key": study.get(acc, ""),
                        })
            coverage_rows.append({
                "accession": acc, "virus": virus, "gene": gene,
                "positions_covered": covered, "reference_length": len(ref),
                "study_key": study.get(acc, ""), "stratum": stratum,
            })

        n_total = len(seqs)
        for pos in sorted(counts):
            ref_aa = ref[pos - 1]
            for obs in sorted(counts[pos]):
                if obs == ref_aa:
                    continue                      # record variants only
                for stratum in ("naive_proxy", "treatment_associated"):
                    k = counts[pos][obs][stratum]
                    n = coverage[pos][stratum]
                    if n == 0:
                        continue
                    lo, hi = wilson(k, n)
                    rows.append({
                        "virus": virus, "gene": gene, "position": pos,
                        "ref_aa": ref_aa, "alt_aa": obs,
                        "stratum": stratum,
                        "count": k, "coverage": n,
                        "frequency": round(k / n, 6),
                        "ci_low": round(lo, 6), "ci_high": round(hi, 6),
                    })

        variable = sum(1 for pos in counts
                       if any(a != ref[pos - 1] for a in counts[pos]))
        summary[f"{virus}:{gene}"] = {
            "n_sequences": n_total,
            "n_by_stratum": dict(n_by_stratum),
            "reference_length": len(ref),
            "positions_covered": len(counts),
            "positions_with_variation": variable,
            "mean_coverage": round(
                sum(sum(c.values()) for c in coverage.values())
                / max(len(ref), 1), 1),
        }
        print(f"{virus}:{gene:<6} n={n_total:<5} "
              f"covered={len(counts):<5} variable={variable:<5} "
              f"naive={n_by_stratum['naive_proxy']} "
              f"treated={n_by_stratum['treatment_associated']}")

    if not rows:
        print("no frequency rows produced - has p1_03_extract.py run?")
        return 1

    for name, data in (("allele_frequencies.tsv", rows),
                       ("variant_calls.tsv", calls),
                       ("isolate_coverage.tsv", coverage_rows)):
        with open(OUT / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()),
                               delimiter="\t")
            w.writeheader()
            w.writerows(data)
    print(f"{len(calls)} per-isolate variant calls written")

    (REP / "frequency_summary.json").write_text(
        json.dumps({"per_gene": summary, "n_rows": len(rows)},
                   indent=2, sort_keys=True), encoding="utf-8")
    print(f"\n{len(rows)} frequency rows written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
