#!/usr/bin/env python
"""
Phase 1, step 3 - extract target coding sequences, QC them, de-duplicate.

Reads the batches fetched by p1_01_fetch.py and produces one clean per-isolate
record per gene, plus the isolate metadata table.

QC, applied per sequence:
  - length within 10% of the reference protein
  - no internal stop codon
  - ambiguous residue fraction below 2%

De-duplication, three levels as the protocol specifies:
  1. identical  - same amino-acid sequence across all target genes
  2. serial     - same strain or isolate name, which in HSV submissions is the
                  usual marker of repeat sampling from one patient
  3. study      - same BioProject, collapsed to distinct sequences only

Outputs:
  data/interim/isolate_metadata.tsv
  data/interim/cds_qc.tsv
  data/processed/genes/<virus>_<gene>.faa      QC-passed, de-duplicated
"""

import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_02_gate2 import GENES, match_gene, parse_cds_fasta, qualifiers  # noqa: E402

from Bio.Seq import Seq  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CDS_DIR = ROOT / "data" / "raw" / "cds"
META_DIR = ROOT / "data" / "raw" / "meta"
REF_DIR = ROOT / "data" / "raw" / "reference"
INTERIM = ROOT / "data" / "interim"
GENE_DIR = ROOT / "data" / "processed" / "genes"

LEN_TOLERANCE = 0.10
MAX_AMBIGUOUS = 0.02
# Within-species amino-acid divergence in these genes is well under 2%. A
# sequence beyond 5% from the reference is a recombinant, a misannotated CDS or
# a different virus, and two such records were found to dominate the HSV-2 UL23
# selection signal before this filter was added.
MAX_DIVERGENCE = 0.05

# Laboratory-adapted and reference strains. These are excluded from population
# frequency estimates because they are not independent natural isolates: they
# are passaged stocks shared between laboratories worldwide.
LAB_STRAINS = {
    "17", "17+", "17syn+", "strain 17", "kos", "f", "mp", "hfem", "cl101",
    "hg52", "333", "ms", "g", "sd90e", "angelotti", "186", "cl(101)",
}

META_FIELDS = ["strain", "isolate", "host", "country", "geo_loc_name",
               "collection_date", "isolation_source", "note", "lab_host"]


def parse_metadata():
    """Pull source-feature qualifiers and BioProject out of truncated GenBank."""
    out = {}
    for virus_dir in sorted(META_DIR.glob("*")):
        if not virus_dir.is_dir():
            continue
        virus = virus_dir.name
        for f in sorted(virus_dir.glob("*.txt")):
            for rec in f.read_text(encoding="utf-8").split("\n//\n"):
                m = re.search(r"^VERSION\s+(\S+)", rec, re.M)
                if not m:
                    continue
                acc = m.group(1)
                d = {"accession": acc, "virus": virus}
                dm = re.search(r"^DEFINITION\s+(.+?)(?=^\S)", rec, re.M | re.S)
                d["definition"] = re.sub(r"\s+", " ", dm.group(1)).strip() if dm else ""
                bp = re.search(r"BioProject:\s*(PRJ\w+)", rec)
                d["bioproject"] = bp.group(1) if bp else ""
                pm = re.search(r"^\s+PUBMED\s+(\d+)", rec, re.M)
                d["pubmed"] = pm.group(1) if pm else ""
                jm = re.search(r"^\s{2}JOURNAL\s+(.+?)(?=^\s{0,2}\S)",
                               rec, re.M | re.S)
                d["journal"] = re.sub(r"\s+", " ", jm.group(1)).strip()[:200] if jm else ""
                # A truncated record prints "source          1" with no range,
                # so the location is matched loosely rather than as "a..b".
                src = re.search(
                    r"^\s{5}source\s+\S+\s*\n(.*?)(?=^\s{5}\S|^[A-Z]{2,}|\Z)",
                    rec, re.M | re.S)
                block = src.group(1) if src else ""
                for q in META_FIELDS:
                    qm = re.search(rf'/{q}="([^"]*)"', block, re.S)
                    d[q] = re.sub(r"\s+", " ", qm.group(1)).strip() if qm else ""
                yr = re.search(r"(19|20)\d{2}", d.get("collection_date", ""))
                d["collection_year"] = yr.group(0) if yr else ""
                name = (d.get("strain") or d.get("isolate") or "").strip().lower()
                # Strain names are written inconsistently ("strain 17", "17syn+",
                # "HSV-1 KOS"); normalise before matching the lab-strain list.
                norm = re.sub(r"^(hsv[-\s]?[12]\s*|strain\s+)", "", name).strip()
                norm = re.sub(r"[\s_]+", " ", norm)
                d["strain_key"] = norm
                d["is_lab_strain"] = int(norm in LAB_STRAINS)
                # Study key for within-study de-duplication: the publication if
                # there is one, otherwise the submission's BioProject.
                d["study_key"] = d["pubmed"] or d["bioproject"]
                out[acc] = d
    return out


def load_references():
    refs = {}
    for f in REF_DIR.glob("*_*.faa"):
        stem = f.stem
        if stem.endswith("_cds_aa"):
            continue
        virus, gene = stem.split("_", 1)
        refs[(virus, gene)] = "".join(f.read_text(encoding="utf-8")
                                      .splitlines()[1:])
    return refs


def extract():
    refs = load_references()
    records = defaultdict(dict)          # (virus, acc) -> {gene: protein}
    # Nucleotide CDS for sequences that pass QC, needed for the codon-level
    # selection analysis in Phase 2.
    nucleotides = {}                     # (virus, acc, gene) -> nt
    qc_rows = []

    for virus_dir in sorted(CDS_DIR.glob("*")):
        if not virus_dir.is_dir():
            continue
        virus = virus_dir.name
        for f in sorted(virus_dir.glob("*.txt")):
            for header, nt in parse_cds_fasta(f.read_text(encoding="utf-8")):
                q = qualifiers(header)
                gene = match_gene(q)
                if not gene:
                    continue
                am = re.match(r">lcl\|([A-Za-z0-9_.]+)_cds", header)
                acc = am.group(1) if am else q.get("protein_id", "")
                if not acc:
                    continue
                ref = refs.get((virus, gene))
                if ref is None:
                    continue

                flags = []
                if q.get("pseudo") is not None or "pseudo" in q:
                    flags.append("pseudo")
                if len(nt) % 3:
                    flags.append("not_multiple_of_3")
                aa = str(Seq(nt[:len(nt) - len(nt) % 3]).translate())
                aa = aa[:-1] if aa.endswith("*") else aa
                if "*" in aa:
                    flags.append("internal_stop")
                amb = sum(1 for c in aa if c in "XB Z") / max(len(aa), 1)
                if amb > MAX_AMBIGUOUS:
                    flags.append("ambiguous")
                if abs(len(aa) - len(ref)) / len(ref) > LEN_TOLERANCE:
                    flags.append("length_outlier")
                # Ungapped divergence against the reference, which is adequate
                # to spot gross outliers without aligning at this stage.
                n_cmp = min(len(aa), len(ref))
                div = ((sum(1 for x, y in zip(ref, aa) if x != y)
                        + abs(len(aa) - len(ref))) / len(ref)) if n_cmp else 1.0
                if div > MAX_DIVERGENCE:
                    flags.append("divergence_outlier")

                qc_rows.append({
                    "virus": virus, "accession": acc, "gene": gene,
                    "nt_len": len(nt), "aa_len": len(aa),
                    "ref_aa_len": len(ref),
                    "ambiguous_frac": round(amb, 4),
                    "divergence": round(div, 4),
                    "flags": ";".join(flags) or "pass",
                    "passed": int(not flags),
                })
                if not flags:
                    records[(virus, acc)][gene] = aa
                    nucleotides[(virus, acc, gene)] = nt
    return records, qc_rows, nucleotides


def patient_groups(records, meta):
    """Infer which isolates come from the same patient, from strain names.

    Many studies name serial isolates <patient>-<visit>-<sample>: 10-6-1, 10-6-2,
    10-7-1 are one patient sampled repeatedly, not three independent infections.
    Collapsing identical sequences and identical strain names does not catch
    these, which inflates the denominator and, worse, lets one heavily sampled
    patient create an apparently polymorphic allele.

    Two deliberate restrictions, because over-collapsing would be a worse error
    than the one being fixed:

      numeric prefix only   10-6-1 groups on "10"; INV-25C and IV-7 do not
                            group, because an alphabetic prefix is usually a
                            site or collection code, not a patient
      not a year            2007-22031 and 2009-13458 are year-accession
                            identifiers, not patient-visit. Grouping on them
                            collapsed 29 independent HSV-2 isolates into 3,
                            which is a worse error than the under-collapse this
                            function exists to fix
      short tail segments   a patient-visit-sample id has small components
                            (10-6-1, 914-h2); a 5-digit tail like 2007-22031 is
                            an accession number. This is a second, independent
                            check on the same distinction, deliberately so
      scoped to a study     patient "10" in one study is unrelated to patient
                            "10" in another, so the key includes the study (or
                            the accession's letter prefix when no study is
                            recorded, which keeps a submission block together)

    A prefix seen only once is not a group and is left alone. Returns a map from
    (virus, accession) to a patient key, only for isolates that share one.
    """
    buckets = defaultdict(list)
    for (virus, acc) in sorted(records):
        m = meta.get(acc, {})
        strain = (m.get("strain_key") or "").strip()
        hit = re.match(r"^(\d+)-", strain)
        if not hit:
            continue
        prefix = hit.group(1)
        if len(prefix) == 4 and 1900 <= int(prefix) <= 2100:
            continue
        if any(len(seg) > 3 for seg in strain.split("-")[1:]):
            continue
        scope = m.get("study_key") or "ACC:" + re.match(r"^([A-Za-z]+)", acc).group(1)
        buckets[(virus, scope, prefix)].append((virus, acc))
    return {member: key
            for key, members in buckets.items() if len(members) > 1
            for member in members}


def deduplicate(records, meta):
    """Four-level collapse. Returns kept accessions and a per-level tally."""
    kept, reasons = {}, defaultdict(int)
    seen_seq, seen_strain, seen_study = {}, {}, defaultdict(set)
    patients = patient_groups(records, meta)
    seen_patient = {}
    collapsed = []

    # Deterministic order so the kept representative never depends on dict order.
    for (virus, acc) in sorted(records):
        genes = records[(virus, acc)]
        m = meta.get(acc, {})

        if m.get("is_lab_strain"):
            reasons["lab_strain"] += 1
            continue

        # Frequencies are claims about virus circulating in humans, so isolates
        # recorded from another host are excluded. A blank host is kept, since
        # most unannotated submissions are human.
        host = m.get("host", "").strip().lower()
        if host and "homo sapiens" not in host and "human" not in host:
            reasons["non_human_host"] += 1
            continue

        sig = hashlib.sha1(
            "|".join(f"{g}:{genes[g]}" for g in sorted(genes)).encode()
        ).hexdigest()
        if sig in seen_seq:
            reasons["identical"] += 1
            continue

        skey = m.get("strain_key", "")
        if skey and (virus, skey) in seen_strain:
            reasons["serial_same_strain"] += 1
            continue

        pkey = patients.get((virus, acc))
        if pkey is not None and pkey in seen_patient:
            reasons["same_patient_serial"] += 1
            collapsed.append((acc, skey, "|".join(map(str, pkey)),
                              seen_patient[pkey]))
            continue

        study = m.get("study_key", "")
        if study:
            if sig in seen_study[study]:
                reasons["study_near_duplicate"] += 1
                continue
            seen_study[study].add(sig)

        seen_seq[sig] = acc
        if skey:
            seen_strain[(virus, skey)] = acc
        if pkey is not None:
            seen_patient[pkey] = acc
        kept[(virus, acc)] = genes

    # Every patient collapse is written out, because the grouping is inferred
    # from a naming convention rather than stated by the submitter and has to
    # stay auditable.
    if collapsed:
        with open(INTERIM / "patient_collapse.tsv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["dropped_accession", "strain", "patient_key",
                        "kept_representative"])
            w.writerows(collapsed)
    return kept, dict(reasons)


def main():
    INTERIM.mkdir(parents=True, exist_ok=True)
    GENE_DIR.mkdir(parents=True, exist_ok=True)

    meta = parse_metadata()
    print(f"metadata records parsed: {len(meta)}")

    records, qc_rows, nucleotides = extract()
    print(f"isolates with at least one target gene: {len(records)}")

    with open(INTERIM / "cds_qc.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(qc_rows)

    if meta:
        cols = sorted({k for d in meta.values() for k in d})
        with open(INTERIM / "isolate_metadata.tsv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(meta.values())

    kept, reasons = deduplicate(records, meta)
    print(f"isolates kept after de-duplication: {len(kept)}")
    for k, v in sorted(reasons.items()):
        print(f"    removed {v:>4}  {k}")

    counts = defaultdict(int)
    by_gene = defaultdict(list)
    for (virus, acc), genes in sorted(kept.items()):
        for g, aa in genes.items():
            by_gene[(virus, g)].append((acc, aa))
            counts[f"{virus}:{g}"] += 1

    for (virus, gene), seqs in sorted(by_gene.items()):
        with open(GENE_DIR / f"{virus}_{gene}.faa", "w", encoding="utf-8") as fh:
            for acc, aa in seqs:
                fh.write(f">{acc}\n{aa}\n")
        with open(GENE_DIR / f"{virus}_{gene}.fna", "w", encoding="utf-8") as fh:
            for acc, _ in seqs:
                nt = nucleotides.get((virus, acc, gene))
                if nt:
                    fh.write(f">{acc}\n{nt}\n")

    print("\nsequences per gene after QC and de-duplication:")
    for k, v in sorted(counts.items()):
        print(f"  {k:<12} {v}")

    qc_fail = sum(1 for r in qc_rows if not r["passed"])
    (ROOT / "reports" / "phase1").mkdir(parents=True, exist_ok=True)
    (ROOT / "reports" / "phase1" / "extract_summary.json").write_text(
        json.dumps({
            "metadata_records": len(meta),
            "cds_extracted": len(qc_rows),
            "cds_qc_failed": qc_fail,
            "isolates_with_targets": len(records),
            "isolates_kept": len(kept),
            "removed_by": reasons,
            "sequences_per_gene": dict(counts),
        }, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
