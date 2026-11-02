#!/usr/bin/env python
"""
Phase 2, step 1 - the reclassification table (Aim 1).

Joins every catalogued HSV resistance substitution to its frequency in
circulating virus, and flags the cases where the annotation and the population
frequency disagree. This is the project's primary deliverable.

Each labelled variant also gets a leave-one-study-out range, so a frequency
that rests entirely on one submission batch is visible as such rather than
being read as a population estimate.

Thresholds, declared in the protocol before the data was seen:
  resistance-relevant fold change   > 2x
  notable population frequency      > 1%
  common polymorphism               > 10%

Outputs:
  data/processed/reclassification_table.tsv
  reports/phase2/reclassification_summary.json
"""

import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase2"

FC_RESISTANT = 2.0
FREQ_NOTABLE = 0.01
FREQ_COMMON = 0.10
MIN_COUNT = 2          # singletons are not distinguishable from sequencing error


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def tsv(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def classify(freq, observed, count, fc):
    """Population-frequency band."""
    if not observed:
        return "not_observed"
    if count < MIN_COUNT:
        return "singleton"
    if freq > 0.5:
        return "majority_allele"
    if freq > FREQ_COMMON:
        return "common_polymorphism"
    if freq > FREQ_NOTABLE:
        return "uncommon_polymorphism"
    return "rare"


def disagreement(band, fc, tier):
    """Where the catalogue and the population disagree, and how."""
    strong = fc is not None and fc > FC_RESISTANT
    weak = fc is not None and fc <= FC_RESISTANT
    if band == "majority_allele":
        return ("reference_artifact" if not strong
                else "majority_allele_with_measured_resistance")
    if band in ("common_polymorphism", "uncommon_polymorphism"):
        if weak:
            return "reclassify_as_polymorphism"
        if strong:
            return "circulating_resistance_allele"
        return "unphenotyped_but_common"
    if band in ("rare", "singleton"):
        return "consistent_with_resistance" if strong else "consistent_rare"
    return "not_observed"


def main():
    REP.mkdir(parents=True, exist_ok=True)

    labels = [r for r in tsv(PROC / "hsv_mutations_tiered.tsv")
              if r["variant_type"] == "substitution"]
    freqs = {(r["virus"], r["gene"], int(r["position"]), r["alt_aa"]): r
             for r in tsv(PROC / "allele_frequencies.tsv")
             if r["stratum"] == "naive_proxy"}
    calls = tsv(PROC / "variant_calls.tsv")
    cover = tsv(PROC / "isolate_coverage.tsv")

    # Who carries each variant, and which studies contributed it.
    carriers = defaultdict(set)
    studies_of = defaultdict(set)
    for c in calls:
        k = (c["virus"], c["gene"], int(c["position"]), c["alt_aa"])
        carriers[k].add(c["accession"])
        # Same fallback as below. Counting only non-empty study keys understated
        # how concentrated a variant is: carriers with no study ID simply did
        # not appear in the tally, so "supported by 2 studies" could describe a
        # variant whose carriers came from one study plus unattributed records.
        studies_of[k].add(
            c["study_key"]
            or "ACC:" + re.match(r"^([A-Za-z]+)", c["accession"]).group(1))

    # Gene-level denominators per study, used for the leave-one-study-out range.
    # This is a gene-level approximation of the exact per-position coverage used
    # for the headline frequency; partial CDS make the two differ slightly.
    gene_isolates = defaultdict(set)
    study_of_acc = {}
    for r in cover:
        gene_isolates[(r["virus"], r["gene"])].add(r["accession"])
        # A record with no PubMed ID or BioProject previously got an empty
        # study key, and empty keys were excluded from the dropout loop - so
        # those isolates could never be dropped. A variant carried mostly by one
        # study plus a few study-less records then looked robust to dropout when
        # it was not. Falling back to the accession's letter prefix keeps a
        # submission block together and makes every isolate droppable.
        study_of_acc[r["accession"]] = (
            r["study_key"]
            or "ACC:" + re.match(r"^([A-Za-z]+)", r["accession"]).group(1))

    out, summary = [], defaultdict(int)
    for lab in labels:
        virus, gene = lab["virus"], lab["gene"]
        if not lab["position"]:
            continue
        pos = int(lab["position"])
        k = (virus, gene, pos, lab["alt_aa"])
        f = freqs.get(k)
        observed = f is not None
        count = int(f["count"]) if observed else 0
        covn = int(f["coverage"]) if observed else 0
        freq = float(f["frequency"]) if observed else 0.0
        lo = float(f["ci_low"]) if observed else 0.0
        hi = float(f["ci_high"]) if observed else 0.0
        fc = float(lab["max_fold_change"]) if lab["max_fold_change"] else None

        # Leave-one-study-out over the studies that contribute isolates to this
        # gene, not only those carrying the variant: dropping a large
        # variant-free study raises the frequency, which matters just as much.
        loso = []
        all_acc = gene_isolates.get((virus, gene), set())
        studies = {study_of_acc.get(a, "") for a in all_acc} - {""}
        for s in studies:
            keep = {a for a in all_acc if study_of_acc.get(a, "") != s}
            if len(keep) < 20:
                continue
            kcount = len(carriers[k] & keep)
            loso.append(kcount / len(keep))
        loso_min = round(min(loso), 6) if loso else None
        loso_max = round(max(loso), 6) if loso else None

        # Share of this variant's carriers contributed by its largest study.
        carrier_studies = Counter(study_of_acc.get(a, "")
                                  for a in carriers[k] & all_acc)
        max_share = (round(max(carrier_studies.values()) / sum(carrier_studies.values()), 4)
                     if carrier_studies else "")

        band = classify(freq, observed, count, fc)
        flag = disagreement(band, fc, lab["evidence_tier"])

        # A frequency that collapses below the notable threshold when one study
        # is dropped is not a population estimate.
        robust = (loso_min is not None and loso_min > FREQ_NOTABLE
                  if band in ("uncommon_polymorphism", "common_polymorphism",
                              "majority_allele") else None)

        out.append({
            "virus": virus, "gene": gene, "aa_change": lab["aa_change"],
            "position": pos, "ref_aa": lab["ref_aa"], "alt_aa": lab["alt_aa"],
            "evidence_tier": lab["evidence_tier"],
            "n_refs": lab["n_refs"],
            "drugs_tested": lab["drugs_tested"],
            "max_fold_change": lab["max_fold_change"],
            "observed": int(observed),
            "count": count, "coverage": covn,
            "frequency": round(freq, 6),
            "ci_low": round(lo, 6), "ci_high": round(hi, 6),
            "n_studies_carrying": len(studies_of[k]),
            "loso_min_freq": loso_min, "loso_max_freq": loso_max,
            "robust_to_study_dropout": "" if robust is None else int(robust),
            # robust_to_study_dropout is binary and hides how wide the dropout
            # range is. HSV-1 UL23 R222C passes it at 1.35% after a 3.4x fall
            # from 4.6%, with 8 of 11 carriers from one study - which is what an
            # ascertainment artefact looks like. These two columns expose that:
            # the largest single study's share of carriers, and whether dropping
            # one study more than halves the estimate.
            "max_study_share": max_share,
            "study_dominated": "" if not observed or loso_min is None else int(
                loso_min < 0.5 * freq),
            "frequency_band": band,
            "flag": flag,
        })
        summary[f"band:{band}"] += 1
        summary[f"flag:{flag}"] += 1
        summary[f"tier:{lab['evidence_tier']}"] += 1

    out.sort(key=lambda r: (-r["frequency"], r["virus"], r["gene"],
                            r["position"]))
    with open(PROC / "reclassification_table.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    (REP / "reclassification_summary.json").write_text(
        json.dumps({"n_labels": len(out), "counts": dict(summary),
                    "thresholds": {"fc_resistant": FC_RESISTANT,
                                   "freq_notable": FREQ_NOTABLE,
                                   "freq_common": FREQ_COMMON,
                                   "min_count": MIN_COUNT}},
                   indent=2, sort_keys=True), encoding="utf-8")

    print(f"reclassification table: {len(out)} labelled substitutions\n")
    print("frequency band:")
    for k in sorted(summary):
        if k.startswith("band:"):
            print(f"   {k[5:]:<26}{summary[k]}")
    print("\nflag:")
    for k in sorted(summary):
        if k.startswith("flag:"):
            print(f"   {k[5:]:<42}{summary[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
