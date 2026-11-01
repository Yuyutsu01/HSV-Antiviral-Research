#!/usr/bin/env python
"""
Phase 2, step 3 - pre-therapeutic helicase-primase baseline (Aim 3).

Pritelivir and amenamevir target the UL5/UL52 helicase-primase complex. Neither
is in wide clinical use, so the natural variation present now is a baseline that
cannot be recovered once selection pressure is applied. This step records it.

The headline quantity is the proportion of circulating isolates already carrying
a substitution catalogued as helicase-primase resistance-associated, together
with the overall level of natural variation across both genes.

H3 states that catalogued resistance positions should be more conserved than the
gene as a whole, because they sit in functionally constrained regions. That is
tested here by comparing variation at catalogued positions against the rest of
the gene, with a permutation test over position sets.

Outputs:
  data/processed/hp_baseline.tsv
  reports/phase2/hp_baseline_summary.json
"""

import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase2"

HP_GENES = ("UL5", "UL52")
N_PERMUTATIONS = 10000
SEED = 20261102


def tsv(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def main():
    REP.mkdir(parents=True, exist_ok=True)

    labels = [r for r in tsv(PROC / "hsv_mutations_tiered.tsv")
              if r["gene"] in HP_GENES and r["variant_type"] == "substitution"]
    calls = [r for r in tsv(PROC / "variant_calls.tsv") if r["gene"] in HP_GENES]
    cover = [r for r in tsv(PROC / "isolate_coverage.tsv")
             if r["gene"] in HP_GENES]
    freqs = {(r["virus"], r["gene"], int(r["position"]), r["alt_aa"]): r
             for r in tsv(PROC / "allele_frequencies.tsv")
             if r["gene"] in HP_GENES and r["stratum"] == "naive_proxy"}

    isolates = defaultdict(set)
    for r in cover:
        isolates[(r["virus"], r["gene"])].add(r["accession"])

    # Which isolates carry any catalogued helicase-primase substitution
    labelled_keys = {(r["virus"], r["gene"], int(r["position"]), r["alt_aa"])
                     for r in labels if r["position"]}
    carriers = defaultdict(set)
    for c in calls:
        k = (c["virus"], c["gene"], int(c["position"]), c["alt_aa"])
        if k in labelled_keys:
            carriers[(c["virus"], c["gene"])].add(c["accession"])

    rows = []
    for r in sorted(labels, key=lambda x: (x["virus"], x["gene"],
                                           int(x["position"] or 0))):
        if not r["position"]:
            continue
        k = (r["virus"], r["gene"], int(r["position"]), r["alt_aa"])
        f = freqs.get(k)
        rows.append({
            "virus": r["virus"], "gene": r["gene"],
            "aa_change": r["aa_change"], "position": int(r["position"]),
            "evidence_tier": r["evidence_tier"],
            "drugs_tested": r["drugs_tested"],
            "max_fold_change": r["max_fold_change"],
            "observed": int(f is not None),
            "count": f["count"] if f else 0,
            "coverage": f["coverage"] if f else 0,
            "frequency": f["frequency"] if f else 0,
            "ci_high": f["ci_high"] if f else 0,
        })

    # Variation at catalogued positions against the rest of the gene.
    varpos = defaultdict(lambda: defaultdict(int))
    for c in calls:
        varpos[(c["virus"], c["gene"])][int(c["position"])] += 1

    constraint = {}
    rng = random.Random(SEED)
    for (virus, gene), counts in sorted(varpos.items()):
        n_iso = len(isolates[(virus, gene)])
        cat = sorted({int(r["position"]) for r in labels
                      if r["virus"] == virus and r["gene"] == gene
                      and r["position"]})
        if not cat or n_iso == 0:
            continue
        gene_len = max(max(counts), max(cat))
        observed = sum(counts.get(p, 0) for p in cat) / (len(cat) * n_iso)
        overall = sum(counts.values()) / (gene_len * n_iso)

        # Permutation: draw position sets of the same size at random.
        all_pos = list(range(1, gene_len + 1))
        hits = 0
        for _ in range(N_PERMUTATIONS):
            draw = rng.sample(all_pos, len(cat))
            v = sum(counts.get(p, 0) for p in draw) / (len(cat) * n_iso)
            if v <= observed:
                hits += 1
        constraint[f"{virus}:{gene}"] = {
            "n_catalogued_positions": len(cat),
            "n_isolates": n_iso,
            "variation_at_catalogued": round(observed, 6),
            "variation_gene_wide": round(overall, 6),
            "ratio": round(observed / overall, 3) if overall else None,
            "p_more_conserved_than_random": round(1 - hits / N_PERMUTATIONS, 4),
        }

    baseline = {}
    for (virus, gene), accs in sorted(isolates.items()):
        c = carriers.get((virus, gene), set())
        baseline[f"{virus}:{gene}"] = {
            "n_isolates": len(accs),
            "n_carrying_catalogued_substitution": len(c),
            "proportion": round(len(c) / len(accs), 5) if accs else 0,
        }

    with open(PROC / "hp_baseline.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    (REP / "hp_baseline_summary.json").write_text(
        json.dumps({"baseline_carriage": baseline,
                    "constraint_test": constraint,
                    "n_permutations": N_PERMUTATIONS, "seed": SEED},
                   indent=2, sort_keys=True), encoding="utf-8")

    print("pre-therapeutic carriage of catalogued helicase-primase substitutions:")
    for k, v in sorted(baseline.items()):
        print(f"  {k:<12} {v['n_carrying_catalogued_substitution']:>3}"
              f"/{v['n_isolates']:<4} ({v['proportion']:.2%})")
    print("\nH3, are catalogued positions more conserved than the gene overall:")
    for k, v in sorted(constraint.items()):
        print(f"  {k:<12} catalogued {v['variation_at_catalogued']:.5f} vs "
              f"gene-wide {v['variation_gene_wide']:.5f}  "
              f"ratio {v['ratio']}  p={v['p_more_conserved_than_random']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
