#!/usr/bin/env python
"""
Phase 3, step 3 - assemble the modelling dataset.

Target: log10 of the maximum measured fold change in EC50. Regression is the
primary task because fold change is a quantitative measurement and collapsing
it to resistant/susceptible throws away most of what the experiments actually
report. The binary view at a 2x threshold is kept as a secondary column.

Rows are restricted to single substitutions with an engineered genotype
(Tier A) and no reported co-mutation, so each number is attributable to one
amino-acid change. Tier B (lab-selected) rows are assembled too and marked, for
a sensitivity analysis with more rows and noisier labels.

Feature groups, kept separable so each can be ablated:
  freq_*        population frequency from Aim 1
  cons_*        per-site variability from Aim 2
  prop_*        amino-acid property deltas and BLOSUM62
  struct_*      relative solvent accessibility, distance to bound ligand
  esm_*         ESM-2 zero-shot masked marginal

The gene a variant sits in predicts the label prevalence on its own, so gene is
recorded as a grouping key, not offered as a feature. Grouped cross-validation
in p3_04 splits on gene and on position blocks so neighbouring residues cannot
leak across folds.

Outputs:
  data/processed/model_dataset.tsv
  reports/phase3/dataset_summary.json
"""

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REF_DIR = ROOT / "data" / "raw" / "reference"
REP = ROOT / "reports" / "phase3"

FC_THRESHOLD = 2.0
# Tier C rows carry clinical-isolate phenotypes, so a measured fold change may
# reflect an unreported co-mutation rather than the substitution named. They are
# assembled and marked rather than used in the primary analysis, because only
# 11% of Tier A variants are seen in circulating virus at all, which leaves the
# frequency baseline almost no variance to work with on Tier A alone.
TIERS = {"A_engineered": "A", "B_lab_selected": "B",
         "C_clinical_isolate": "C"}

# Kyte-Doolittle hydropathy.
HYDRO = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
# Residue volume, cubic angstroms (Zamyatnin 1972).
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5}
POLAR = set("STNQCYHKRDE")


def load_tsv(path):
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def read_fasta_one(path):
    return "".join(l.strip() for l in
                   path.read_text(encoding="utf-8").splitlines()[1:] if l.strip())


def main():
    REP.mkdir(parents=True, exist_ok=True)
    blosum = substitution_matrices.load("BLOSUM62")

    muts = load_tsv(PROC / "hsv_mutations_tiered.tsv")
    freq = {(r["virus"], r["gene"], r["position"], r["alt_aa"]): r
            for r in load_tsv(PROC / "allele_frequencies.tsv")
            if r["stratum"] == "naive_proxy"}
    site = {(r["virus"], r["gene"], r["position"]): r
            for r in load_tsv(PROC / "selection_by_site.tsv")}
    struct = {(r["virus"], r["gene"], r["position"]): r
              for r in load_tsv(PROC / "structure_features.tsv")}
    esm_rows = load_tsv(PROC / "esm2_scores.tsv")
    if not esm_rows:
        print("esm2_scores.tsv missing or empty - run p3_01_esm.py first")
        return 1
    esm = {(r["virus"], r["gene"], r["position"], r["alt_aa"]): r
           for r in esm_rows}

    lengths = {}
    for p in REF_DIR.glob("*_*.faa"):
        if p.stem.endswith("_cds_aa"):
            continue
        v, g = p.stem.split("_", 1)
        lengths[(v, g)] = len(read_fasta_one(p))

    rows, dropped = [], Counter()
    for m in muts:
        if m["variant_type"] != "substitution":
            dropped["not_substitution"] += 1
            continue
        tier = TIERS.get(m["evidence_tier"])
        if tier is None:
            dropped["tier_C_clinical"] += 1
            continue
        fc = fnum(m["max_fold_change"])
        if fc is None or fc <= 0:
            dropped["no_numeric_fold_change"] += 1
            continue
        if m["co_mutation_flag"] != "0":
            dropped["co_mutation"] += 1
            continue
        key3 = (m["virus"], m["gene"], m["position"])
        if key3[:2] not in lengths:
            dropped["no_reference_protein"] += 1
            continue
        key4 = (m["virus"], m["gene"], m["position"], m["alt_aa"])
        e = esm.get(key4)
        if e is None:
            dropped["no_esm_score"] += 1
            continue

        ref_aa, alt_aa = m["ref_aa"], m["alt_aa"]
        f = freq.get(key4)
        s = site.get(key3)
        st = struct.get(key3)

        cov = fnum(s["coverage"], 0) if s else 0
        nonsyn = fnum(s["n_nonsyn"], 0) if s else 0

        row = {
            "virus": m["virus"], "gene": m["gene"], "position": int(m["position"]),
            "ref_aa": ref_aa, "alt_aa": alt_aa, "aa_change": m["aa_change"],
            "tier": tier,
            "drugs_tested": m["drugs_tested"], "n_refs": m["n_refs"],
            # targets
            "fold_change": fc,
            "log10_fc": round(math.log10(fc), 5),
            "resistant": int(fc > FC_THRESHOLD),
            # Aim 1 features
            "freq_frequency": round(fnum(f["frequency"], 0.0), 6) if f else 0.0,
            "freq_count": int(fnum(f["count"], 0)) if f else 0,
            "freq_observed": int(f is not None),
            "freq_ci_high": round(fnum(f["ci_high"], 0.0), 6) if f else 0.0,
            # Aim 2 features
            "cons_site_variability": round(nonsyn / cov, 6) if cov else "",
            "cons_site_coverage": int(cov),
            # property features
            "prop_blosum62": int(blosum[ref_aa, alt_aa])
            if (ref_aa in blosum.alphabet and alt_aa in blosum.alphabet) else "",
            "prop_d_hydropathy": round(HYDRO[alt_aa] - HYDRO[ref_aa], 3)
            if ref_aa in HYDRO and alt_aa in HYDRO else "",
            "prop_d_volume": round(VOLUME[alt_aa] - VOLUME[ref_aa], 2)
            if ref_aa in VOLUME and alt_aa in VOLUME else "",
            "prop_d_charge": round(CHARGE.get(alt_aa, 0) - CHARGE.get(ref_aa, 0), 2),
            "prop_polarity_switch": int((ref_aa in POLAR) != (alt_aa in POLAR)),
            "prop_rel_position": round(int(m["position"]) / lengths[key3[:2]], 5),
            "prop_to_proline": int(alt_aa == "P"),
            "prop_from_glycine": int(ref_aa == "G"),
            # structural features
            "struct_rsa": st["rsa"] if st else "",
            "struct_dist_to_ligand": st["dist_to_ligand"] if st else "",
            "struct_site_resolved": st["site_resolved"] if st else "",
            "struct_inherited": st["inherited"] if st else "",
            "struct_covered": int(st is not None),
            # ESM-2
            "esm_masked_marginal": e["esm2_masked_marginal"],
            "esm_ref_logprob": e["esm2_ref_logprob"],
        }
        rows.append(row)

    rows.sort(key=lambda r: (r["virus"], r["gene"], r["position"], r["alt_aa"]))
    with open(PROC / "model_dataset.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def describe(sel):
        fcs = sorted(r["fold_change"] for r in sel)
        return {
            "n": len(sel),
            "n_resistant": sum(r["resistant"] for r in sel),
            "n_negative": sum(1 - r["resistant"] for r in sel),
            "fc_min": fcs[0] if fcs else None,
            "fc_median": fcs[len(fcs) // 2] if fcs else None,
            "fc_max": fcs[-1] if fcs else None,
            "by_gene": dict(sorted(Counter(
                f"{r['virus']}:{r['gene']}" for r in sel).items())),
            "struct_covered": sum(r["struct_covered"] for r in sel),
            "freq_observed": sum(r["freq_observed"] for r in sel),
        }

    tier_a = [r for r in rows if r["tier"] == "A"]
    summary = {
        "fc_threshold": FC_THRESHOLD,
        "dropped": dict(dropped),
        "primary_tier_A": describe(tier_a),
        "sensitivity_tier_AB": describe([r for r in rows if r["tier"] != "C"]),
        "sensitivity_all_tiers": describe(rows),
    }
    (REP / "dataset_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"{len(rows)} rows written ({len(tier_a)} Tier A primary)")
    print("\ndropped:")
    for k, v in dropped.most_common():
        print(f"  {k:28s} {v}")
    a = summary["primary_tier_A"]
    print(f"\nTier A: n={a['n']}  resistant={a['n_resistant']}  "
          f"negative={a['n_negative']}")
    print(f"  fold change {a['fc_min']} .. {a['fc_max']} (median {a['fc_median']})")
    print(f"  structure covered {a['struct_covered']}/{a['n']}, "
          f"observed in population {a['freq_observed']}/{a['n']}")
    print("  by gene:", a["by_gene"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
