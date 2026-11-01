#!/usr/bin/env python
"""
Phase 0 audit of the HerpesDRG database: establish the HSV label budget.

Reproduces every number in reports/phase0/PHASE0_REPORT.md from the pinned
raw file in data/raw/. Run from the project root:

    python scripts/phase0_audit.py

Gates are expressed as assertions. If a gate fails the script exits non-zero
and the project plan changes; it does not warn and continue.
"""

import csv
import collections
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "herpesdrg-db.tsv"
OUT = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase0"

# Gate thresholds, declared before the data was examined.
GATE_TIER_A_HSV1 = 50

# README: "it is important you filter for rows where the column status == 'A'"
ACTIVE = "A"

HSV = ("HSV1", "HSV2")

# strain_generation -> evidence tier. Values are lowercased before lookup
# because the raw column mixes "Marker Transfer" and "Marker transfer".
TIER_MAP = {
    "recombinant bac": "A_engineered",
    "marker transfer": "A_engineered",
    "lab induced strain": "B_lab_selected",
    "lab derived mutant": "B_lab_selected",
    "drug induced lab mutant": "B_lab_selected",
    "isolated strain": "C_clinical_isolate",
}
TIER_ORDER = ["A_engineered", "B_lab_selected", "C_clinical_isolate", "Z_unknown"]

DRUGS = [
    "Ganciclovir", "Aciclovir", "Cidofovir", "Foscarnet", "Brincidofovir",
    "Letermovir", "Brivudine", "Penciclovir", "Tomeglovir", "Maribavir",
    "Cyclopropavir", "Amenamevir", "Pritelivir",
]

RE_SUB = re.compile(r"^([A-Z])(\d+)([A-Z])$")
RE_STOP = re.compile(r"^([A-Z])(\d+)\*$")
RE_INDEL = re.compile(r"fs|del|ins|dup|frameshift|_", re.I)


def tier_of(strain_generation: str) -> str:
    return TIER_MAP.get(strain_generation.strip().lower(), "Z_unknown")


def variant_type(aa_change: str) -> str:
    a = aa_change.strip()
    if RE_SUB.match(a):
        return "substitution"
    if RE_STOP.match(a):
        return "nonsense"
    if RE_INDEL.search(a):
        return "indel_or_frameshift"
    if a.isdigit():
        return "bare_position"
    return "other"


def parse_substitution(aa_change: str):
    """Return (ref_aa, position, alt_aa) or None for non-substitutions."""
    m = RE_SUB.match(aa_change.strip())
    return (m.group(1), int(m.group(2)), m.group(3)) if m else None


def load_active():
    with open(RAW, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    active = [r for r in rows if r["status"].strip() == ACTIVE]
    return rows, active


def main():
    all_rows, active = load_active()
    hsv_rows = [r for r in active if r["virus"].strip() in HSV]

    # Collapse rows to unique mutations. A mutation appears once per
    # drug/reference/assay combination, so row counts overstate label counts.
    muts = {}
    for r in hsv_rows:
        key = (r["virus"].strip(), r["gene"].strip(), r["aa_change"].strip())
        rec = muts.setdefault(key, {
            "tiers": set(), "n_rows": 0, "refs": set(),
            "drugs": set(), "fold_changes": {}, "co_mut": False,
        })
        rec["tiers"].add(tier_of(r["strain_generation"]))
        rec["n_rows"] += 1
        if r["ref_doi"].strip():
            rec["refs"].add(r["ref_doi"].strip())
        if r["co_gene"].strip() or r["co_aa"].strip():
            rec["co_mut"] = True
        for d in DRUGS:
            v = r.get(d, "").strip()
            if v:
                rec["drugs"].add(d)
                try:
                    rec["fold_changes"].setdefault(d, []).append(float(v))
                except ValueError:
                    pass

    for rec in muts.values():
        rec["best_tier"] = next(t for t in TIER_ORDER if t in rec["tiers"])

    # ---------------------------------------------------------------- summary
    summary = {
        "rows_total": len(all_rows),
        "rows_active": len(active),
        "rows_hsv_active": len(hsv_rows),
        "unique_hsv_mutations": len(muts),
        "status_counts": dict(collections.Counter(r["status"].strip() for r in all_rows)),
        "virus_counts": dict(collections.Counter(r["virus"].strip() for r in active)),
        "tier_counts": {},
        "tier_by_virus": {},
        "tier_by_gene": {},
        "variant_types": dict(collections.Counter(variant_type(k[2]) for k in muts)),
        "n_independent_refs": {},
    }

    tier_all = collections.Counter()
    tier_vir = collections.defaultdict(collections.Counter)
    tier_gene = collections.defaultdict(collections.Counter)
    for (vir, gene, aa), rec in muts.items():
        t = rec["best_tier"]
        tier_all[t] += 1
        tier_vir[vir][t] += 1
        tier_gene[f"{vir}:{gene}"][t] += 1

    summary["tier_counts"] = dict(tier_all)
    summary["tier_by_virus"] = {v: dict(c) for v, c in tier_vir.items()}
    summary["tier_by_gene"] = {g: dict(c) for g, c in tier_gene.items()}

    # Independent-reference count per tier: the real unit of evidence, since
    # one paper can contribute many mutations.
    for t in TIER_ORDER:
        refs = set()
        for rec in muts.values():
            if rec["best_tier"] == t:
                refs |= rec["refs"]
        summary["n_independent_refs"][t] = len(refs)

    # ------------------------------------------------------------------ gates
    tier_a_hsv1 = tier_vir["HSV1"]["A_engineered"]
    tier_a_both = tier_all["A_engineered"]
    summary["gate_1"] = {
        "threshold": GATE_TIER_A_HSV1,
        "hsv1_tier_a": tier_a_hsv1,
        "hsv1_hsv2_tier_a": tier_a_both,
        "passed": tier_a_hsv1 >= GATE_TIER_A_HSV1,
    }

    # Substitutions are the modelling scope; everything else is a separate stream.
    subs = {k: v for k, v in muts.items() if variant_type(k[2]) == "substitution"}
    subs_tier_a_hsv1 = sum(
        1 for (vir, _, _), rec in subs.items()
        if vir == "HSV1" and rec["best_tier"] == "A_engineered"
    )
    summary["gate_1_substitutions_only"] = {
        "hsv1_tier_a_substitutions": subs_tier_a_hsv1,
        "passed": subs_tier_a_hsv1 >= GATE_TIER_A_HSV1,
    }

    # -------------------------------------------------------------- write out
    OUT.mkdir(parents=True, exist_ok=True)
    REP.mkdir(parents=True, exist_ok=True)

    with open(OUT / "hsv_mutations_tiered.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "virus", "gene", "aa_change", "variant_type", "ref_aa", "position",
            "alt_aa", "evidence_tier", "n_rows", "n_refs", "co_mutation_flag",
            "drugs_tested", "max_fold_change",
        ])
        for (vir, gene, aa), rec in sorted(muts.items()):
            p = parse_substitution(aa)
            fcs = [x for vals in rec["fold_changes"].values() for x in vals]
            w.writerow([
                vir, gene, aa, variant_type(aa),
                p[0] if p else "", p[1] if p else "", p[2] if p else "",
                rec["best_tier"], rec["n_rows"], len(rec["refs"]),
                int(rec["co_mut"]),
                ";".join(sorted(rec["drugs"])),
                f"{max(fcs):.3g}" if fcs else "",
            ])

    with open(REP / "phase0_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    # ------------------------------------------------------------------ print
    print(f"rows total/active/HSV : {len(all_rows)} / {len(active)} / {len(hsv_rows)}")
    print(f"unique HSV mutations  : {len(muts)}")
    print()
    for t in TIER_ORDER:
        if tier_all[t]:
            print(f"  {t:<22} {tier_all[t]:>4}   "
                  f"HSV1={tier_vir['HSV1'][t]:<4} HSV2={tier_vir['HSV2'][t]:<4} "
                  f"refs={summary['n_independent_refs'][t]}")
    print()
    print(f"GATE 1  Tier A HSV-1 >= {GATE_TIER_A_HSV1}: "
          f"{tier_a_hsv1} -> {'PASS' if summary['gate_1']['passed'] else 'FAIL'}")
    print(f"        substitutions only: {subs_tier_a_hsv1} -> "
          f"{'PASS' if summary['gate_1_substitutions_only']['passed'] else 'FAIL'}")

    assert summary["gate_1"]["passed"], (
        f"GATE 1 FAILED: only {tier_a_hsv1} Tier A HSV-1 mutations "
        f"(need {GATE_TIER_A_HSV1}). The supervised arm (Aim 4) is dropped."
    )
    assert summary["gate_1_substitutions_only"]["passed"], (
        f"GATE 1 FAILED on substitutions: {subs_tier_a_hsv1} "
        f"(need {GATE_TIER_A_HSV1})."
    )
    print("\nAll Phase 0 gates passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
