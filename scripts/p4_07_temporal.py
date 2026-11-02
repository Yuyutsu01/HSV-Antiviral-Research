"""Phase 4, step 7 — the temporal hold-out.

The protocol: "Train on mutations characterised before a cut-off year and test
on those characterised after it. This approximates prospective validation --
can the model, knowing only what was published by a given date, identify what
was subsequently discovered? It is cheap, and it is the single most persuasive
validation in the plan."

It was blocked until now because publication year was recoverable from the DOI
string for only 35 of 89 studies, clustered by publisher. Crossref resolved all
89 (`p4_06_doi_years.py`), so the split is no longer confounded with journal.

A mutation's year is the **earliest** year among its references: when the field
first characterised it. Using the latest would ask the model to predict things
already known.

**The cut-off year was not pre-declared.** Choosing one on test performance
would undo the point of the exercise, so every cut-off from 2005 to 2015 is
reported, and the primary is 2012 on the outcome-independent grounds that it
splits the primary set most evenly (59 train / 60 test). Read the sensitivity
table, not the primary line, if the two disagree.

Hyperparameters are the fixed declared values. There is no inner tuning loop:
the split is a single hold-out, and tuning against it would be selecting on the
test set.
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p3_04_baselines import (                                   # noqa: E402
    FEATURE_SETS, SEED, STAGES, load_dataset, make_models, matrix,
)
from p4_01_metrics import (                                     # noqa: E402
    average_precision, confusion_at_2x, precision_at_recall,
)
from p4_06_doi_years import normalise                           # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "reports" / "phase4"
N_BOOTSTRAP = 2000

PRIMARY_CUTOFF = 2012
CUTOFFS = [2005, 2008, 2010, 2012, 2015]

MODELS = [
    ("frequency_only", FEATURE_SETS["frequency_only"], "ridge",
     {"m__alpha": 1.0}),
    ("conservation_only", FEATURE_SETS["conservation_only"], "ridge",
     {"m__alpha": 1.0}),
    ("structure_only", FEATURE_SETS["structure_only"], "ridge",
     {"m__alpha": 1.0}),
    ("structure_only_forest", FEATURE_SETS["structure_only"], "forest",
     {"m__min_samples_leaf": 5, "m__max_features": 0.5, "m__n_estimators": 300}),
    ("esm_only", FEATURE_SETS["esm_only"], "ridge", {"m__alpha": 1.0}),
    ("full_model", STAGES["stage5_plus_esm"], "ridge", {"m__alpha": 1.0}),
    ("full_model_forest", STAGES["stage5_plus_esm"], "forest",
     {"m__min_samples_leaf": 5, "m__max_features": 0.5, "m__n_estimators": 300}),
]

LABEL_SETS = {"primary_tier_A": ("A",), "secondary_all_tiers": ("A", "B", "C")}


def first_year_by_mutation():
    with open(ROOT / "data" / "raw" / "doi_years.json", encoding="utf-8") as fh:
        cache = json.load(fh)
    years = {d: v["year"] for d, v in cache["years"].items() if v.get("year")}
    with open(ROOT / "data" / "raw" / "herpesdrg-db.tsv", encoding="utf-8") as fh:
        db = [r for r in csv.DictReader(fh, delimiter="\t")
              if r["status"] == "A" and r["virus"] in ("HSV1", "HSV2")]
    seen = defaultdict(set)
    for r in db:
        y = years.get(normalise(r["ref_doi"]))
        if y:
            seen[(r["virus"], r["gene"], r["aa_change"])].add(y)
    return {k: min(v) for k, v in seen.items()}, cache.get("retrieved")


def score(y_bin, y_cont, pred, rng):
    conf = confusion_at_2x(y_bin, pred)
    ap = average_precision(y_bin, pred)
    boot = []
    idx = np.arange(len(y_bin))
    for _ in range(N_BOOTSTRAP):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(set(y_bin[b])) < 2:
            continue
        boot.append(average_precision(y_bin[b], pred[b]))
    out = {
        "n_test": int(len(y_bin)),
        "prevalence": round(float(y_bin.mean()), 4),
        "average_precision": round(float(ap), 4),
        "ap_lift_over_prevalence": round(float(ap - y_bin.mean()), 4),
        "precision_at_95_recall": round(
            float(precision_at_recall(y_bin, pred)), 4),
        "spearman": round(float(spearmanr(y_cont, pred).statistic), 4),
        "threshold_degenerate": bool(
            conf["n_predicted_resistant"] in (0, len(y_bin))),
    }
    if boot:
        out["ap_ci95"] = [round(float(v), 4)
                          for v in np.percentile(boot, [2.5, 97.5])]
    return out


def split_run(rows, y, binary, year_of, cutoff, rng):
    yrs = np.array([year_of.get((r["virus"], r["gene"], r["aa_change"]), -1)
                    for r in rows])
    tr = np.where((yrs > 0) & (yrs <= cutoff))[0]
    te = np.where(yrs > cutoff)[0]
    res = {"cutoff": cutoff, "n_train": int(len(tr)), "n_test": int(len(te)),
           "n_without_year": int(np.sum(yrs <= 0))}
    if len(tr) < 20 or len(te) < 15:
        res["note"] = "split too small"
        return res
    if len(set(binary[te])) < 2 or len(set(binary[tr])) < 2:
        res["note"] = "a side is single-class"
        return res
    res["train_prevalence"] = round(float(binary[tr].mean()), 4)
    res["test_prevalence"] = round(float(binary[te].mean()), 4)
    res["models"] = {}
    for name, feats, est_name, params in MODELS:
        X = matrix(rows, feats)
        est, _ = make_models()[est_name]
        est = est.set_params(**params)
        est.fit(X[tr], y[tr])
        pred = est.predict(X[te])
        if np.std(pred) < 1e-12:
            res["models"][name] = {"note": "constant prediction"}
            continue
        res["models"][name] = score(binary[te], y[te], pred, rng)
    ranked = [(k, v["average_precision"]) for k, v in res["models"].items()
              if "average_precision" in v]
    res["best_by_average_precision"] = (max(ranked, key=lambda kv: kv[1])[0]
                                        if ranked else None)
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    year_of, retrieved = first_year_by_mutation()
    print(f"DOI years retrieved {retrieved}; "
          f"{len(year_of)} mutations carry a first-characterisation year\n",
          flush=True)

    summary = {"seed": SEED, "n_bootstrap": N_BOOTSTRAP,
               "primary_metric": "average_precision",
               "primary_cutoff": PRIMARY_CUTOFF,
               "cutoff_chosen_by": "most even split of the primary set; "
                                   "not pre-declared, not chosen on performance",
               "mutation_year_rule": "earliest year among the mutation's references",
               "doi_years_retrieved": retrieved,
               "label_sets": {}}

    for label, tiers in LABEL_SETS.items():
        rows = load_dataset(tiers)
        y = np.array([float(r["log10_fc"]) for r in rows])
        binary = np.array([int(r["resistant"]) for r in rows])
        covered = sum(1 for r in rows
                      if (r["virus"], r["gene"], r["aa_change"]) in year_of)
        print(f"################ {label}  n={len(rows)}  "
              f"with a year: {covered}", flush=True)

        entry = {"n": len(rows), "n_with_year": covered, "cutoffs": {}}
        for cutoff in CUTOFFS:
            r = split_run(rows, y, binary, year_of, cutoff, rng)
            entry["cutoffs"][str(cutoff)] = r
            tag = "  <-- primary" if cutoff == PRIMARY_CUTOFF else ""
            if "note" in r:
                print(f"\n=== cutoff {cutoff}: {r['note']}{tag}", flush=True)
                continue
            print(f"\n=== cutoff {cutoff}: train {r['n_train']} "
                  f"(prev {r['train_prevalence']:.2f}) -> test {r['n_test']} "
                  f"(prev {r['test_prevalence']:.2f}){tag}", flush=True)
            for name, m in r["models"].items():
                if "average_precision" not in m:
                    print(f"    {name:24s} {m.get('note','')}", flush=True)
                    continue
                print(f"    {name:24s} AP={m['average_precision']:.3f} "
                      f"{m.get('ap_ci95')} lift={m['ap_lift_over_prevalence']:+.3f} "
                      f"rho={m['spearman']:+.3f}"
                      f"{'  [degenerate at 2x]' if m['threshold_degenerate'] else ''}",
                      flush=True)
            print(f"    best by AP: {r['best_by_average_precision']}", flush=True)
        summary["label_sets"][label] = entry
        print("", flush=True)

    path = OUT / "temporal_holdout.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
