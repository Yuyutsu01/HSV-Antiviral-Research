"""Phase 4, step 4 — cross-virus transfer, the second pre-declared external
validation.

Train on HSV-1 entries, test on HSV-2, and the reverse. The protocol calls this
"biologically meaningful, since the orthologous targets differ modestly in
sequence, and an honest test of generalisation", conditional on the HSV-2 label
count being adequate.

It is a single held-out split, not cross-validation, so there is no fold
structure and no inner tuning loop: hyperparameters are the fixed, declared
values. Tuning them would mean selecting on the test virus.

The companion validation, the temporal hold-out, is NOT implemented here.
HerpesDRG supplies a DOI for every row but publication year is only recoverable
offline for 35 of 89 distinct studies, and the recoverable ones cluster by
publisher -- ASM and Elsevier encode the year in the DOI, Oxford and Wiley do
not. Splitting on a year available for one group of journals and not another
would confound era with journal, so it is left undone and declared as such
rather than approximated. It needs a DOI-to-year resolution step.

One asymmetry to keep in mind when reading the result: HSV-2 structural
features are inherited from HSV-1 structures by pairwise alignment, because no
HSV-2 structures exist for these genes. So the HSV-1 -> HSV-2 direction tests
transfer of a model onto features that were themselves derived from the
training virus's geometry. That makes a good score less impressive than it
looks, and it is the direction with the larger test set.
"""

import json
import sys
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

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "reports" / "phase4"
N_BOOTSTRAP = 2000

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


def score(y_bin, y_cont, pred, rng):
    if len(set(y_bin)) < 2:
        return {"note": "test set is single-class"}
    conf = confusion_at_2x(y_bin, pred)
    ap = average_precision(y_bin, pred)
    boot_ap, boot_rho = [], []
    idx = np.arange(len(y_bin))
    for _ in range(N_BOOTSTRAP):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(set(y_bin[b])) < 2:
            continue
        boot_ap.append(average_precision(y_bin[b], pred[b]))
        r = spearmanr(y_cont[b], pred[b]).statistic
        if not np.isnan(r):
            boot_rho.append(float(r))
    out = {
        "n_test": int(len(y_bin)),
        "prevalence": round(float(y_bin.mean()), 4),
        "average_precision": round(float(ap), 4),
        "ap_lift_over_prevalence": round(float(ap - y_bin.mean()), 4),
        "precision_at_95_recall": round(float(precision_at_recall(y_bin, pred)), 4),
        "spearman": round(float(spearmanr(y_cont, pred).statistic), 4),
        "confusion_at_2x": conf,
        "threshold_degenerate": bool(
            conf["n_predicted_resistant"] in (0, len(y_bin))),
    }
    if boot_ap:
        out["ap_ci95"] = [round(float(v), 4)
                          for v in np.percentile(boot_ap, [2.5, 97.5])]
    if boot_rho:
        out["spearman_ci95"] = [round(float(v), 4)
                                for v in np.percentile(boot_rho, [2.5, 97.5])]
    return out


def run(rows, y, binary, train_virus, test_virus, rng):
    tr = np.array([i for i, r in enumerate(rows) if r["virus"] == train_virus])
    te = np.array([i for i, r in enumerate(rows) if r["virus"] == test_virus])
    res = {"n_train": int(len(tr)), "n_test": int(len(te))}
    if len(tr) < 20 or len(te) < 10:
        res["note"] = "too few labels for a transfer test"
        return res
    if len(set(binary[te])) < 2:
        res["note"] = "test virus is single-class"
        return res

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
    best = max((k for k, v in res["models"].items()
                if "average_precision" in v),
               key=lambda k: res["models"][k]["average_precision"], default=None)
    res["best_by_average_precision"] = best
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    summary = {"seed": SEED, "n_bootstrap": N_BOOTSTRAP,
               "primary_metric": "average_precision",
               "temporal_holdout": "not implemented; see module docstring",
               "label_sets": {}}

    for label, tiers in LABEL_SETS.items():
        rows = load_dataset(tiers)
        y = np.array([float(r["log10_fc"]) for r in rows])
        binary = np.array([int(r["resistant"]) for r in rows])
        counts = {v: sum(1 for r in rows if r["virus"] == v)
                  for v in ("HSV1", "HSV2")}
        print(f"\n################ {label}  n={len(rows)}  {counts}",
              flush=True)

        entry = {"n": len(rows), "by_virus": counts, "directions": {}}
        for a, b in (("HSV1", "HSV2"), ("HSV2", "HSV1")):
            r = run(rows, y, binary, a, b, rng)
            entry["directions"][f"{a}_to_{b}"] = r
            print(f"\n=== train {a} (n={r['n_train']}) "
                  f"-> test {b} (n={r['n_test']})", flush=True)
            if "note" in r:
                print(f"    {r['note']}", flush=True)
                continue
            print(f"    prevalence in test (AP floor): "
                  f"{r['models'][MODELS[0][0]]['prevalence']:.3f}", flush=True)
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

    path = OUT / "cross_virus.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
