"""Phase 4, step 2 — the matched-negative gap as a headline number, and the
feature correlation matrix published in full.

The protocol asks for both of these in the main text, the gap because it is
"the quantitative measure of how much apparent performance was circular".

Two differences from the Phase 3 audit, which computed something close to this
but not the thing the protocol specified:

1. **The control is size-matched.** The audit compared a conservation-matched
   subset of ~64 rows against the full 119. That difference mixes two causes:
   the matching, and the halved sample. Here each matched draw is paired with a
   random draw taking the same number of positives and negatives while ignoring
   conservation, so the only thing that differs between the two arms is whether
   the conservation distributions were equalised.
2. **It is scored on average precision**, the pre-declared primary metric,
   rather than Spearman and AUROC.

The gap is reported as random minus matched, so a positive gap means apparent
performance was inflated by the conservation shortcut.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p3_04_baselines import (                                   # noqa: E402
    FEATURE_SETS, SEED, STAGES, groups_for, load_dataset, make_models, matrix,
)
from p4_01_metrics import average_precision, precision_at_recall  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "reports" / "phase4"
PROC = ROOT / "data" / "processed"

N_DRAWS = 100
CONS_FEATURE = "cons_site_variability"
SCHEMES = ["leave_one_gene_out", "position_blocks"]

# Fixed hyperparameters, as in the audit: the two arms must be produced by an
# identical procedure or the difference between them is not interpretable.
MODELS = [
    ("structure_only", FEATURE_SETS["structure_only"], "ridge",
     {"m__alpha": 1.0}),
    ("structure_only_forest", FEATURE_SETS["structure_only"], "forest",
     {"m__min_samples_leaf": 5, "m__max_features": 0.5, "m__n_estimators": 100}),
    ("full_model", STAGES["stage5_plus_esm"], "ridge", {"m__alpha": 1.0}),
    ("full_model_forest", STAGES["stage5_plus_esm"], "forest",
     {"m__min_samples_leaf": 5, "m__max_features": 0.5, "m__n_estimators": 100}),
]


def oof_predict(X, y, groups, model_name, params, scheme):
    outer = (LeaveOneGroupOut() if scheme == "leave_one_gene_out"
             else GroupKFold(n_splits=min(5, len(set(groups)))))
    est, _ = make_models()[model_name]
    est = est.set_params(**params)
    pred = np.full(len(y), np.nan)
    for tr, te in outer.split(X, y, groups):
        est.fit(X[tr], y[tr])
        pred[te] = est.predict(X[te])
    ok = ~np.isnan(pred)
    if ok.sum() < 10 or np.std(pred[ok]) < 1e-12:
        return None
    return pred, ok


def score(y_cont, binary, pred, ok):
    yb = binary[ok]
    if len(set(yb)) < 2:
        return None
    return {
        "average_precision": float(average_precision(yb, pred[ok])),
        "precision_at_95_recall": float(precision_at_recall(yb, pred[ok])),
        "spearman": float(spearmanr(y_cont[ok], pred[ok]).statistic),
        "prevalence": float(yb.mean()),
    }


def bin_structure(rows, binary):
    cons = np.array([float(r[CONS_FEATURE]) if r[CONS_FEATURE] != "" else np.nan
                     for r in rows])
    cons = np.where(np.isnan(cons), np.nanmedian(cons), cons)
    pos = np.where(binary == 1)[0]
    neg = np.where(binary == 0)[0]
    edges = np.unique(np.quantile(cons, [0.25, 0.5, 0.75]))
    pos_bin = np.digitize(cons[pos], edges)
    neg_bin = np.digitize(cons[neg], edges)
    n_bins = len(edges) + 1
    per_bin = {b: min(int(np.sum(pos_bin == b)), int(np.sum(neg_bin == b)))
               for b in range(n_bins)}
    return cons, pos, neg, pos_bin, neg_bin, per_bin


def gap_for(rows, y, binary, scheme, model_name, feats, est_name, params, rng):
    cons, pos, neg, pos_bin, neg_bin, per_bin = bin_structure(rows, binary)
    if sum(per_bin.values()) < 15:
        return {"note": "conservation bins do not overlap enough to match"}

    X_all = matrix(rows, feats)
    n_pos_target = sum(per_bin.values())
    n_neg_target = n_pos_target
    arms = {"matched": [], "random": []}
    sizes = []

    for _ in range(N_DRAWS):
        # ---- arm 1: conservation-matched, equal counts per bin per class
        keep_m = []
        for b, k in per_bin.items():
            if k == 0:
                continue
            keep_m += list(rng.choice(pos[pos_bin == b], size=k, replace=False))
            keep_m += list(rng.choice(neg[neg_bin == b], size=k, replace=False))
        keep_m = np.array(sorted(set(keep_m)))

        # ---- arm 2: same class counts, drawn without regard to conservation
        keep_r = np.array(sorted(
            list(rng.choice(pos, size=min(n_pos_target, len(pos)), replace=False))
            + list(rng.choice(neg, size=min(n_neg_target, len(neg)), replace=False))))

        for arm, keep in (("matched", keep_m), ("random", keep_r)):
            sub = [rows[i] for i in keep]
            if len(set(groups_for(sub, scheme))) < 3 or len(set(binary[keep])) < 2:
                continue
            got = oof_predict(X_all[keep], y[keep], groups_for(sub, scheme),
                              est_name, params, scheme)
            if got is None:
                continue
            pred, ok = got
            s = score(y[keep], binary[keep], pred, ok)
            if s:
                arms[arm].append(s)
                if arm == "matched":
                    sizes.append(len(keep))

    if not arms["matched"] or not arms["random"]:
        return {"note": "no usable draws"}

    def agg(xs, key):
        a = np.array([x[key] for x in xs])
        return {"mean": round(float(a.mean()), 4),
                "ci95": [round(float(v), 4)
                         for v in np.percentile(a, [2.5, 97.5])]}

    out = {"n_draws_matched": len(arms["matched"]),
           "n_draws_random": len(arms["random"]),
           "median_subset_size": int(np.median(sizes)) if sizes else None,
           "per_bin_counts": per_bin}
    for key in ("average_precision", "precision_at_95_recall", "spearman"):
        m, r = agg(arms["matched"], key), agg(arms["random"], key)
        out[key] = {"matched": m, "random": r,
                    "gap_random_minus_matched": round(r["mean"] - m["mean"], 4)}
    # Prevalence is 0.5 by construction in both arms, so AP is comparable
    # between them and against a 0.5 floor rather than the 0.731 of the
    # unmatched set. Recorded so that is not left implicit.
    out["prevalence_both_arms"] = agg(arms["matched"], "prevalence")["mean"]
    return out


def correlation_matrix(rows):
    feats = sorted({f for s in FEATURE_SETS.values() for f in s})
    X = matrix(rows, feats)
    n = len(feats)
    rho = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if np.std(X[:, i]) < 1e-12 or np.std(X[:, j]) < 1e-12:
                continue
            rho[i, j] = spearmanr(X[:, i], X[:, j]).statistic
    path = PROC / "feature_correlation_matrix.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["feature"] + feats)
        for i, f in enumerate(feats):
            w.writerow([f] + ["" if np.isnan(v) else round(float(v), 4)
                              for v in rho[i]])
    # The pairs that matter for the circularity argument.
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if not np.isnan(rho[i, j]) and abs(rho[i, j]) >= 0.5:
                pairs.append({"a": feats[i], "b": feats[j],
                              "spearman": round(float(rho[i, j]), 4)})
    pairs.sort(key=lambda p: -abs(p["spearman"]))
    return {"n_features": n, "path": str(path.relative_to(ROOT)),
            "pairs_above_0.5": pairs}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows = load_dataset(("A",))
    y = np.array([float(r["log10_fc"]) for r in rows])
    binary = np.array([int(r["resistant"]) for r in rows])
    print(f"primary set: n={len(rows)} positives={int(binary.sum())}", flush=True)

    summary = {"seed": SEED, "n_draws": N_DRAWS, "n": len(rows),
               "primary_metric": "average_precision",
               "control": "size- and class-matched random draw",
               "correlation_matrix": correlation_matrix(rows),
               "schemes": {}}

    print("\nfeature pairs with |rho| >= 0.5:", flush=True)
    for p in summary["correlation_matrix"]["pairs_above_0.5"]:
        print(f"  {p['a']:28s} {p['b']:28s} {p['spearman']:+.3f}", flush=True)

    for scheme in SCHEMES:
        print(f"\n################ {scheme} ################", flush=True)
        summary["schemes"][scheme] = {}
        for name, feats, est_name, params in MODELS:
            g = gap_for(rows, y, binary, scheme, name, feats, est_name,
                        params, rng)
            summary["schemes"][scheme][name] = g
            if "note" in g:
                print(f"  {name:24s} {g['note']}", flush=True)
                continue
            ap = g["average_precision"]
            print(f"  {name:24s} n={g['median_subset_size']:3d}  "
                  f"AP random {ap['random']['mean']:.3f} "
                  f"vs matched {ap['matched']['mean']:.3f}  "
                  f"GAP {ap['gap_random_minus_matched']:+.3f}", flush=True)

    path = OUT / "circularity_gap.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
