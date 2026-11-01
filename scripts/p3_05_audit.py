#!/usr/bin/env python
"""
Phase 3, step 5 - circularity audit.

A model can score well here for reasons that have nothing to do with
predicting resistance. Three specific ways, each with a test:

permutation_tests
    Shuffle the target within cross-validation groups, many times, and build
    the null distribution of the score. Two things come out of it. First, a
    leakage check: if noise alone scores positively, information is crossing
    the fold boundary and every other number in Phase 3 is void. Second, the
    reference the real scores should be judged against, which is NOT zero -
    leave-one-gene-out predicts about the training mean for a held-out gene,
    and removing a gene shifts that mean away from it, so the scheme produces
    systematic negative correlation from noise. Measured here at about -0.21,
    which is why several baselines looked "worse than chance" when scored
    against zero. Every reported score therefore carries a permutation p.

conservation_matched
    Resistance mutations are catalogued at conserved, buried, active-site
    positions, and so is every other functionally important residue. A model
    that has only learned "this position is conserved" will look predictive.
    Negatives are resampled to match the positives' conservation distribution,
    removing that shortcut. If performance collapses, the apparent skill was
    conservation all along.

esm_residualised
    ESM-2 is trained on evolutionary data, so its score partly restates the
    conservation features. The ESM score is regressed on conservation and
    structure, and only the residual is given to the model. What survives is
    ESM's independent contribution, if any.

Outputs:
  reports/phase3/audit_summary.json
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p3_04_baselines import (                                   # noqa: E402
    FEATURE_SETS, STAGES, SEED, evaluate, groups_for, load_dataset,
    make_models, matrix, nested_predict,
)

ROOT = Path(__file__).resolve().parent.parent
REP = ROOT / "reports" / "phase3"

N_MATCHED_DRAWS = 100
# 100 replicates puts the smallest attainable p at 1/101, which is below 0.05
# with room to spare. Going to 200 only sharpens a p-value that is already at
# its floor for the feature sets that matter, and quadruples the runtime of the
# forest targets.
N_PERMUTATIONS = 100
# The null does not need a large ensemble; 100 trees is enough to characterise
# it, and observed and permuted statistics use the same setting either way.
PERM_FOREST_TREES = 100
CONS_FEATURE = "cons_site_variability"

# Feature sets worth a permutation test, with a single fixed hyperparameter
# each. Ridge for the linear view, a mildly regularised forest for the
# non-linear one; no tuning, because the null must come from the same procedure
# as the observed statistic.
_RIDGE = ("ridge", {"m__alpha": 1.0})
_FOREST = ("forest", {"m__min_samples_leaf": 5, "m__max_features": 0.5,
                      "m__n_estimators": PERM_FOREST_TREES})
PERM_TARGETS = [
    ("frequency_only", FEATURE_SETS["frequency_only"], *_RIDGE),
    ("conservation_only", FEATURE_SETS["conservation_only"], *_RIDGE),
    ("structure_only", FEATURE_SETS["structure_only"], *_RIDGE),
    ("structure_only_forest", FEATURE_SETS["structure_only"], *_FOREST),
    ("esm_only", FEATURE_SETS["esm_only"], *_RIDGE),
    ("esm_only_forest", FEATURE_SETS["esm_only"], *_FOREST),
    ("full_model", STAGES["stage5_plus_esm"], *_RIDGE),
    ("full_model_forest", STAGES["stage5_plus_esm"], *_FOREST),
]


def run_once(rows, feats, y, binary, scheme, rng, model="ridge", X=None):
    groups = groups_for(rows, scheme)
    if X is None:
        X = matrix(rows, feats)
    est, grid = make_models()[model]
    pred, _ = nested_predict(X, y, groups, est, grid, scheme)
    return evaluate(y, pred, binary, rng)


def cv_rho(X, y, groups, scheme, model_name, params):
    """Grouped out-of-fold Spearman with FIXED hyperparameters.

    No inner tuning loop: choosing hyperparameters against a permuted target is
    meaningless, and a permutation test is only valid if the observed statistic
    and the null are produced by exactly the same procedure. So both sides use
    this function."""
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
    r = spearmanr(y[ok], pred[ok]).statistic
    return None if np.isnan(r) else float(r)


def permutation_test(rows, y, scheme, rng, feats, model_name, params, n_rep):
    """Shuffle the target WITHIN each cross-validation group.

    Within-group shuffling preserves each group's mean label, so the null keeps
    whatever bias the scheme itself produces. That matters here: leave-one-gene-
    out predicts roughly the training mean for a held-out gene, and removing a
    gene shifts that mean away from it, which generates systematic negative
    correlation from noise alone. A null centred on zero would wrongly treat
    that bias as signal."""
    groups = groups_for(rows, scheme)
    X = matrix(rows, feats)
    obs = cv_rho(X, y, groups, scheme, model_name, params)
    if obs is None:
        return {"note": "observed statistic degenerate"}
    null = []
    for _ in range(n_rep):
        yp = y.copy()
        for g in set(groups):
            idx = np.where(groups == g)[0]
            yp[idx] = rng.permutation(yp[idx])
        r = cv_rho(X, yp, groups, scheme, model_name, params)
        if r is not None:
            null.append(r)
    if len(null) < 20:
        return {"note": f"only {len(null)} usable replicates"}
    arr = np.array(null)
    # One-sided: how often does noise reach the observed value or beyond.
    p = (1 + int(np.sum(arr >= obs))) / (1 + len(arr))
    return {
        "model": f"{model_name}{params}",
        "observed_spearman": round(obs, 4),
        "n_replicates": len(null),
        "null_mean": round(float(arr.mean()), 4),
        "null_sd": round(float(arr.std()), 4),
        "null_ci95": [round(float(x), 4)
                      for x in np.percentile(arr, [2.5, 97.5])],
        "null_p95": round(float(np.percentile(arr, 95)), 4),
        "permutation_p": round(p, 5),
        "significant_vs_null": bool(p < 0.05),
        # Beating the null is not the same as predicting anything. Under
        # leave-one-gene-out the null sits near -0.25, so a feature set can be
        # "significant" at rho = -0.21 while still ranking variants backwards.
        # Only a positive observed correlation that also clears its null is
        # usable, and the two must never be reported as one claim.
        "usable": bool(p < 0.05 and obs > 0),
        "excess_over_null_mean": round(obs - float(arr.mean()), 4),
    }


def conservation_matched(rows, y, binary, scheme, rng):
    """Equalise the conservation distribution of positives and negatives, then
    re-evaluate on the matched subset.

    Both classes are subsampled to the smaller count in each conservation bin.
    Drawing only negatives to match positives does not work here: negatives are
    the minority (32 against 87), so every request exceeds the pool and the
    draw silently returns the whole dataset, which tests nothing."""
    cons = np.array([float(r[CONS_FEATURE]) if r[CONS_FEATURE] != "" else np.nan
                     for r in rows])
    med = np.nanmedian(cons)
    cons = np.where(np.isnan(cons), med, cons)
    pos = np.where(binary == 1)[0]
    neg = np.where(binary == 0)[0]
    if len(neg) < 10 or len(pos) < 10:
        return {"note": f"too few rows to match ({len(pos)}+/{len(neg)}-)"}

    # Bins over the whole conservation range, so neither class defines the grid.
    qs = np.quantile(cons, [0, 0.25, 0.5, 0.75, 1.0])
    edges = np.unique(qs[1:-1])
    pos_bin = np.digitize(cons[pos], edges)
    neg_bin = np.digitize(cons[neg], edges)
    n_bins = len(edges) + 1
    # Per bin, take the same number from each class.
    per_bin = {b: min(int(np.sum(pos_bin == b)), int(np.sum(neg_bin == b)))
               for b in range(n_bins)}
    if sum(per_bin.values()) < 15:
        return {"note": "conservation bins do not overlap enough to match",
                "per_bin": per_bin}

    feats = STAGES["stage5_plus_esm"]
    X_all = matrix(rows, feats)
    rhos, aucs, sizes = [], [], []
    for _ in range(N_MATCHED_DRAWS):
        keep = []
        for b, k in per_bin.items():
            if k == 0:
                continue
            keep += list(rng.choice(pos[pos_bin == b], size=k, replace=False))
            keep += list(rng.choice(neg[neg_bin == b], size=k, replace=False))
        keep = np.array(sorted(set(keep)))
        sub = [rows[i] for i in keep]
        if len(set(groups_for(sub, scheme))) < 3 or len(set(binary[keep])) < 2:
            continue
        ev = run_once(sub, feats, y[keep], binary[keep], scheme, rng,
                      X=X_all[keep])
        if ev.get("spearman") is not None:
            rhos.append(ev["spearman"])
            sizes.append(len(keep))
        if ev.get("auc") is not None:
            aucs.append(ev["auc"])
    if not rhos:
        return {"note": "no usable matched draws"}
    out = {
        "n_draws": len(rhos),
        "median_subset_size": int(np.median(sizes)),
        "spearman_mean": round(float(np.mean(rhos)), 4),
        "spearman_ci": [round(float(x), 4)
                        for x in np.percentile(rhos, [2.5, 97.5])],
    }
    if aucs:
        out["auc_mean"] = round(float(np.mean(aucs)), 4)
        out["auc_ci"] = [round(float(x), 4)
                         for x in np.percentile(aucs, [2.5, 97.5])]
    return out


def esm_residualised(rows, y, binary, scheme, rng):
    """Replace the ESM features with their residual after regressing on
    conservation and structure, so only ESM's independent part remains."""
    esm_feats = FEATURE_SETS["esm_only"]
    confound = (FEATURE_SETS["conservation_only"]
                + FEATURE_SETS["structure_only"])
    E = matrix(rows, esm_feats)
    C = matrix(rows, confound)
    imp = SimpleImputer(strategy="median")
    C = imp.fit_transform(C)
    E = SimpleImputer(strategy="median").fit_transform(E)

    resid = np.zeros_like(E)
    r2 = []
    for j in range(E.shape[1]):
        lr = LinearRegression().fit(C, E[:, j])
        resid[:, j] = E[:, j] - lr.predict(C)
        r2.append(round(float(lr.score(C, E[:, j])), 4))

    full = STAGES["stage5_plus_esm"]
    X_full = matrix(rows, full)
    X_resid = X_full.copy()
    for k, f in enumerate(full):
        if f in esm_feats:
            X_resid[:, k] = resid[:, esm_feats.index(f)]

    ev_full = run_once(rows, full, y, binary, scheme, rng, X=X_full)
    ev_resid = run_once(rows, full, y, binary, scheme, rng, X=X_resid)
    ev_esm_only = run_once(rows, esm_feats, y, binary, scheme, rng, X=E)
    ev_esm_resid = run_once(rows, esm_feats, y, binary, scheme, rng, X=resid)
    return {
        "esm_variance_explained_by_confounds_r2": dict(zip(esm_feats, r2)),
        "full_model": {k: ev_full.get(k) for k in ("spearman", "auc", "rmse")},
        "full_model_esm_residualised": {k: ev_resid.get(k)
                                        for k in ("spearman", "auc", "rmse")},
        "esm_alone": {k: ev_esm_only.get(k) for k in ("spearman", "auc")},
        "esm_alone_residualised": {k: ev_esm_resid.get(k)
                                   for k in ("spearman", "auc")},
    }


def main():
    REP.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows = load_dataset(("A",))
    y = np.array([float(r["log10_fc"]) for r in rows])
    binary = np.array([int(r["resistant"]) for r in rows])
    print(f"audit on {len(rows)} Tier A rows "
          f"({binary.sum()} resistant, {len(rows) - binary.sum()} negative)\n")

    out = {"n": len(rows), "seed": SEED, "n_permutations": N_PERMUTATIONS}
    for scheme in ("leave_one_gene_out", "position_blocks"):
        print(f"=== {scheme} ===", flush=True)
        s = {"permutation_tests": {}}
        # The two circularity tests run before the permutation nulls: they are
        # the ones the phase's conclusion turns on, and they are much cheaper.
        print("  conservation-matched classes ...", flush=True)
        s["conservation_matched"] = conservation_matched(
            rows, y, binary, scheme, rng)
        print("   ", s["conservation_matched"], flush=True)
        print("  ESM residualised ...", flush=True)
        s["esm_residualised"] = esm_residualised(rows, y, binary, scheme, rng)
        print("   ", json.dumps(s["esm_residualised"]["full_model"]),
              "->", json.dumps(s["esm_residualised"]["full_model_esm_residualised"]),
              flush=True)
        for fname, feats, model, params in PERM_TARGETS:
            r = permutation_test(rows, y, scheme, rng, feats, model, params,
                                 N_PERMUTATIONS)
            s["permutation_tests"][fname] = r
            if "observed_spearman" in r:
                print(f"  perm {fname:22s} obs={r['observed_spearman']:+.3f} "
                      f"null={r['null_mean']:+.3f}+-{r['null_sd']:.3f} "
                      f"p95={r['null_p95']:+.3f} "
                      f"p={r['permutation_p']:.4f} "
                      f"{'SIGNIFICANT' if r['significant_vs_null'] else 'not significant'}",
                      flush=True)
            else:
                print(f"  perm {fname:22s} {r}", flush=True)
        out[scheme] = s
        print(flush=True)

    # Leakage shows up as a null centred ABOVE zero: noise alone scoring
    # positively means information is crossing the fold boundary. A null below
    # zero is the scheme's own pessimistic bias, which is not leakage and is
    # handled by scoring against the null instead of against zero.
    leak = {}
    for sc in ("leave_one_gene_out", "position_blocks"):
        means = [t["null_mean"] for t in out[sc]["permutation_tests"].values()
                 if "null_mean" in t]
        leak[sc] = {
            "max_null_mean": round(max(means), 4) if means else None,
            "min_null_mean": round(min(means), 4) if means else None,
            "leakage_suspected": bool(means and max(means) > 0.10),
            "scheme_biased_negative": bool(means and min(means) < -0.10),
        }
        print(f"{sc}: null means {leak[sc]['min_null_mean']} .. "
              f"{leak[sc]['max_null_mean']} -> "
              f"{'LEAKAGE SUSPECTED' if leak[sc]['leakage_suspected'] else 'no leakage'}"
              f"{', scheme biased negative' if leak[sc]['scheme_biased_negative'] else ''}")
    out["leakage_check"] = leak

    schemes = ("leave_one_gene_out", "position_blocks")
    sig = {sc: sorted(k for k, t in out[sc]["permutation_tests"].items()
                      if t.get("significant_vs_null")) for sc in schemes}
    use = {sc: sorted(k for k, t in out[sc]["permutation_tests"].items()
                      if t.get("usable")) for sc in schemes}
    out["significant_feature_sets"] = sig
    out["usable_feature_sets"] = use
    print("\nsignificant against their own permutation null:")
    for sc in schemes:
        print(f"  {sc}: {sig[sc] if sig[sc] else 'NONE'}")
    print("\nUSABLE (positive correlation AND clears its null):")
    for sc in schemes:
        print(f"  {sc}: {use[sc] if use[sc] else 'NONE'}")
        only_sig = sorted(set(sig[sc]) - set(use[sc]))
        if only_sig:
            print(f"    beat their null but predict backwards: {only_sig}")

    (REP / "audit_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
