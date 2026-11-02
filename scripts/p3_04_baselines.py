#!/usr/bin/env python
"""
Phase 3, step 4 - the six baselines and the staged models.

Aim 4 asks whether a sequence-based model predicts measured resistance better
than knowing how common the allele is. The frequency-only baseline is therefore
not a throwaway control, it is the thing to beat, and it is evaluated under the
same cross-validation as everything else.

Two cross-validation schemes, because they answer different questions and the
easier one is the one usually reported without saying so:

  leave_one_gene_out   train on other genes, predict a held-out gene. Tests
                       whether the model learned anything transferable. Hard,
                       and the honest test for a method claiming to generalise.
  position_blocks      group by 50-residue blocks within a gene, so the model
                       sees the same gene in training but never a neighbouring
                       residue. Tests prediction of new variants in a gene that
                       is already characterised.

Both are grouped, so no variant's immediate neighbours ever appear in the
training fold. Hyperparameters are chosen in an inner grouped loop, so the
reported numbers are not tuned on the test fold.

Confidence intervals come from a bootstrap over outer-fold predictions. With
n around 120 they are wide, and that width is a finding rather than a nuisance.

Outputs:
  data/processed/model_predictions.tsv
  reports/phase3/baselines_summary.json
"""

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase3"

SEED = 20261102
N_BOOTSTRAP = 2000
BLOCK = 50

FEATURE_SETS = {
    "frequency_only": ["freq_frequency", "freq_count", "freq_observed"],
    "conservation_only": ["cons_site_variability", "cons_site_coverage"],
    "blosum_only": ["prop_blosum62"],
    "properties_only": ["prop_blosum62", "prop_d_hydropathy", "prop_d_volume",
                        "prop_d_charge", "prop_polarity_switch",
                        "prop_rel_position", "prop_to_proline",
                        "prop_from_glycine"],
    "structure_only": ["struct_rsa", "struct_dist_to_ligand"],
    "esm_only": ["esm_masked_marginal", "esm_ref_logprob"],
}
# Staged models, each adding a group to the one before.
STAGES = {
    "stage1_freq": FEATURE_SETS["frequency_only"],
    "stage2_freq_cons": FEATURE_SETS["frequency_only"]
    + FEATURE_SETS["conservation_only"],
    "stage3_freq_cons_prop": FEATURE_SETS["frequency_only"]
    + FEATURE_SETS["conservation_only"] + FEATURE_SETS["properties_only"],
    "stage4_plus_structure": FEATURE_SETS["frequency_only"]
    + FEATURE_SETS["conservation_only"] + FEATURE_SETS["properties_only"]
    + FEATURE_SETS["structure_only"],
    "stage5_plus_esm": FEATURE_SETS["frequency_only"]
    + FEATURE_SETS["conservation_only"] + FEATURE_SETS["properties_only"]
    + FEATURE_SETS["structure_only"] + FEATURE_SETS["esm_only"],
}


def load_dataset(tiers=("A",)):
    with open(PROC / "model_dataset.tsv", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["tier"] in tiers]
    return rows


def matrix(rows, features):
    X = np.full((len(rows), len(features)), np.nan)
    for i, r in enumerate(rows):
        for j, f in enumerate(features):
            v = r.get(f, "")
            if v != "" and v is not None:
                try:
                    X[i, j] = float(v)
                except ValueError:
                    pass
    return X


def groups_for(rows, scheme):
    if scheme == "leave_one_gene_out":
        return np.array([f"{r['virus']}:{r['gene']}" for r in rows])
    return np.array([f"{r['virus']}:{r['gene']}:{int(r['position']) // BLOCK}"
                     for r in rows])


def make_models():
    """Two estimators: a linear one that cannot overfit 120 rows badly, and a
    flexible one that can, so the comparison shows whether flexibility helps."""
    return {
        "ridge": (Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("sc", StandardScaler()),
            ("m", Ridge(random_state=None)),
        ]), {"m__alpha": [0.1, 1.0, 10.0, 100.0]}),
        "forest": (Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            # n_jobs=1 deliberately: at ~120 rows the worker start-up for each
            # of several hundred nested fits costs far more than the trees do.
            ("m", RandomForestRegressor(n_estimators=300, random_state=SEED,
                                        n_jobs=1)),
        ]), {"m__min_samples_leaf": [1, 5], "m__max_features": [0.5, 1.0]}),
    }


def nested_predict(X, y, groups, estimator, grid, scheme):
    """Outer loop gives honest predictions; inner loop picks hyperparameters
    using only the training fold's groups."""
    outer = (LeaveOneGroupOut() if scheme == "leave_one_gene_out"
             else GroupKFold(n_splits=min(5, len(set(groups)))))
    pred = np.full(len(y), np.nan)
    chosen = []
    for tr, te in outer.split(X, y, groups):
        g_tr = groups[tr]
        n_inner = min(3, len(set(g_tr)))
        best, best_score = None, -np.inf
        if n_inner < 2:
            best = {k: v[0] for k, v in grid.items()}
        else:
            inner = GroupKFold(n_splits=n_inner)
            keys = list(grid)
            import itertools
            for combo in itertools.product(*(grid[k] for k in keys)):
                params = dict(zip(keys, combo))
                scores = []
                for itr, ite in inner.split(X[tr], y[tr], g_tr):
                    est = estimator.set_params(**params)
                    est.fit(X[tr][itr], y[tr][itr])
                    p = est.predict(X[tr][ite])
                    if np.std(p) < 1e-12 or len(set(y[tr][ite])) < 2:
                        scores.append(0.0)
                    else:
                        rho = spearmanr(y[tr][ite], p).statistic
                        scores.append(0.0 if np.isnan(rho) else rho)
                s = float(np.mean(scores))
                if s > best_score:
                    best, best_score = params, s
        estimator.set_params(**best)
        estimator.fit(X[tr], y[tr])
        pred[te] = estimator.predict(X[te])
        chosen.append(best)
    return pred, chosen


def evaluate(y, pred, binary, rng):
    ok = ~np.isnan(pred)
    y, pred, binary = y[ok], pred[ok], binary[ok]
    out = {"n": int(ok.sum())}
    if np.std(pred) < 1e-12:
        out.update(spearman=None, rmse=float(np.sqrt(np.mean((y - pred) ** 2))),
                   auc=None, ap=None, degenerate=True)
        return out
    out["spearman"] = round(float(spearmanr(y, pred).statistic), 4)
    out["rmse"] = round(float(np.sqrt(np.mean((y - pred) ** 2))), 4)
    if len(set(binary)) == 2:
        out["auc"] = round(float(roc_auc_score(binary, pred)), 4)
        out["ap"] = round(float(average_precision_score(binary, pred)), 4)
    else:
        out["auc"] = out["ap"] = None

    boots = {"spearman": [], "auc": []}
    n = len(y)
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, n)
        if np.std(pred[idx]) < 1e-12:
            continue
        rho = spearmanr(y[idx], pred[idx]).statistic
        if not np.isnan(rho):
            boots["spearman"].append(rho)
        if len(set(binary[idx])) == 2:
            boots["auc"].append(roc_auc_score(binary[idx], pred[idx]))
    for k, v in boots.items():
        if len(v) > 50:
            lo, hi = np.percentile(v, [2.5, 97.5])
            out[f"{k}_ci"] = [round(float(lo), 4), round(float(hi), 4)]
    return out


def run_label_set(label_name, tiers, rng):
    rows = load_dataset(tiers)
    if len(rows) < 40:
        print(f"only {len(rows)} rows in {label_name} - not enough to model")
        return None, None, None
    y = np.array([float(r["log10_fc"]) for r in rows])
    binary = np.array([int(r["resistant"]) for r in rows])

    n_obs = sum(int(r["freq_observed"]) for r in rows)
    print(f"\n################ {label_name}: n={len(rows)}, "
          f"resistant={binary.sum()}, negative={len(rows) - binary.sum()}, "
          f"observed in population={n_obs} ################")
    by_gene = Counter(r["virus"] + ":" + r["gene"] for r in rows)
    print("genes:", dict(sorted(by_gene.items())), "\n")

    results, pred_store = {}, {}
    all_sets = {**{f"baseline:{k}": v for k, v in FEATURE_SETS.items()},
                **{f"staged:{k}": v for k, v in STAGES.items()}}

    for scheme in ("leave_one_gene_out", "position_blocks"):
        groups = groups_for(rows, scheme)
        print(f"=== {scheme} ({len(set(groups))} groups) ===", flush=True)
        results[scheme] = {}
        for name, feats in all_sets.items():
            X = matrix(rows, feats)
            if np.all(np.isnan(X)):
                continue
            results[scheme][name] = {}
            for mname, (est, grid) in make_models().items():
                pred, chosen = nested_predict(X, y, groups, est, grid, scheme)
                ev = evaluate(y, pred, binary, rng)
                ev["hyperparams"] = [{k: v for k, v in c.items()} for c in chosen]
                results[scheme][name][mname] = ev
                pred_store[(scheme, name, mname)] = pred
                rho = ev.get("spearman")
                auc = ev.get("auc")
                ci = ev.get("spearman_ci")
                print(f"  {name:28s} {mname:7s} "
                      f"rho={rho if rho is None else f'{rho:+.3f}'} "
                      f"{'' if not ci else f'[{ci[0]:+.2f},{ci[1]:+.2f}]'} "
                      f"auc={auc if auc is None else f'{auc:.3f}'} "
                      f"rmse={ev['rmse']:.3f}", flush=True)
        print(flush=True)

    # Gate 4: does any model beat frequency-only, outside the bootstrap noise?
    gate = {}
    for scheme in results:
        base = max(
            (results[scheme]["baseline:frequency_only"][m].get("spearman") or -1
             for m in ("ridge", "forest")))
        best_name, best_val = None, -np.inf
        for name, bymodel in results[scheme].items():
            if name == "baseline:frequency_only":
                continue
            for m, ev in bymodel.items():
                v = ev.get("spearman")
                if v is not None and v > best_val:
                    best_name, best_val = f"{name}/{m}", v
        # Paired bootstrap on the difference, same resampled rows for both.
        bestkey = None
        for (sc, nm, mm), p in pred_store.items():
            if sc == scheme and best_name == f"{nm}/{mm}":
                bestkey = (sc, nm, mm)
        diffs = []
        if bestkey:
            pb = pred_store[bestkey]
            pf_candidates = [pred_store[(scheme, "baseline:frequency_only", m)]
                             for m in ("ridge", "forest")]
            pf = max(pf_candidates,
                     key=lambda p: (spearmanr(y[~np.isnan(p)], p[~np.isnan(p)])
                                    .statistic if np.std(p[~np.isnan(p)]) > 0
                                    else -1))
            ok = ~np.isnan(pb) & ~np.isnan(pf)
            yv, pbv, pfv = y[ok], pb[ok], pf[ok]
            for _ in range(N_BOOTSTRAP):
                idx = rng.integers(0, len(yv), len(yv))
                if np.std(pbv[idx]) < 1e-12 or np.std(pfv[idx]) < 1e-12:
                    continue
                a = spearmanr(yv[idx], pbv[idx]).statistic
                b = spearmanr(yv[idx], pfv[idx]).statistic
                if not (np.isnan(a) or np.isnan(b)):
                    diffs.append(a - b)
        lo = hi = None
        if len(diffs) > 50:
            lo, hi = (round(float(x), 4)
                      for x in np.percentile(diffs, [2.5, 97.5]))
        gate[scheme] = {
            "frequency_only_spearman": round(float(base), 4),
            "best_other": best_name,
            "best_other_spearman": round(float(best_val), 4),
            "difference": round(float(best_val - base), 4),
            "difference_ci": [lo, hi],
            "beats_frequency_baseline": bool(
                lo is not None and lo > 0),
        }
        print(f"Gate 4 [{scheme}]: frequency rho={base:+.3f}, "
              f"best={best_name} rho={best_val:+.3f}, "
              f"diff={best_val - base:+.3f} "
              f"CI[{lo}, {hi}] -> "
              f"{'BEATS' if gate[scheme]['beats_frequency_baseline'] else 'DOES NOT BEAT'}"
              " frequency baseline")

    return rows, {"n": len(rows), "n_observed_in_population": n_obs,
                  "gate_4": gate, "results": results}, pred_store


def main():
    REP.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    label_sets = {
        "primary_tier_A": ("A",),
        "secondary_all_tiers": ("A", "B", "C"),
    }
    out, stores = {}, {}
    for name, tiers in label_sets.items():
        rows, res, pred_store = run_label_set(name, tiers, rng)
        if res is None:
            continue
        out[name] = res
        stores[name] = (rows, pred_store)

    for name, (rows, pred_store) in stores.items():
        with open(PROC / f"model_predictions_{name}.tsv", "w", newline="",
                  encoding="utf-8") as fh:
            cols = ["virus", "gene", "position", "aa_change", "tier",
                    "log10_fc", "resistant", "freq_observed"]
            keys = sorted(pred_store)
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols + [f"{s}|{n}|{m}" for s, n, m in keys])
            for i, r in enumerate(rows):
                w.writerow([r[c] for c in cols]
                           + ["" if np.isnan(pred_store[k][i])
                              else round(float(pred_store[k][i]), 5)
                              for k in keys])

    (REP / "baselines_summary.json").write_text(
        json.dumps({"seed": SEED, "label_sets": out}, indent=2),
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
