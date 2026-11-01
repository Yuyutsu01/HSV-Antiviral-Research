#!/usr/bin/env python
"""
Phase 3, step 6 - Gate 4, evaluated properly.

p3_04 asks only whether a model beats the frequency-only baseline. That turns
out to be the wrong question on its own: under leave-one-gene-out the frequency
baseline scores below zero, so a model can "beat" it while still having no
predictive skill at all. A difference against a broken reference is not a
finding.

Gate 4 therefore requires two things of the same model, both as 95% intervals
over the saved out-of-fold predictions:

  skill       its own Spearman interval excludes zero, i.e. it predicts the
              measured fold change better than chance
  advantage   a paired bootstrap of (model - frequency baseline) excludes zero

Only a model passing both can be said to add anything. Passing the second
alone is reported explicitly as not passing.

Reads the out-of-fold predictions written by p3_04, so nothing is refitted and
no fold assignment can differ from the one that produced them.

Outputs:
  reports/phase3/gate4_summary.json
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase3"

SEED = 20261102
N_BOOT = 5000
FREQ_BASELINE = "baseline:frequency_only"


def ci(vals, lo=2.5, hi=97.5):
    return [round(float(x), 4) for x in np.percentile(vals, [lo, hi])]


def analyse(path, rng):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    y = np.array([float(r["log10_fc"]) for r in rows])
    binary = np.array([int(r["resistant"]) for r in rows])
    cols = [c for c in rows[0] if "|" in c]
    preds = {}
    for c in cols:
        scheme, name, model = c.split("|")
        v = np.array([np.nan if r[c] == "" else float(r[c]) for r in rows])
        preds[(scheme, name, model)] = v

    out = {}
    for scheme in sorted({k[0] for k in preds}):
        # Pick the better of the two estimators for the frequency baseline, so
        # the reference is the strongest version of it rather than a weak one.
        freq_keys = [k for k in preds
                     if k[0] == scheme and k[1] == FREQ_BASELINE]
        def rho_of(k):
            p = preds[k]
            ok = ~np.isnan(p)
            if ok.sum() < 10 or np.std(p[ok]) < 1e-12:
                return -np.inf
            r = spearmanr(y[ok], p[ok]).statistic
            return -np.inf if np.isnan(r) else r
        fkey = max(freq_keys, key=rho_of) if freq_keys else None
        pf = preds[fkey] if fkey else None

        entries = {}
        for k, p in preds.items():
            if k[0] != scheme:
                continue
            ok = ~np.isnan(p)
            if ok.sum() < 10 or np.std(p[ok]) < 1e-12:
                continue
            yv, pv, bv = y[ok], p[ok], binary[ok]
            rho = spearmanr(yv, pv).statistic
            auc = (roc_auc_score(bv, pv) if len(set(bv)) == 2 else None)

            # Bootstrap this model's own skill, and its paired advantage over
            # the frequency baseline on identical resampled rows.
            okb = ok & ~np.isnan(pf) if pf is not None else ok
            yb, pb = y[okb], p[okb]
            fb = pf[okb] if pf is not None else None
            rs, ds, aucs = [], [], []
            n = len(yb)
            for _ in range(N_BOOT):
                idx = rng.integers(0, n, n)
                if np.std(pb[idx]) < 1e-12:
                    continue
                a = spearmanr(yb[idx], pb[idx]).statistic
                if np.isnan(a):
                    continue
                rs.append(a)
                bb = binary[okb][idx]
                if len(set(bb)) == 2:
                    aucs.append(roc_auc_score(bb, pb[idx]))
                if fb is not None and np.std(fb[idx]) > 1e-12:
                    c = spearmanr(yb[idx], fb[idx]).statistic
                    if not np.isnan(c):
                        ds.append(a - c)
            e = {
                "n": int(ok.sum()),
                "spearman": round(float(rho), 4),
                "spearman_ci": ci(rs) if len(rs) > 50 else None,
                "auc": None if auc is None else round(float(auc), 4),
                "auc_ci": ci(aucs) if len(aucs) > 50 else None,
                "rmse": round(float(np.sqrt(np.mean((yv - pv) ** 2))), 4),
            }
            e["diff_vs_frequency"] = (round(float(rho - rho_of(fkey)), 4)
                                      if fkey else None)
            e["diff_ci"] = ci(ds) if len(ds) > 50 else None
            e["has_skill"] = bool(e["spearman_ci"] and e["spearman_ci"][0] > 0)
            e["beats_frequency"] = bool(e["diff_ci"] and e["diff_ci"][0] > 0)
            e["passes_gate_4"] = bool(e["has_skill"] and e["beats_frequency"])
            entries[f"{k[1]}/{k[2]}"] = e

        passing = {k: v for k, v in entries.items() if v["passes_gate_4"]}
        with_skill = {k: v for k, v in entries.items() if v["has_skill"]}
        best = max(entries.items(), key=lambda kv: kv[1]["spearman"],
                   default=(None, None))
        out[scheme] = {
            "frequency_baseline": (f"{fkey[1]}/{fkey[2]}" if fkey else None),
            "frequency_baseline_spearman": (round(float(rho_of(fkey)), 4)
                                            if fkey else None),
            "n_models_evaluated": len(entries),
            "n_with_skill_above_chance": len(with_skill),
            "n_passing_gate_4": len(passing),
            "models_passing_gate_4": sorted(passing),
            "best_by_spearman": best[0],
            "best_spearman": best[1]["spearman"] if best[1] else None,
            "best_spearman_ci": best[1]["spearman_ci"] if best[1] else None,
            "gate_4_passed": len(passing) > 0,
            "models": entries,
        }
    return out


def main():
    REP.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    result = {}
    for path in sorted(PROC.glob("model_predictions_*.tsv")):
        label = path.stem.replace("model_predictions_", "")
        print(f"\n################ {label} ################", flush=True)
        res = analyse(path, rng)
        result[label] = res
        for scheme, r in res.items():
            print(f"\n=== {scheme} ===")
            print(f"  frequency baseline: {r['frequency_baseline']} "
                  f"rho={r['frequency_baseline_spearman']:+.3f}")
            print(f"  {'model':34s} {'rho':>7s} {'95% CI':>16s} "
                  f"{'auc':>6s} {'skill':>6s} {'beats':>6s} {'GATE':>5s}")
            for name, e in sorted(r["models"].items(),
                                  key=lambda kv: -kv[1]["spearman"]):
                c = e["spearman_ci"] or [float("nan")] * 2
                auc_txt = "n/a" if e["auc"] is None else format(e["auc"], ".3f")
                gate_txt = "PASS" if e["passes_gate_4"] else "-"
                print(f"  {name:34s} {e['spearman']:+7.3f} "
                      f"[{c[0]:+.3f},{c[1]:+.3f}] "
                      f"{auc_txt:>6s} "
                      f"{str(e['has_skill']):>6s} {str(e['beats_frequency']):>6s} "
                      f"{gate_txt:>5s}")
            print(f"  -> Gate 4 {'PASSED' if r['gate_4_passed'] else 'NOT PASSED'}"
                  f" ({r['n_with_skill_above_chance']}/{r['n_models_evaluated']}"
                  f" models above chance, {r['n_passing_gate_4']} passing)")

    (REP / "gate4_summary.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
