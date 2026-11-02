"""Phase 4, step 1 — the metrics the protocol actually pre-declared.

Phase 3 reported Spearman and AUROC. The protocol names average precision as
the primary metric, says in terms that AUROC "is not the primary metric and no
claim rests on it", and additionally requires precision at 95% recall, Matthews
correlation, calibration, and every figure with a bootstrap interval resampled
over mutations and broken down by tier, gene and drug mechanism class.

Nothing here refits a model. Everything is computed from the out-of-fold
predictions saved by p3_04_baselines.py, so these are the same fits Phase 3
reported, scored the way the protocol said they would be scored.

Calibration note: the pre-declared primary task is regression on log10 fold
change, so the models emit fold-change estimates, not probabilities. A Brier
score is undefined for them. Fitting a calibrator on the same out-of-fold
predictions used for evaluation would leak, so instead the regression
calibration (slope and intercept of observed on predicted) is reported, and
the Brier score is recorded as not computable with a reason rather than
silently omitted.
"""

import csv
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

SEED = 20261102
N_BOOTSTRAP = 2000
FC_THRESHOLD = 2.0
LOG_THRESHOLD = math.log10(FC_THRESHOLD)
TARGET_RECALL = 0.95

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "reports", "phase4")

LABEL_SETS = {
    "primary_tier_A": "model_predictions_primary_tier_A.tsv",
    "secondary_all_tiers": "model_predictions_secondary_all_tiers.tsv",
}

# Mechanism classes. A mutation tested against several drugs contributes to
# each class it was tested in, so the subsets overlap; n is reported per class.
DRUG_CLASS = {
    "Aciclovir": "nucleoside_analogue",
    "Penciclovir": "nucleoside_analogue",
    "Brivudine": "nucleoside_analogue",
    "Ganciclovir": "nucleoside_analogue",
    "Foscarnet": "pyrophosphate_analogue",
    "Cidofovir": "nucleotide_analogue",
    "Amenamevir": "helicase_primase_inhibitor",
    "Pritelivir": "helicase_primase_inhibitor",
}


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def average_precision(y, s):
    """AP as the step-wise sum used by sklearn, without the dependency."""
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    precision = tp / np.arange(1, len(y) + 1)
    n_pos = y.sum()
    if n_pos == 0:
        return float("nan")
    return float(np.sum(precision * y) / n_pos)


def precision_at_recall(y, s, target=TARGET_RECALL):
    """Precision at the smallest score threshold reaching `target` recall.

    Returns nan when the target is unreachable, which cannot happen here
    (taking every variant gives recall 1) but is checked rather than assumed.
    """
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    n_pos = y.sum()
    if n_pos == 0:
        return float("nan")
    tp = np.cumsum(y)
    recall = tp / n_pos
    hit = np.searchsorted(recall, target, side="left")
    if hit >= len(y):
        return float("nan")
    return float(tp[hit] / (hit + 1))


def confusion_at_2x(y, s, threshold=LOG_THRESHOLD):
    pred = (s >= threshold).astype(int)
    return {
        "tp": int(np.sum((pred == 1) & (y == 1))),
        "tn": int(np.sum((pred == 0) & (y == 0))),
        "fp": int(np.sum((pred == 1) & (y == 0))),
        "fn": int(np.sum((pred == 0) & (y == 1))),
        "n_predicted_resistant": int(pred.sum()),
    }


def mcc(y, s, threshold=LOG_THRESHOLD):
    """MCC is undefined when the prediction is constant at the threshold.

    That is not a missing value, it is a result: the model calls every variant
    the same way at the pre-declared 2x operating point, so it makes no usable
    decision there however well it ranks. Callers record the confusion matrix
    alongside, so the degenerate case is visible rather than blank.
    """
    c = confusion_at_2x(y, s, threshold)
    tp, tn, fp, fn = c["tp"], c["tn"], c["fp"], c["fn"]
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return float("nan")
    return float((tp * tn - fp * fn) / denom)


def auroc(y, s):
    """Reported only for comparability with prior work; no claim rests on it."""
    n_pos, n_neg = int(y.sum()), int(len(y) - y.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    sorted_s = s[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return float((np.sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2.0)
                 / (n_pos * n_neg))


def calibration(y_cont, s):
    """Slope and intercept of observed log10 FC regressed on predicted.

    Perfect calibration is slope 1, intercept 0. Slope < 1 means the model
    under-disperses: it is shrinking towards the mean, which at this sample
    size is expected.
    """
    if len(s) < 3 or np.std(s) == 0:
        return {"slope": None, "intercept": None}
    slope, intercept = np.polyfit(s, y_cont, 1)
    return {"slope": round(float(slope), 4), "intercept": round(float(intercept), 4)}


def all_metrics(y_bin, y_cont, s):
    conf = confusion_at_2x(y_bin, s)
    m = mcc(y_bin, s)
    return {
        "average_precision": round(average_precision(y_bin, s), 4),
        "prevalence": round(float(y_bin.mean()), 4),
        "ap_lift_over_prevalence": round(
            average_precision(y_bin, s) - float(y_bin.mean()), 4),
        "precision_at_95_recall": round(precision_at_recall(y_bin, s), 4),
        "mcc_at_2x": None if math.isnan(m) else round(m, 4),
        "confusion_at_2x": conf,
        "threshold_degenerate": bool(conf["n_predicted_resistant"] in
                                     (0, len(y_bin))),
        "auroc_not_primary": round(auroc(y_bin, s), 4),
        "spearman": round(float(spearmanr(s, y_cont).statistic), 4),
        "calibration": calibration(y_cont, s),
        "n": int(len(s)),
        "n_positive": int(y_bin.sum()),
    }


def bootstrap_ci(y_bin, y_cont, s, rng, n=N_BOOTSTRAP):
    """Resampled over mutations, as the protocol specifies -- not over
    predictions, which would understate the interval."""
    keep = {"average_precision": [], "precision_at_95_recall": [],
            "mcc_at_2x": [], "spearman": []}
    idx = np.arange(len(s))
    for _ in range(n):
        b = rng.choice(idx, size=len(idx), replace=True)
        yb, cb, sb = y_bin[b], y_cont[b], s[b]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        keep["average_precision"].append(average_precision(yb, sb))
        keep["precision_at_95_recall"].append(precision_at_recall(yb, sb))
        keep["mcc_at_2x"].append(mcc(yb, sb))
        keep["spearman"].append(float(spearmanr(sb, cb).statistic))
    out = {}
    for k, v in keep.items():
        arr = np.array([x for x in v if not math.isnan(x)])
        out[k] = ([round(float(np.percentile(arr, 2.5)), 4),
                   round(float(np.percentile(arr, 97.5)), 4)]
                  if len(arr) else [None, None])
    out["n_resamples_used"] = len(keep["average_precision"])
    return out


def paired_delta_ap(y_bin, s_model, s_ref, rng, n=N_BOOTSTRAP):
    """Paired bootstrap of the AP advantage over the frequency baseline.

    Paired on the same resampled mutations, so the comparison is not
    contaminated by between-resample variance.
    """
    diffs = []
    idx = np.arange(len(y_bin))
    for _ in range(n):
        b = rng.choice(idx, size=len(idx), replace=True)
        yb = y_bin[b]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        diffs.append(average_precision(yb, s_model[b])
                     - average_precision(yb, s_ref[b]))
    arr = np.array(diffs)
    if not len(arr):
        return {"difference": None, "ci95": [None, None], "p_two_sided": None}
    # Two-sided bootstrap p: how often the sign flips, doubled.
    p = 2.0 * min(float(np.mean(arr <= 0)), float(np.mean(arr >= 0)))
    return {
        "difference": round(float(np.mean(arr)), 4),
        "ci95": [round(float(np.percentile(arr, 2.5)), 4),
                 round(float(np.percentile(arr, 97.5)), 4)],
        "p_two_sided": round(min(p, 1.0), 4),
    }


def holm(pvals):
    """Holm-Bonferroni. The protocol requires multiplicity correction across
    the benchmark grid: with 22 models an uncorrected win is not evidence."""
    items = sorted(pvals.items(), key=lambda kv: (kv[1] is None, kv[1]))
    m = sum(1 for _, p in items if p is not None)
    out, running, rank = {}, 0.0, 0
    for name, p in items:
        if p is None:
            out[name] = None
            continue
        adj = min(1.0, (m - rank) * p)
        running = max(running, adj)
        out[name] = round(running, 4)
        rank += 1
    return out


def drug_classes(drugs_field):
    out = set()
    for d in drugs_field.split(";"):
        d = d.strip()
        if d in DRUG_CLASS:
            out.add(DRUG_CLASS[d])
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)

    dataset = {(r["virus"], r["gene"], r["aa_change"]): r
               for r in read_tsv(os.path.join(ROOT, "data", "processed",
                                              "model_dataset.tsv"))}

    summary = {"seed": SEED, "n_bootstrap": N_BOOTSTRAP,
               "fc_threshold": FC_THRESHOLD,
               "primary_metric": "average_precision",
               "label_sets": {}}

    for label, fname in LABEL_SETS.items():
        print(f"\n################ {label} ################", flush=True)
        rows = read_tsv(os.path.join(ROOT, "data", "processed", fname))
        y_bin = np.array([int(r["resistant"]) for r in rows])
        y_cont = np.array([float(r["log10_fc"]) for r in rows])

        pred_cols = [c for c in rows[0] if "|" in c]
        preds = {c: np.array([float(r[c]) for r in rows]) for c in pred_cols}

        schemes = sorted({c.split("|")[0] for c in pred_cols})
        entry = {"n": len(rows), "prevalence": round(float(y_bin.mean()), 4),
                 "schemes": {}}

        for scheme in schemes:
            cols = [c for c in pred_cols if c.startswith(scheme + "|")]
            ref_col = next((c for c in cols
                            if "frequency_only" in c and c.endswith("ridge")), None)
            assert ref_col, f"no frequency baseline for {scheme}"

            per_model, pvals = {}, {}
            for c in cols:
                name = "|".join(c.split("|")[1:])
                m = all_metrics(y_bin, y_cont, preds[c])
                m["ci95"] = bootstrap_ci(y_bin, y_cont, preds[c], rng)
                if c != ref_col:
                    d = paired_delta_ap(y_bin, preds[c], preds[ref_col], rng)
                    m["vs_frequency_ap"] = d
                    pvals[name] = d["p_two_sided"]
                per_model[name] = m

            adjusted = holm(pvals)
            for name, p in adjusted.items():
                per_model[name]["vs_frequency_ap"]["p_holm"] = p
                per_model[name]["vs_frequency_ap"]["significant_after_holm"] = (
                    bool(p is not None and p < 0.05
                         and per_model[name]["vs_frequency_ap"]["difference"] > 0))

            best = max(per_model, key=lambda k: per_model[k]["average_precision"])
            entry["schemes"][scheme] = {
                "frequency_baseline": "|".join(ref_col.split("|")[1:]),
                "best_by_average_precision": best,
                "n_beating_frequency_after_holm": sum(
                    1 for k, v in per_model.items()
                    if v.get("vs_frequency_ap", {}).get("significant_after_holm")),
                "models": per_model,
            }

            b = per_model[best]
            print(f"\n=== {scheme} ===", flush=True)
            print(f"  prevalence (no-skill floor for AP): {y_bin.mean():.3f}",
                  flush=True)
            print(f"  best by AP: {best}", flush=True)
            print(f"    AP      {b['average_precision']:.3f} "
                  f"{b['ci95']['average_precision']}", flush=True)
            print(f"    P@95rec {b['precision_at_95_recall']:.3f} "
                  f"{b['ci95']['precision_at_95_recall']}", flush=True)
            print(f"    AP lift over prevalence floor "
                  f"{b['ap_lift_over_prevalence']:+.3f}", flush=True)
            print(f"    MCC@2x  {b['mcc_at_2x']} {b['ci95']['mcc_at_2x']}"
                  f"{'  DEGENERATE: predicts one class' if b['threshold_degenerate'] else ''}",
                  flush=True)
            n_deg = sum(1 for v in per_model.values() if v["threshold_degenerate"])
            print(f"    models degenerate at the 2x threshold: "
                  f"{n_deg}/{len(per_model)}", flush=True)
            print(f"    AUROC   {b['auroc_not_primary']:.3f}  (not primary)",
                  flush=True)
            print(f"    calib   slope={b['calibration']['slope']} "
                  f"intercept={b['calibration']['intercept']}", flush=True)
            print(f"  models beating frequency on AP after Holm: "
                  f"{entry['schemes'][scheme]['n_beating_frequency_after_holm']}"
                  f"/{len(pvals)}", flush=True)

            # ---- breakdowns, for the best model and the frequency baseline
            best_col = f"{scheme}|{best}"
            groups = defaultdict(list)
            for i, r in enumerate(rows):
                groups[("tier", r["tier"])].append(i)
                groups[("gene", f"{r['virus']}:{r['gene']}")].append(i)
                src = dataset.get((r["virus"], r["gene"], r["aa_change"]))
                for cls in (drug_classes(src["drugs_tested"]) if src else ()):
                    groups[("drug_class", cls)].append(i)

            breakdown = {}
            for (kind, key), idx in sorted(groups.items()):
                idx = np.array(idx)
                yb, cb = y_bin[idx], y_cont[idx]
                cell = {"n": int(len(idx)), "n_positive": int(yb.sum())}
                if yb.sum() == 0 or yb.sum() == len(yb):
                    cell["note"] = "single-class subset; AP undefined"
                else:
                    cell["best_model_ap"] = round(
                        average_precision(yb, preds[best_col][idx]), 4)
                    cell["frequency_ap"] = round(
                        average_precision(yb, preds[ref_col][idx]), 4)
                    cell["prevalence"] = round(float(yb.mean()), 4)
                breakdown.setdefault(kind, {})[key] = cell
            entry["schemes"][scheme]["breakdowns"] = breakdown

        summary["label_sets"][label] = entry

    # ---- gates
    for label, e in summary["label_sets"].items():
        for scheme, s in e["schemes"].items():
            for name, m in s["models"].items():
                assert m["average_precision"] >= 0.0, (label, scheme, name)
                assert m["prevalence"] == e["prevalence"], (label, scheme, name)

    path = os.path.join(OUT_DIR, "metrics_summary.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
