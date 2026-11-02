#!/usr/bin/env python
"""Phase 4, step 8 - the manuscript figures.

Three figures, chosen because each carries a claim that a table states less
well:

  Figure 1  The catalogue against the population. Panel A is the frequency
            distribution of all 644 catalogued substitutions; panel B is the
            joint distribution of frequency and measured fold change, which is
            where the project's premise either shows or does not. The empty
            upper-right region - common and strongly resistant - is the
            finding, and no table makes an absence visible.

  Figure 2  pN/pS by gene with bootstrap intervals. The Aim 2 contrast is a
            comparison of intervals, and intervals belong on an axis.

  Figure 3  Where the structural signal lives (panel A) and whether it holds
            prospectively (panel B). Panel A exists to bound the Aim 4 claim
            rather than to advertise it: two of five genes fall to the
            prevalence floor.

Nothing is drawn that is not already in a machine-readable summary or a
release table, so every figure is reproducible from the committed data and
regenerating them cannot introduce a number that no gate checks.

Aim 3 deliberately has no figure. Its result is four counts, three of them
zero, and a bar chart of zeros communicates less than the sentence does.

Outputs, as both vector and raster at 300 dpi:
  reports/phase4/figures/figure1_catalogue_vs_population.{pdf,png}
  reports/phase4/figures/figure2_selection.{pdf,png}
  reports/phase4/figures/figure3_structural_signal.{pdf,png}
  reports/phase4/figures/figure_data.json   every plotted value, flat
"""

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402
import numpy as np                                               # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
REL = ROOT / "data" / "release"
REP = ROOT / "reports"
FIG = REP / "phase4" / "figures"

SEED = 20261102
FC_THRESHOLD = 2.0
# The fold change above which no common variant is observed at all. Established
# from the data in this script rather than assumed, and checked by p4_05_verify,
# because the figure's central claim is an absence and an absence is only as
# good as the boundary drawn around it.
EMPTY_FC = 5.0

# Okabe-Ito, which stays distinguishable under the common forms of colour
# vision deficiency and survives greyscale printing in rank order.
C = {
    "structure": "#0072B2",
    "frequency": "#D55E00",
    "conservation": "#009E73",
    "esm": "#CC79A7",
    "full": "#E69F00",
    "floor": "#666666",
    "ink": "#1A1A1A",
    "pale": "#BBBBBB",
}

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.linewidth": 0.6,
    "axes.edgecolor": C["ink"],
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.frameon": False,
    "legend.fontsize": 7,
    "lines.linewidth": 1.2,
})

plotted = {}


def rows(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def num(x):
    """Return float(x), or None for the empty strings the tables use."""
    if x in ("", "NA", None):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def panel_label(ax, text):
    ax.text(-0.14, 1.04, text, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left", color=C["ink"])


# --------------------------------------------------------------- figure 1
def figure1():
    rc = rows(REL / "reclassification_table.tsv")
    mv = rows(REL / "master_variant_table.tsv")

    fig = plt.figure(figsize=(7.4, 3.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.40)
    axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # ---- panel A: the 644 by frequency band
    order = [("not_observed", "Not observed"), ("singleton", "Singleton"),
             ("rare", "Rare (<1%)"), ("uncommon_polymorphism", "Uncommon (1-10%)"),
             ("common_polymorphism", "Common (10-50%)"),
             ("majority_allele", "Majority (>50%)")]
    counts = {}
    for r in rc:
        counts[r["frequency_band"]] = counts.get(r["frequency_band"], 0) + 1
    labels = [lab for _, lab in order]
    vals = [counts.get(key, 0) for key, _ in order]
    # The two bands that carry the paper's claim are the two that are coloured;
    # shading everything would say that every band matters equally.
    colours = [C["pale"]] * 4 + [C["frequency"], C["structure"]]

    y = np.arange(len(vals))
    axA.barh(y, vals, color=colours, height=0.68, edgecolor="none")
    axA.set_yticks(y, labels)
    axA.invert_yaxis()
    axA.set_xlabel("Catalogued substitutions (n = %d)" % len(rc))
    axA.set_xlim(0, max(vals) * 1.17)
    for yi, v in zip(y, vals):
        axA.text(v + max(vals) * 0.015, yi, str(v), va="center", fontsize=7,
                 color=C["ink"])
    axA.set_title("Catalogued mutations by frequency", loc="left")
    panel_label(axA, "A")
    plotted["figure1A_bands"] = dict(zip(labels, vals))

    # ---- panel B: frequency against measured fold change
    phen = [r for r in mv if num(r["max_fold_change"]) is not None]
    obs = [r for r in phen if r["observed"] == "1"]
    unobs = [r for r in phen if r["observed"] != "1"]

    rng = np.random.default_rng(SEED)
    # Unobserved variants have no frequency to plot, and dropping them would
    # hide that they are the overwhelming majority of the phenotyped set. They
    # go in a jittered strip left of the axis break instead.
    strip_x = 10 ** (rng.normal(-3.08, 0.055, len(unobs)))
    strip_y = [num(r["max_fold_change"]) for r in unobs]
    axB.scatter(strip_x, strip_y, s=9, color=C["pale"], edgecolor="none",
                alpha=0.75, zorder=2)

    ox = [num(r["frequency"]) for r in obs]
    oy = [num(r["max_fold_change"]) for r in obs]
    axB.scatter(ox, oy, s=20, color=C["structure"], edgecolor="white",
                linewidth=0.4, zorder=4)

    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(1.4e-3, 1.4)
    axB.set_ylim(0.2, 1.2e4)

    axB.axhline(FC_THRESHOLD, color=C["ink"], lw=0.7, ls="--", zorder=1)
    axB.text(1.6e-3, FC_THRESHOLD * 1.18, "declared 2x threshold", ha="left",
             fontsize=6.3, color=C["ink"])
    axB.axvline(0.05, color=C["floor"], lw=0.6, ls=":", zorder=1)
    axB.text(0.0465, 6.5e3, "5%", fontsize=6.3, color=C["floor"], ha="right")

    # The claim the figure exists to make is an absence, so the boundary has to
    # be the one that is actually empty. Above 5% frequency nothing reaches 5x.
    # The two entries that do clear the declared 2x threshold at that frequency
    # are both marginal (2.3x and 2.4x), which is why the shaded region starts
    # at 5x rather than at the threshold.
    axB.add_patch(plt.Rectangle((0.05, EMPTY_FC), 1.4 - 0.05, 1.2e4,
                                facecolor=C["frequency"], alpha=0.07,
                                edgecolor="none", zorder=0))
    axB.axhline(EMPTY_FC, color=C["frequency"], lw=0.6, ls=":", zorder=1)
    axB.text(0.235, 9e2, "common (>5%) and\nresistant above 5x\nn = 0",
             ha="center", va="center", fontsize=6.5, color=C["frequency"],
             linespacing=1.35)

    # Label the variants the text discusses, and only those. R222C is labelled
    # leftwards on purpose: it sits at 4.6%, just outside the shaded region, and
    # a label running rightwards would read as though it were inside.
    # R41H and P920S both sit within a few points of the 2x line, so they are
    # pushed clear of it in opposite directions rather than labelled in place:
    # R41H above, where the gap to the 5x line is widest, and P920S below.
    place = {"R176Q": ((7, -1), "left"), "R222C": ((-7, -1), "right"),
             "S724N": ((7, -2), "left"), "R41H": ((-6, 9), "right"),
             "G39E": ((-7, 5), "right"), "P920S": ((2, -11), "center")}
    for r in obs:
        spec = place.get(r["aa_change"])
        if spec:
            (dx, dy), ha = spec
            # A label pushed clear of a crowded line is far enough from its
            # point to be read as belonging to a neighbour, so those two get
            # a hairline leader. The rest sit against their point and do not.
            arrow = (dict(arrowstyle="-", lw=0.4, color=C["floor"],
                          shrinkA=0.5, shrinkB=2.5)
                     if abs(dy) >= 9 else None)
            axB.annotate("%s %s" % (r["gene"], r["aa_change"]),
                         (num(r["frequency"]), num(r["max_fold_change"])),
                         textcoords="offset points", xytext=(dx, dy),
                         fontsize=6.2, color=C["ink"], ha=ha, zorder=5,
                         arrowprops=arrow)

    axB.set_xlabel("Frequency in circulating virus (log)")
    axB.set_ylabel("Maximum measured fold change (log)")
    axB.set_title("Frequency against measured resistance", loc="left")
    ticks = [1e-2, 1e-1, 1]
    axB.set_xticks([10 ** -3.08] + ticks)
    axB.set_xticklabels(["not\nobserved", "1%", "10%", "100%"])
    axB.get_xticklabels()[0].set_color(C["floor"])
    axB.get_xticklabels()[0].set_fontsize(6.2)

    # Below the axis, because every in-axes corner holds data.
    axB.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=4, color=C["structure"],
               label="Observed in circulation (n = %d)" % len(obs)),
        Line2D([], [], marker="o", ls="", ms=3.4, color=C["pale"],
               label="Not observed (n = %d)" % len(unobs)),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=2,
        columnspacing=1.4, handletextpad=0.35)
    panel_label(axB, "B")

    med_o = float(np.median(oy))
    med_u = float(np.median(strip_y))
    plotted["figure1B"] = {
        "n_phenotyped": len(phen), "n_observed": len(obs),
        "n_not_observed": len(unobs),
        "median_fold_change_observed": med_o,
        "median_fold_change_not_observed": med_u,
        "n_common_at_2x": sum(
            1 for r in obs
            if num(r["frequency"]) >= 0.05
            and num(r["max_fold_change"]) >= FC_THRESHOLD),
        "n_common_above_5x": sum(
            1 for r in obs
            if num(r["frequency"]) >= 0.05
            and num(r["max_fold_change"]) >= EMPTY_FC),
        "max_fold_change_at_or_above_5pc": max(
            num(r["max_fold_change"]) for r in obs
            if num(r["frequency"]) >= 0.05),
    }
    save(fig, "figure1_catalogue_vs_population")


# --------------------------------------------------------------- figure 2
def figure2():
    sel = load(REP / "phase2" / "selection_summary.json")
    pg = sel["per_gene"]
    a2 = sel["aim2_contrast"]["HSV1"]

    # Reliable genes first, HSV-1 above HSV-2, each block by point estimate.
    # The unreliable genes are shown greyed rather than dropped, because a
    # reader is entitled to see the estimates the reliability rule excluded.
    def key(k):
        return (0 if pg[k]["reliable"] else 1, 0 if k.startswith("HSV1") else 1,
                -pg[k]["pN_pS"])
    keys = sorted(pg, key=key)

    fig, ax = plt.subplots(figsize=(6.2, 3.3))
    for i, k in enumerate(keys):
        g = pg[k]
        y = len(keys) - 1 - i
        rel = g["reliable"]
        col = C["structure"] if rel else C["pale"]
        lo, hi = g["ci"]
        ax.plot([lo, hi], [y, y], color=col, lw=1.4, solid_capstyle="round",
                zorder=2)
        ax.plot([g["pN_pS"]], [y], marker="o", ms=5, color=col,
                markeredgecolor="white", markeredgewidth=0.5, zorder=3)
    # n goes in the tick label rather than beside each interval, where it
    # collided with the widest CI and ran off the axis.
    names = {"HSV1:UL23": "HSV-1 UL23  thymidine kinase",
             "HSV1:UL30": "HSV-1 UL30  DNA polymerase",
             "HSV1:UL5": "HSV-1 UL5  helicase",
             "HSV1:UL52": "HSV-1 UL52  primase",
             "HSV2:UL30": "HSV-2 UL30",
             "HSV2:UL23": "HSV-2 UL23  (unreliable)",
             "HSV2:UL5": "HSV-2 UL5  (unreliable)",
             "HSV2:UL52": "HSV-2 UL52  (unreliable)"}
    ax.set_yticks(range(len(keys)),
                  ["%s  (n = %d)" % (names[k], pg[k]["n_isolates"])
                   for k in reversed(keys)])
    ax.set_xlabel("pN/pS against the population consensus (95% bootstrap CI)")
    ax.axvline(1.0, color=C["ink"], lw=0.6, ls="--")
    ax.text(1.03, len(keys) - 0.55, "neutrality", fontsize=6.3, color=C["ink"],
            va="top")
    ax.set_xlim(0, 4.0)
    ax.set_ylim(-0.7, len(keys) - 0.3)

    # The contrast is stated rather than drawn between the two rows: UL5 sits
    # between UL23 and UL30 in rank order, so a connector would cross a gene
    # the comparison has nothing to do with.
    ax.text(1.45, len(keys) - 1.6,
            "HSV-1 UL23 − UL30\ndifference %+.4f\n[%+.4f, %+.4f]\n"
            "UL23 higher in 100%% of resamples"
            % (a2["difference"], a2["ci_low"], a2["ci_high"]),
            fontsize=6.8, color=C["frequency"], va="center", linespacing=1.5)

    ax.set_title("Purifying selection by gene", loc="left")
    plotted["figure2"] = {
        k: {"pN_pS": pg[k]["pN_pS"], "ci": pg[k]["ci"],
            "reliable": pg[k]["reliable"], "n_isolates": pg[k]["n_isolates"]}
        for k in keys}
    plotted["figure2"]["aim2_contrast"] = {
        "difference": a2["difference"],
        "ci": [a2["ci_low"], a2["ci_high"]],
        "fraction_resamples": a2["fraction_of_resamples_UL23_higher"]}
    save(fig, "figure2_selection")


# --------------------------------------------------------------- figure 3
def figure3():
    met = load(REP / "phase4" / "metrics_summary.json")
    tem = load(REP / "phase4" / "temporal_holdout.json")
    pb = met["label_sets"]["primary_tier_A"]["schemes"]["position_blocks"]
    gene = pb["breakdowns"]["gene"]

    fig = plt.figure(figsize=(7.4, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.30)
    axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # ---- panel A: per-gene average precision against the no-skill floor
    keys = [k for k in ("HSV1:UL23", "HSV2:UL23", "HSV1:UL30", "HSV1:UL5",
                        "HSV1:UL52") if "best_model_ap" in gene.get(k, {})]
    x = np.arange(len(keys))
    w = 0.36
    struct = [gene[k]["best_model_ap"] for k in keys]
    freq = [gene[k]["frequency_ap"] for k in keys]
    floor = [gene[k]["prevalence"] for k in keys]

    axA.bar(x - w / 2, struct, w, color=C["structure"], label="Structure only",
            edgecolor="none", zorder=2)
    axA.bar(x + w / 2, freq, w, color=C["frequency"], label="Frequency only",
            edgecolor="none", zorder=2)
    for xi, f in zip(x, floor):
        axA.plot([xi - w - 0.06, xi + w + 0.06], [f, f], color=C["ink"],
                 lw=1.0, ls=(0, (2.5, 1.6)), zorder=3)

    # Which genes have a ligand-resolved site is the whole point of the panel,
    # so it is marked on the axis itself rather than left to the caption. A
    # dagger on the tick keeps it legible where two such genes sit adjacent.
    nosite = {"HSV1:UL5", "HSV1:UL52"}
    axA.set_xticks(x, [k.replace("HSV1:", "HSV-1\n").replace("HSV2:", "HSV-2\n")
                       + (" †" if k in nosite else "") for k in keys])
    axA.text(0.0, -0.30, "† no ligand-resolved active site",
             transform=axA.transAxes, fontsize=6.3, color=C["floor"])
    axA.set_ylabel("Average precision")
    axA.set_ylim(0, 1.08)
    axA.set_title("Where the structural signal lives", loc="left")
    for i, k in enumerate(keys):
        axA.annotate("n = %d" % gene[k]["n"],
                     (x[i], max(struct[i], freq[i], floor[i]) + 0.03),
                     ha="center", fontsize=6.2, color=C["floor"])
    axA.legend(handles=[
        Line2D([], [], marker="s", ls="", ms=5, color=C["structure"],
               label="Structure only"),
        Line2D([], [], marker="s", ls="", ms=5, color=C["frequency"],
               label="Frequency only"),
        Line2D([], [], ls=(0, (2.5, 1.6)), lw=1.0, color=C["ink"],
               label="No-skill floor"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3,
        columnspacing=1.1, handletextpad=0.5)
    panel_label(axA, "A")
    plotted["figure3A"] = {k: {"n": gene[k]["n"], "prevalence": gene[k]["prevalence"],
                               "structure_ap": gene[k]["best_model_ap"],
                               "frequency_ap": gene[k]["frequency_ap"]}
                           for k in keys}

    # ---- panel B: the temporal hold-out
    tp = tem["label_sets"]["primary_tier_A"]["cutoffs"]
    cutoffs = sorted(tp, key=int)
    series = [("structure", "Structure", C["structure"], "o"),
              ("conservation_only", "Conservation", C["conservation"], "s"),
              ("esm_only", "ESM-2", C["esm"], "^"),
              ("full_model", "Full model", C["full"], "D"),
              ("frequency_only", "Frequency", C["frequency"], "v")]
    plotted["figure3B"] = {}
    for key, label, col, mk in series:
        ys = []
        for co in cutoffs:
            m = tp[co]["models"]
            if key == "structure":
                # The pre-declared comparison is the structure feature group,
                # whose better estimator differs by cut-off; taking the worse
                # one would understate it and cherry-picking per model would
                # not be the group's performance.
                v = max(m["structure_only"]["ap_lift_over_prevalence"],
                        m["structure_only_forest"]["ap_lift_over_prevalence"])
            else:
                v = m[key]["ap_lift_over_prevalence"]
            ys.append(v)
        lw = 1.8 if key == "structure" else 1.0
        z = 4 if key == "structure" else 2
        axB.plot([int(c) for c in cutoffs], ys, marker=mk, ms=3.6, color=col,
                 lw=lw, label=label, zorder=z)
        plotted["figure3B"][key] = dict(zip(cutoffs, [round(v, 4) for v in ys]))

    axB.axhline(0, color=C["ink"], lw=0.7, ls="--", zorder=1)
    axB.text(2015.2, 0.006, "no skill", fontsize=6.3, color=C["ink"],
             ha="right")
    axB.set_xticks([int(c) for c in cutoffs],
                   ["%s\n%d→%d" % (c, tp[c]["n_train"], tp[c]["n_test"])
                    for c in cutoffs])
    axB.set_xlabel("Training cut-off year, and train → test split")
    axB.set_ylabel("Average precision, lift over test-set floor")
    axB.set_title("Prospective validation across five cut-offs", loc="left")
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3,
               columnspacing=1.1, handletextpad=0.5)
    axB.set_ylim(-0.13, 0.27)
    panel_label(axB, "B")

    save(fig, "figure3_structural_signal")


def save(fig, stem):
    """Write both formats, byte-reproducibly.

    A PDF carries a CreationDate by default, so re-running this script would
    rewrite all three PDFs with identical content and a new timestamp, dirty
    the working tree and break any checksum taken over them. Passing the date
    as None omits it. The PNGs need no equivalent: matplotlib writes no
    timestamp into them, which is why a re-run leaves them untouched.
    """
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        meta = {"CreationDate": None} if ext == "pdf" else {}
        fig.savefig(FIG / f"{stem}.{ext}", metadata=meta)
    plt.close(fig)
    print(f"  {stem}.pdf, {stem}.png")


def main():
    print("figures ->", FIG)
    figure1()
    figure2()
    figure3()
    (FIG / "figure_data.json").write_text(
        json.dumps(plotted, indent=2, sort_keys=True), encoding="utf-8")
    print("  figure_data.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
