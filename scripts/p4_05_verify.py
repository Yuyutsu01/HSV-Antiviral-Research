"""Phase 4, step 5 — check the write-up against the outputs.

Every number quoted in MANUSCRIPT.md, README.md, PRE_DECLARATION.md and
DATA_DICTIONARY.md that came from a machine-readable summary is re-read from
that summary here and compared. A report can drift from its data silently;
nothing else in the pipeline would catch a mistyped figure in prose.

This is a gate, not a report: it exits non-zero on the first mismatch.
"""

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REL = ROOT / "data" / "release"
REP = ROOT / "reports"

failures = []
checks = 0


def load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def rows(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check(label, got, want, tol=0.0):
    global checks
    checks += 1
    ok = (abs(got - want) <= tol) if isinstance(want, (int, float)) \
        else got == want
    if not ok:
        failures.append(f"{label}: report says {want!r}, data says {got!r}")


def main():
    global checks
    metrics = load(REP / "phase4" / "metrics_summary.json")
    gap = load(REP / "phase4" / "circularity_gap.json")
    cross = load(REP / "phase4" / "cross_virus.json")
    release = load(REP / "phase4" / "release_summary.json")
    audit = load(REP / "phase3" / "audit_summary.json")
    selection = load(REP / "phase2" / "selection_summary.json")
    hp = load(REP / "phase2" / "hp_baseline_summary.json")

    # ---------------------------------------------------------- Aim 1 counts
    rc = rows(REL / "reclassification_table.tsv")
    check("catalogued substitutions", len(rc), 644)
    check("observed in circulation", sum(1 for r in rc if r["observed"] == "1"), 128)
    check("not observed", sum(1 for r in rc if r["observed"] == "0"), 516)

    bands = {}
    flags = {}
    for r in rc:
        bands[r["frequency_band"]] = bands.get(r["frequency_band"], 0) + 1
        flags[r["flag"]] = flags.get(r["flag"], 0) + 1
    check("majority alleles", bands.get("majority_allele"), 11)
    check("common polymorphisms", bands.get("common_polymorphism"), 15)
    check("uncommon polymorphisms", bands.get("uncommon_polymorphism"), 36)
    check("singletons", bands.get("singleton"), 45)
    check("rare", bands.get("rare"), 21)
    check("reference_artifact", flags.get("reference_artifact"), 10)
    check("reclassify_as_polymorphism", flags.get("reclassify_as_polymorphism"), 7)
    check("unphenotyped_but_common", flags.get("unphenotyped_but_common"), 41)
    check("consistent_rare", flags.get("consistent_rare"), 60)
    check("consistent_with_resistance", flags.get("consistent_with_resistance"), 6)
    check("circulating_resistance_allele",
          flags.get("circulating_resistance_allele"), 3)
    check("study_dominated among observed",
          sum(1 for r in rc if r["study_dominated"] == "1"), 60)

    # the three alleles, as tabulated in the manuscript and data dictionary
    want = {"R222C": (0.046, 11, 3, 0.7273, 0.013393, "98"),
            "R176Q": (0.017, 4, 1, 1.0, 0.0, "200"),
            "R41H": (0.051, 12, 5, 0.3333, 0.038095, "2.3")}
    for r in rc:
        if r["gene"] == "UL23" and r["virus"] == "HSV1" \
                and r["aa_change"] in want:
            f, n, ns, share, loso, fc = want[r["aa_change"]]
            a = r["aa_change"]
            check(f"{a} frequency", float(r["frequency"]), f, 0.001)
            check(f"{a} carriers", int(r["count"]), n)
            check(f"{a} studies", int(r["n_studies_carrying"]), ns)
            check(f"{a} max_study_share", float(r["max_study_share"]), share, 0.001)
            check(f"{a} loso_min", float(r["loso_min_freq"]), loso, 0.0001)
            check(f"{a} fold change", r["max_fold_change"], fc)

    # ------------------------------------------------------------ Aims 2 & 3
    a2 = selection["aim2_contrast"]["HSV1"]
    check("Aim 2 difference", a2["difference"], 0.1980, 0.0001)
    check("Aim 2 CI low", a2["ci_low"], 0.1664, 0.0001)
    check("Aim 2 CI high", a2["ci_high"], 0.2361, 0.0001)
    check("Aim 2 bootstrap support", a2["fraction_of_resamples_UL23_higher"], 1.0)
    pg = selection["per_gene"]
    check("UL23 pN/pS", pg["HSV1:UL23"]["pN_pS"], 0.371, 0.001)
    check("UL30 pN/pS", pg["HSV1:UL30"]["pN_pS"], 0.173, 0.001)
    check("UL23 Sd", pg["HSV1:UL23"]["Sd"], 620.0)
    check("UL30 Sd", pg["HSV1:UL30"]["Sd"], 1535.5)
    check("HSV2 UL23 unreliable", pg["HSV2:UL23"]["reliable"], False)

    bc = hp["baseline_carriage"]
    check("HSV1 UL5 carriers", bc["HSV1:UL5"]["n_carrying_catalogued_substitution"], 0)
    check("HSV1 UL5 isolates", bc["HSV1:UL5"]["n_isolates"], 251)
    check("HSV1 UL52 carriers", bc["HSV1:UL52"]["n_carrying_catalogued_substitution"], 0)
    check("HSV1 UL52 isolates", bc["HSV1:UL52"]["n_isolates"], 299)
    check("HSV2 UL5 carriers", bc["HSV2:UL5"]["n_carrying_catalogued_substitution"], 2)
    check("HSV2 UL5 isolates", bc["HSV2:UL5"]["n_isolates"], 247)
    check("HSV2 UL52 isolates", bc["HSV2:UL52"]["n_isolates"], 225)
    ct = hp["constraint_test"]
    check("HSV1 UL5 constraint p", ct["HSV1:UL5"]["p_more_conserved_than_random"],
          0.54, 0.01)
    check("HSV1 UL52 constraint p", ct["HSV1:UL52"]["p_more_conserved_than_random"],
          0.375, 0.001)
    # Reported as non-significant: the manuscript's claim is the side of 0.05,
    # not the value, so the check is the inequality as well as the figure.
    check("HSV2 UL5 constraint p", ct["HSV2:UL5"]["p_more_conserved_than_random"],
          0.052, 0.001)
    if ct["HSV2:UL5"]["p_more_conserved_than_random"] < 0.05:
        failures.append(
            "manuscript reports HSV-2 UL5 constraint as non-significant, but "
            f"p = {ct['HSV2:UL5']['p_more_conserved_than_random']}")
    check("HSV2 UL5 variability ratio", ct["HSV2:UL5"]["ratio"], 2.21, 0.001)

    # -------------------------------------------------------------- Aim 4
    pb = metrics["label_sets"]["primary_tier_A"]["schemes"]["position_blocks"]
    check("Tier A n", metrics["label_sets"]["primary_tier_A"]["n"], 119)
    check("Tier A prevalence", metrics["label_sets"]["primary_tier_A"]["prevalence"],
          0.731, 0.001)
    # The abstract, §5.2, §7 and §8 all turn on this one count, and nothing
    # checked it until a re-run moved it from 13 to 14 silently.
    dataset = load(REP / "phase3" / "dataset_summary.json")
    check("Tier A variants observed in circulation",
          dataset["primary_tier_A"]["freq_observed"], 14)
    m = pb["models"]
    check("structure forest AP", m["baseline:structure_only|forest"]["average_precision"],
          0.916, 0.001)
    check("structure forest lift",
          m["baseline:structure_only|forest"]["ap_lift_over_prevalence"], 0.185, 0.001)
    check("structure ridge AP", m["baseline:structure_only|ridge"]["average_precision"],
          0.905, 0.001)
    check("frequency AP", m["baseline:frequency_only|ridge"]["average_precision"],
          0.723, 0.001)
    check("frequency lift below floor",
          m["baseline:frequency_only|ridge"]["ap_lift_over_prevalence"], -0.008, 0.001)
    check("conservation AP", m["baseline:conservation_only|ridge"]["average_precision"],
          0.717, 0.001)
    check("conservation lift below floor",
          m["baseline:conservation_only|ridge"]["ap_lift_over_prevalence"],
          -0.014, 0.001)
    check("frequency P@95",
          m["baseline:frequency_only|ridge"]["precision_at_95_recall"], 0.762, 0.001)
    check("esm AP", m["baseline:esm_only|ridge"]["average_precision"], 0.848, 0.001)
    check("structure forest P@95",
          m["baseline:structure_only|forest"]["precision_at_95_recall"], 0.790, 0.001)
    check("structure forest MCC",
          m["baseline:structure_only|forest"]["mcc_at_2x"], 0.152, 0.001)
    check("models surviving Holm", pb["n_beating_frequency_after_holm"], 3)
    check("degenerate at 2x",
          sum(1 for v in m.values() if v["threshold_degenerate"]), 8)
    # §5.2 names the three survivors, and §5.4 quotes their calibration slopes.
    for name, slope in (("baseline:structure_only|forest", 0.77),
                        ("baseline:structure_only|ridge", 1.20),
                        ("staged:stage5_plus_esm|forest", 0.95)):
        check(f"{name} calibration slope", m[name]["calibration"]["slope"],
              slope, 0.005)
    check("uncorrected wins over frequency",
          sum(1 for k, v in m.items()
              if k != pb["frequency_baseline"]
              and (v.get("vs_frequency_ap") or {}).get("p_two_sided") is not None
              and v["vs_frequency_ap"]["p_two_sided"] < 0.05), 7)

    gene = pb["breakdowns"]["gene"]
    for key, n, ap in (("HSV1:UL23", 31, 0.988), ("HSV2:UL23", 12, 1.000),
                       ("HSV1:UL30", 45, 0.823), ("HSV1:UL5", 7, 0.567),
                       ("HSV1:UL52", 6, 0.411)):
        check(f"{key} n", gene[key]["n"], n)
        check(f"{key} structure AP", gene[key]["best_model_ap"], ap, 0.001)
    dc = pb["breakdowns"]["drug_class"]
    check("HPI class AP", dc["helicase_primase_inhibitor"]["best_model_ap"],
          0.582, 0.001)
    check("HPI class n", dc["helicase_primase_inhibitor"]["n"], 14)

    # ------------------------------------------------------ circularity audit
    for scheme in ("leave_one_gene_out", "position_blocks"):
        for model, g in gap["schemes"][scheme].items():
            if "note" in g:
                continue
            v = g["average_precision"]["gap_random_minus_matched"]
            # The bound is one decimal place looser than the extreme value it
            # has to contain. Stating it as -0.015 once before let a gap of
            # -0.0152 sit outside a range that its own rounding appeared to
            # satisfy; the prose quotes -0.016 for that reason.
            checks_ok = -0.016 <= v <= 0.003
            if not checks_ok:
                failures.append(
                    f"matched-negative gap {scheme}/{model} = {v}, outside the "
                    f"reported range -0.016..+0.003")
    pairs = {(p["a"], p["b"]): p["spearman"]
             for p in gap["correlation_matrix"]["pairs_above_0.5"]}
    check("conservation vs frequency correlation",
          pairs.get(("cons_site_variability", "freq_frequency")), 0.888, 0.001)
    check("n correlated pairs above 0.5", len(pairs), 7)
    assert not any("struct" in a or "struct" in b for a, b in pairs), \
        "manuscript claims no structural feature correlates above 0.5"

    pbk = audit["position_blocks"]["esm_residualised"]
    check("esm alone, blocks", pbk["esm_alone"]["spearman"], 0.3053, 0.001)
    check("esm alone residualised, blocks",
          pbk["esm_alone_residualised"]["spearman"], -0.0109, 0.001)
    check("esm variance explained by confounds",
          pbk["esm_variance_explained_by_confounds_r2"]["esm_masked_marginal"],
          0.336, 0.001)
    logo = audit["leave_one_gene_out"]["esm_residualised"]
    check("esm alone, logo", logo["esm_alone"]["spearman"], 0.1378, 0.001)
    check("esm alone residualised, logo",
          logo["esm_alone_residualised"]["spearman"], -0.1134, 0.001)

    # §6.1: the two schemes' permutation nulls, which the report quotes as
    # ranges and uses to warn against scoring against zero.
    lg_nulls = [v["null_mean"]
                for v in audit["leave_one_gene_out"]["permutation_tests"].values()]
    check("logo null min", round(min(lg_nulls), 2), -0.29, 0.001)
    check("logo null max", round(max(lg_nulls), 2), -0.16, 0.001)
    pb_nulls = audit["position_blocks"]["permutation_tests"]
    signal = ("structure_only", "structure_only_forest", "esm_only",
              "esm_only_forest", "full_model", "full_model_forest")
    sig_nulls = [pb_nulls[k]["null_mean"] for k in signal]
    check("blocks null min, signal groups", round(min(sig_nulls), 3), -0.014, 0.001)
    check("blocks null max, signal groups", round(max(sig_nulls), 3), 0.079, 0.001)
    check("blocks null, frequency", pb_nulls["frequency_only"]["null_mean"],
          -0.167, 0.001)
    check("blocks null, conservation", pb_nulls["conservation_only"]["null_mean"],
          -0.115, 0.001)

    # ---------------------------------------------------- cross-virus transfer
    cv = cross["label_sets"]["primary_tier_A"]["directions"]
    h2h1 = cv["HSV2_to_HSV1"]
    check("HSV2->HSV1 test n", h2h1["n_test"], 89)
    check("HSV2->HSV1 train n", h2h1["n_train"], 30)
    check("HSV2->HSV1 structure AP",
          h2h1["models"]["structure_only"]["average_precision"], 0.904, 0.001)
    check("HSV2->HSV1 floor",
          h2h1["models"]["structure_only"]["prevalence"], 0.685, 0.001)
    check("HSV2->HSV1 full model AP",
          h2h1["models"]["full_model"]["average_precision"], 0.632, 0.001)

    # -------------------------------------------------------- temporal hold-out
    temporal = load(REP / "phase4" / "temporal_holdout.json")
    check("primary cutoff", temporal["primary_cutoff"], 2012)
    doi = load(ROOT / "data" / "raw" / "doi_years.json")["years"]
    resolved = sum(1 for v in doi.values() if v.get("year"))
    check("DOIs resolved to a year", resolved, 89)
    check("DOIs queried", len(doi), 89)
    years = [v["year"] for v in doi.values() if v.get("year")]
    check("earliest study year", min(years), 1987)
    check("latest study year", max(years), 2025)
    check("publishers represented",
          len({v.get("publisher") for v in doi.values() if v.get("year")}), 10)

    tp = temporal["label_sets"]["primary_tier_A"]
    check("primary variants with a year", tp["n_with_year"], 119)
    check("secondary variants with a year",
          temporal["label_sets"]["secondary_all_tiers"]["n_with_year"], 235)

    # the lift table quoted in the manuscript, primary Tier A
    want_lift = {
        "2005": {"structure_only": 0.195, "frequency_only": -0.077,
                 "conservation_only": 0.085, "esm_only": 0.118,
                 "full_model": 0.110, "_split": (41, 78)},
        "2008": {"structure_only": 0.215, "frequency_only": 0.019,
                 "conservation_only": 0.100, "esm_only": 0.052,
                 "full_model": 0.126, "_split": (47, 72)},
        "2010": {"structure_only": 0.218, "frequency_only": 0.020,
                 "conservation_only": 0.062, "esm_only": -0.011,
                 "full_model": 0.123, "_split": (52, 67)},
        "2012": {"structure_only_forest": 0.229, "frequency_only": 0.030,
                 "conservation_only": 0.042, "esm_only": 0.050,
                 "full_model": 0.090, "_split": (59, 60)},
        "2015": {"structure_only_forest": 0.170, "frequency_only": -0.006,
                 "conservation_only": 0.162, "esm_only": 0.101,
                 "full_model": 0.066, "_split": (70, 49)},
    }
    for cutoff, wants in want_lift.items():
        c = tp["cutoffs"][cutoff]
        tr, te = wants.pop("_split")
        check(f"{cutoff} train n", c["n_train"], tr)
        check(f"{cutoff} test n", c["n_test"], te)
        for model, lift in wants.items():
            check(f"{cutoff} {model} lift",
                  c["models"][model]["ap_lift_over_prevalence"], lift, 0.001)

    # the claims that carry the section
    for cutoff, c in tp["cutoffs"].items():
        best = c["best_by_average_precision"]
        if not best.startswith("structure_only"):
            failures.append(
                f"manuscript claims structure wins every primary cut-off, but "
                f"{cutoff} was won by {best}")
        s = c["models"]["structure_only"]["ap_lift_over_prevalence"]
        sf = c["models"]["structure_only_forest"]["ap_lift_over_prevalence"]
        if max(s, sf) < 0.13:
            failures.append(
                f"manuscript claims structure lift never below +0.13; "
                f"{cutoff} gives {max(s, sf)}")
        # §5.6 no longer claims the alternatives change sign — they do not. It
        # claims instead that none of them reaches +0.163 at any split, which
        # is the comparison the table is actually making.
        others = max(m["ap_lift_over_prevalence"]
                     for name, m in c["models"].items()
                     if not name.startswith("structure_only"))
        if others > 0.163:
            failures.append(
                f"manuscript claims no non-structure lift exceeds +0.163; "
                f"{cutoff} gives {others}")
        if c["train_prevalence"] <= c["test_prevalence"]:
            failures.append(
                f"manuscript claims prevalence falls from train to test; "
                f"{cutoff} has {c['train_prevalence']} -> {c['test_prevalence']}")
    ts = temporal["label_sets"]["secondary_all_tiers"]["cutoffs"]
    # The quoted range spans both structure variants across all cut-offs, so
    # the bounds are over the flattened set, not over a per-cut-off maximum.
    sec = [c["models"][m]["ap_lift_over_prevalence"] for c in ts.values()
           for m in ("structure_only", "structure_only_forest")]
    check("secondary structure lift min", round(min(sec), 3), 0.245, 0.001)
    check("secondary structure lift max", round(max(sec), 3), 0.349, 0.001)
    n_sec_struct = sum(1 for c in ts.values()
                       if c["best_by_average_precision"].startswith("structure"))
    check("secondary cut-offs won by structure", n_sec_struct, 5)
    sec_freq = [c["models"]["frequency_only"]["ap_lift_over_prevalence"]
                for c in ts.values()]
    check("secondary frequency lift min", round(min(sec_freq), 3), 0.001, 0.001)
    check("secondary frequency lift max", round(max(sec_freq), 3), 0.025, 0.001)

    # ----------------------------------------------------------- release stats
    st = release["master_table"]
    check("master table rows", st["n_rows"], 1263)
    check("master in catalogue", st["n_in_catalogue"], 644)
    check("master observed", st["n_observed"], 747)
    check("master both", st["n_both"], 128)
    check("master catalogue only", st["n_catalogue_only"], 516)
    check("master observed only", st["n_observed_only"], 619)
    check("master in model dataset", st["n_in_model_dataset"], 235)
    fates = release["accession_fates"]
    check("accessions analysed", fates["analysed"], 599)
    check("accessions total", sum(fates.values()), 1143)
    check("no cds extracted", fates["excluded_no_cds_extracted"], 168)
    check("downstream filter", fates["excluded_downstream_filter"], 187)
    check("cds qc", fates["excluded_cds_qc"], 127)
    check("lab strain", fates["excluded_laboratory_strain"], 30)
    check("patient duplicate", fates["excluded_patient_duplicate"], 18)
    check("non-human host", fates["excluded_non_human_host"], 14)

    # ---- the variants-of-unknown-significance pool, quoted in the
    # pre-declaration and the manuscript's limitations as the reason the
    # positive-unlabelled extension is not attempted.
    mv = rows(REL / "master_variant_table.tsv")
    vus = [r for r in mv
           if r["in_model_dataset"] != "1" and r["observed"] == "1"]
    check("VUS pool", len(vus), 717)
    check("VUS never catalogued", sum(1 for r in vus if r["in_catalogue"] == "0"), 619)
    check("VUS site-resolved",
          sum(1 for r in vus if r["struct_site_resolved"] == "1"), 317)
    check("VUS in UL30 and UL52",
          sum(1 for r in vus if r["gene"] in ("UL30", "UL52")), 517)
    check("master rows with structural features",
          sum(1 for r in mv if r["struct_rsa"] not in ("", "NA")), 966)
    check("master rows with ESM-2 scores",
          sum(1 for r in mv if r["esm2_masked_marginal"] not in ("", "NA")), 1262)
    check("VUS in UL23", sum(1 for r in vus if r["gene"] == "UL23"), 81)
    check("labelled set", sum(1 for r in mv if r["in_model_dataset"] == "1"), 235)

    # §6.4's "the frequency signal is one bit" is computed in the prose and
    # nowhere else, so it is recomputed here rather than read from a summary.
    from scipy.stats import spearmanr
    import math
    phen = [r for r in mv if r["max_fold_change"] not in ("", "NA")]
    check("variants with a numeric phenotype", len(phen), 241)
    logfc = [math.log10(float(r["max_fold_change"])) for r in phen]
    freq = [float(r["frequency"]) if r["frequency"] not in ("", "NA") else 0.0
            for r in phen]
    obs = [1.0 if r["observed"] == "1" else 0.0 for r in phen]
    check("frequency vs fold change", round(spearmanr(freq, logfc).statistic, 3),
          -0.248, 0.001)
    check("observed flag vs fold change",
          round(spearmanr(obs, logfc).statistic, 3), -0.243, 0.001)

    # release files must exist and match their recorded checksums
    import hashlib
    for rel, meta in release["files"].items():
        p = ROOT / rel
        if not p.exists():
            failures.append(f"release file missing: {rel}")
            continue
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        check(f"sha256 {rel}", h, meta["sha256"])

    # The checksums above are computed from the working tree. A reviewer sees
    # the committed blob, and git's line-ending conversion is per-platform, so
    # the two can differ and a fresh clone would fail this very check. The
    # repository sets `* -text` to prevent that; this confirms it holds.
    import subprocess
    try:
        for rel, meta in release["files"].items():
            blob = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT,
                                  capture_output=True, timeout=60)
            if blob.returncode != 0:
                continue
            got = hashlib.sha256(blob.stdout).hexdigest()
            if got != meta["sha256"]:
                failures.append(
                    f"committed blob for {rel} hashes to {got[:12]}... but "
                    f"{meta['sha256'][:12]}... is recorded; a fresh clone would "
                    f"fail verification (check .gitattributes)")
            else:
                checks += 1
    except (OSError, subprocess.SubprocessError):
        print("  (git unavailable; skipped committed-blob checksum check)")

    corr = PROC / "feature_correlation_matrix.tsv"
    if not corr.exists():
        failures.append("feature correlation matrix not published")
    for doc in ("reports/phase4/MANUSCRIPT.md",
                "reports/phase4/PRE_DECLARATION.md",
                "data/release/DATA_DICTIONARY.md", "README.md"):
        if not (ROOT / doc).exists():
            failures.append(f"missing document: {doc}")

    print(f"{checks} checks run")
    if failures:
        print(f"\n{len(failures)} MISMATCH(ES):", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    print("all reported figures match their source summaries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
