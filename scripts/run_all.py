"""Reproducibility runner — the whole pipeline, in order, from a clean checkout.

    python scripts/run_all.py              # everything
    python scripts/run_all.py --from p3    # resume at a phase
    python scripts/run_all.py --only p4    # one phase
    python scripts/run_all.py --list       # show stages and what they need
    python scripts/run_all.py --dry-run

Two stages are not cheap and not pure: `p1_01_fetch` goes to the network, and
`p3_01_esm` wants a GPU. Both are skipped when their outputs already exist,
unless `--force` is passed, so a reader who has the repository's data can
reproduce every number downstream of them without a network connection or a
graphics card. Everything else is deterministic given SEED = 20261102.

Exit status is non-zero if any stage fails, and the runner stops at the first
failure rather than carrying on with stale inputs.
"""

import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (phase, script, needs, produces-a-sentinel, note)
STAGES = [
    ("p0", "phase0_audit.py", "network",
     "data/processed/hsv_mutations_tiered.tsv",
     "Pins HerpesDRG, tiers mutations, Gate 1"),
    ("p1", "p1_01_fetch.py", "network",
     "data/raw/accessions/HSV1.json",
     "Frozen GenBank query and genome download"),
    ("p1", "p1_02_gate2.py", "-",
     "reports/phase1/gate2_summary.json",
     "Gate 2: catalogue positions against the reference"),
    ("p1", "p1_03_extract.py", "-",
     "data/processed/variant_calls.tsv",
     "Gene extraction, QC, three-level de-duplication"),
    ("p1", "p1_04_frequencies.py", "-",
     "data/processed/allele_frequencies.tsv",
     "Allele frequencies with Wilson intervals (Aim 1)"),
    ("p2", "p2_01_reclassification.py", "-",
     "data/processed/reclassification_table.tsv",
     "Reclassification table, LOSO, ascertainment flags"),
    ("p2", "p2_02_selection.py", "-",
     "data/processed/selection_by_gene.tsv",
     "pN/pS against population consensus (Aim 2)"),
    ("p2", "p2_03_hp_baseline.py", "-",
     "reports/phase2/hp_baseline_summary.json",
     "Pre-therapeutic helicase-primase baseline (Aim 3)"),
    ("p3", "p3_01_esm.py", "gpu",
     "data/processed/esm2_scores.tsv",
     "ESM-2 650M masked marginals, all 8 proteins"),
    ("p3", "p3_02_structure.py", "network",
     "data/processed/structure_features.tsv",
     "PDB selection, RSA, ligand distances"),
    ("p3", "p3_03_dataset.py", "-",
     "data/processed/model_dataset.tsv",
     "Assembles the labelled modelling set"),
    ("p3", "p3_04_baselines.py", "slow",
     "reports/phase3/baselines_summary.json",
     "Nested grouped CV, baselines and staged models (~30 min)"),
    ("p3", "p3_05_audit.py", "slow",
     "reports/phase3/audit_summary.json",
     "Circularity audit: permutations, matching, residualisation (~40 min)"),
    ("p3", "p3_06_gate4.py", "slow",
     "reports/phase3/gate4_summary.json",
     "Gate 4 under the two-part criterion (~30 min)"),
    ("p4", "p4_01_metrics.py", "-",
     "reports/phase4/metrics_summary.json",
     "The pre-declared metrics, from saved predictions"),
    ("p4", "p4_02_circularity_gap.py", "slow",
     "reports/phase4/circularity_gap.json",
     "Matched-negative gap against a size-matched control (~10 min)"),
    ("p4", "p4_03_release.py", "-",
     "data/release/master_variant_table.tsv",
     "Frozen accession list, master variant table, deposit"),
    ("p4", "p4_04_cross_virus.py", "-",
     "reports/phase4/cross_virus.json",
     "Cross-virus transfer validation"),
    ("p4", "p4_05_verify.py", "-",
     "reports/phase4/MANUSCRIPT.md",
     "Gate: every figure in the write-up against its source summary"),
]

SKIPPABLE = {"network", "gpu"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", help="phase to resume at, e.g. p3")
    ap.add_argument("--only", help="run a single phase")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="re-run network and GPU stages even if outputs exist")
    args = ap.parse_args()

    if args.list:
        print(f"{'phase':6s} {'script':28s} {'needs':8s} note")
        for ph, s, needs, _, note in STAGES:
            print(f"{ph:6s} {s:28s} {needs:8s} {note}")
        return 0

    phases = [p for p, *_ in STAGES]
    todo = STAGES
    if args.only:
        todo = [s for s in STAGES if s[0] == args.only]
    elif args.start:
        if args.start not in phases:
            print(f"unknown phase {args.start!r}", file=sys.stderr)
            return 2
        first = min(i for i, s in enumerate(STAGES) if s[0] == args.start)
        todo = STAGES[first:]
    if not todo:
        print("nothing to run", file=sys.stderr)
        return 2

    print(f"{len(todo)} stage(s)\n")
    t_all = time.time()
    for ph, script, needs, sentinel, note in todo:
        path = os.path.join(ROOT, sentinel)
        if needs in SKIPPABLE and os.path.exists(path) and not args.force:
            print(f"[{ph}] {script:28s} SKIP  ({needs} stage, output present)",
                  flush=True)
            continue
        label = f"[{ph}] {script:28s}"
        if args.dry_run:
            print(f"{label} would run  ({note})", flush=True)
            continue
        print(f"{label} running ...", flush=True)
        t0 = time.time()
        r = subprocess.run([sys.executable, "-u",
                            os.path.join("scripts", script)], cwd=ROOT)
        dt = time.time() - t0
        if r.returncode != 0:
            print(f"{label} FAILED after {dt/60:.1f} min "
                  f"(exit {r.returncode}) -- stopping", file=sys.stderr)
            return r.returncode
        print(f"{label} ok  {dt/60:.1f} min", flush=True)

    if not args.dry_run:
        print(f"\nall stages ok in {(time.time() - t_all)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
