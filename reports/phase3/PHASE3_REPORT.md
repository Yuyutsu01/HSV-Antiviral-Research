# Phase 3 report — Aim 4, modelling and the circularity audit

**Date:** 2026-11-02
**Reproduce with:** `p3_01_esm.py` → `p3_02_structure.py` → `p3_03_dataset.py` →
`p3_04_baselines.py` → `p3_05_audit.py` → `p3_06_gate4.py`
**Seed:** 20261102 throughout

---

## Verdict

**Gate 4 passed, but not in the way the protocol expected, and the headline
number is not the interesting part.**

Sequence- and structure-based features do predict measured fold change better
than population frequency does — under every cross-validation scheme and both
label sets, with permutation *p* ≤ 0.02 and most of the usable feature sets at
0.0099, the floor for 100 replicates.

But the protocol framed Aim 4 as "do models beat a frequency baseline?", and
that framing turned out to be built on an assumption that does not hold:
**only 14 of 119 phenotyped Tier A variants are observed in circulating virus
at all.** Frequency is the right instrument for Aim 1 and a nearly empty one
here. Beating it is a low bar, and saying a model "beat the frequency baseline"
without saying that would be misleading.

The findings that survive scrutiny are narrower and more useful:

1. Distance to the bound ligand and relative solvent accessibility are the
   strongest predictors available, and they survive conservation matching.
2. ESM-2 looks predictive on its own, but **none of its signal survives
   residualisation against conservation and structure** — the circularity
   objection holds for ESM-2 specifically.
3. Population frequency does **not** predict fold change on the phenotyped
   subset, and its signal is one bit — presence, not magnitude.
4. On the primary set, adding frequency, conservation and property features to
   structure **does not reliably help, and which way it goes depends on the
   estimator**: worse under ridge, better under random forest, and worse again
   on average precision. At 119 labels that comparison is at its resolution
   limit. The staged models do pull ahead on the 235-row secondary set, so the
   feature groups are not useless — there were not enough labels to fit them.

---

## The label budget

| Filter | Lost | Remaining |
| --- | --- | --- |
| HSV substitutions in the catalogue | — | 644 |
| Numeric fold change present | 403 | 241 |
| Tier A or B evidence (engineered or lab-selected) | — | — |
| No reported co-mutation | 5 | — |
| Reference protein available | 1 | 235 |
| **Tier A only (primary)** | | **119** |

Primary set: 119 single engineered substitutions, 87 resistant / 32 negative at
a 2× threshold, fold change 0.9× to 5000×. Secondary set adds Tier B and Tier C
for 235 rows.

Target is log₁₀ fold change. Regression, not classification: fold change is a
measurement, and collapsing a 5000× effect and a 2.1× effect into one class
discards most of what the experiments report. The binary view is retained as a
secondary column only.

**Gene coverage is uneven and limits what can be claimed:** HSV-1 UL30 45,
HSV-1 UL23 31, HSV-2 UL30 18, HSV-2 UL23 12, HSV-1 UL5 7, HSV-1 UL52 6. Two
genes carry fewer than ten labels each.

---

## Features

| Group | Features | Source |
| --- | --- | --- |
| `freq_*` | frequency, count, observed flag, CI upper | Aim 1 allele table |
| `cons_*` | per-site nonsynonymous rate, coverage | Aim 2 site table |
| `prop_*` | BLOSUM62, Δhydropathy, Δvolume, Δcharge, polarity switch, relative position, to-proline, from-glycine | amino-acid scales |
| `struct_*` | relative solvent accessibility, distance to bound ligand | PDB |
| `esm_*` | masked-marginal score, reference log-probability | ESM-2 650M |

**Gene identity is deliberately excluded as a feature.** It predicts label
prevalence directly — UL30 contributes 45 labels with a different resistance
base rate than UL52's 6 — so including it would let the model recover
prevalence instead of biology. Gene is used only as a grouping key.

### Structures, and why the active site was not hand-curated

| Gene | PDB | Chain | Coverage | Ligand contacts |
| --- | --- | --- | --- | --- |
| UL23 | 1vtk | A | 83% | 33 |
| UL30 | 8exx | A | 88% | 92 |
| UL5 | 9vlq | A | 86% | 40 |
| UL52 | 9nnp | B | 91% | **0** |

Entries were auto-selected by alignment coverage against the reference protein,
not chosen by name, and features are reported in reference coordinates after
pairwise alignment because PDB numbering routinely differs.

The active site is defined by **the ligands actually present in each
structure** rather than a list of catalytic residues. This matters for the
circularity question: a hand-written list of "known active-site residues" would
import the published claim about where resistance mutations sit, which is
precisely what Aim 4 is testing. Taking the geometry from bound ligands keeps
the feature independent of the catalogue.

It validated itself. Without any residue being named, the thymidine-pocket
residues of HSV-1 TK came out buried and in contact range: E83 2.79 Å, Y101
2.67 Å, Q125 2.85 Å, R163 2.45 Å, R176 4.96 Å, R222 3.27 Å, all with relative
accessibility ≈ 0.

**UL52 has no resolved site of its own.** Zero residues within 5 Å of any
ligand, median distance 42 Å — its nearest ligand belongs to a neighbouring
subunit of the helicase–primase complex. Rows carry `site_resolved = 0` and the
distance must not be read as an active-site feature for that gene.

HSV-2 has no structures for these genes and inherits HSV-1 geometry through
pairwise alignment; those rows are flagged `inherited = 1`.

---

## The methodological finding that changes how everything here is read

**Leave-one-gene-out manufactures negative correlation from pure noise.**

Permuting the target within cross-validation groups — which preserves each
gene's mean label and destroys only the feature-label link — gives a null
centred at **−0.16 to −0.29**, not zero, depending on feature set.

The mechanism is not subtle once seen: the model predicts approximately the
training-set mean for a held-out gene, and removing a gene shifts that mean
*away* from the gene's own mean. Whenever genes differ in average resistance,
and they do, the scheme produces anti-correlation by construction.

Consequences:

- Scored against zero, frequency-only (−0.150), conservation-only (−0.228) and
  BLOSUM-only (−0.254) look "worse than chance". They are not. They sit **on
  the noise floor** for this scheme.
- Every score under this scheme must be referenced to its own permutation null.
  A bootstrap confidence interval around zero is the wrong test and was
  initially reported as though it were the right one.
- **Position-blocks does not have this problem for the feature sets that
  matter.** Its null runs −0.014 to +0.079 across structure, ESM-2 and the
  full model, so those numbers can be read directly, and it is the sounder of
  the two schemes and carries the headline results. It is **not** a zero null
  for frequency (−0.167) or conservation (−0.115): those two have a negative
  null under both schemes, so neither scheme licenses scoring them against
  zero. They are also the two a reader most wants to compare against, which is
  why this warning is in the body rather than a footnote.

This is a leakage check too, and it passes: leakage would push the null
*above* zero. No feature set in either scheme has a null above +0.08.

---

## Results

### Out-of-fold performance, Tier A (n = 119)

| Feature set | LOGO ρ | Position-blocks ρ | Best AUC |
| --- | --- | --- | --- |
| frequency only | −0.15 / −0.15 | −0.08 / −0.10 | 0.50 |
| conservation only | −0.23 / **+0.25** | −0.07 / +0.10 | 0.69 |
| BLOSUM only | −0.25 / −0.09 | −0.16 / +0.13 | 0.60 |
| properties only | −0.03 / −0.05 | −0.07 / −0.13 | 0.48 |
| **structure only** | **+0.29 / +0.31** | **+0.48 / +0.41** | **0.81** |
| **ESM-2 only** | +0.14 / **+0.34** | **+0.31 / +0.36** | 0.66 |
| staged: freq+cons | −0.20 / **+0.34** | −0.08 / +0.13 | 0.72 |
| staged: +properties | −0.05 / −0.03 | +0.01 / +0.01 | 0.56 |
| staged: +structure | +0.15 / +0.04 | +0.30 / +0.33 | 0.75 |
| staged: +ESM (full) | +0.21 / +0.34 | +0.39 / **+0.46** | 0.77 |

Each cell is ridge / random forest.

### Gate 4, under the corrected two-part criterion

A model passes only if **both** its own Spearman interval excludes zero *and* a
paired bootstrap of its advantage over frequency excludes zero. The original
single-part test would have passed a model with no skill at all, because the
reference it was compared against scores below zero.

| Label set | Scheme | Above chance | Best model | ρ (95% CI) |
| --- | --- | --- | --- | --- |
| Tier A | LOGO | 7/22 | full model, forest | +0.336 [+0.143, +0.518] |
| Tier A | position blocks | 8/22 | structure only, ridge | +0.481 [+0.313, +0.635] |
| All tiers | LOGO | 12/22 | structure only, forest | +0.550 [+0.450, +0.634] |
| All tiers | position blocks | 14/22 | staged +structure, ridge | +0.492 [+0.383, +0.587] |

**Gate 4 passed in all four combinations.**

### Permutation tests, 100 replicates

Position-blocks, where the null is near zero and the test is cleanest:

| Feature set | Observed | Null | *p* | Usable |
| --- | --- | --- | --- | --- |
| frequency only | −0.104 | −0.167 ± 0.076 | 0.188 | no |
| conservation only | −0.069 | −0.115 ± 0.041 | 0.228 | no |
| structure only | **+0.431** | +0.020 ± 0.124 | **0.0099** | **yes** |
| structure only, forest | **+0.412** | +0.079 ± 0.140 | **0.0099** | **yes** |
| ESM-2 only | **+0.303** | +0.027 ± 0.115 | **0.0099** | **yes** |
| ESM-2 only, forest | **+0.380** | +0.075 ± 0.111 | **0.0099** | **yes** |
| full model | **+0.333** | −0.014 ± 0.112 | **0.0198** | **yes** |
| full model, forest | **+0.446** | +0.048 ± 0.138 | **0.0099** | **yes** |

**"Significant" and "usable" are separate columns on purpose.** Under
leave-one-gene-out, conservation-only reaches *p* = 0.0099 while its observed
correlation is **−0.228**. It beats a null that is more negative still
(−0.294), which means it carries a little information, and it ranks variants
backwards, which means it is useless. Reporting that as significant without
the distinction would be a real misrepresentation. Only the six structure /
ESM / full entries are usable, and the same six are usable under both schemes.

The *p* values for the two non-predictive feature sets are unstable across
re-runs at 100 replicates — frequency-only has come out anywhere between
*p* = 0.04 and *p* = 0.19 — because the statistic sits inside the null
distribution and the estimate is Monte-Carlo noise. That instability is itself
the finding: there is nothing there to be stable about. The six usable sets are
all pinned at the 1/101 floor in every run.

---

## Circularity audit

### Conservation-matched classes

Resistance mutations are catalogued at conserved, buried, active-site
positions — and so is every other functionally important residue. A model that
has only learned "this position is conserved" would look predictive. Positives
and negatives were therefore subsampled to equal counts within conservation
quartiles, removing that shortcut.

| Scheme | Unmatched ρ | Matched ρ (100 draws, n = 64) | Matched AUC |
| --- | --- | --- | --- |
| LOGO | +0.211 | +0.242 [+0.019, +0.505] | 0.643 [0.507, 0.780] |
| position blocks | +0.389 | **+0.405 [+0.081, +0.631]** | 0.736 [0.566, 0.859] |

**The signal survives.** Under both schemes the matched estimate is at least as
high as the unmatched one, with a lower bound clear of zero.

Honest statement of what this does and does not establish: the apparent skill
is not merely a conservation proxy, but matching costs roughly half the rows
(119 → 64), so this version of the test bounds the conservation contribution
loosely rather than eliminating it as an explanation. It should not be quoted
as a clean exoneration. **Phase 4 closes that gap** by pairing each matched
draw with a random draw of identical class counts, which separates matching
from sample size; the gap then comes out between −0.016 and +0.003 across all
eight model/scheme combinations — that is, at zero.

An earlier version of this test was **broken and reported a null result as a
pass** — see Bugs below.

### ESM-2 residualised against conservation and structure

ESM-2 is trained on evolutionary data, so the objection that it merely restates
the conservation features is reasonable *a priori*. **It is what happens.**

Each ESM feature was regressed on the conservation and structure features and
replaced by its residual. Conservation plus structure linearly explain only
*r*² = 0.34 of the masked-marginal score and 0.32 of the reference
log-probability — so on a variance accounting, two thirds of ESM-2 is something
else. But that two thirds carries no predictive signal:

| Scheme | ESM-2 alone | ESM-2 alone, residualised |
| --- | --- | --- |
| LOGO | ρ +0.138, AUC 0.596 | **ρ −0.113, AUC 0.411** |
| position blocks | ρ +0.305, AUC 0.660 | **ρ −0.011, AUC 0.389** |

Removing the 34% of ESM-2 that conservation and structure can explain destroys
essentially 100% of its ability to rank fold change — the residual lands either
side of zero, at −0.011 and −0.113. Whatever ESM-2 knows that
the other features do not, it is not predictive of resistance here.

The full model is more robust, because residualising ESM does not remove the
structure features themselves — the model simply falls back on them:

| Scheme | Full model | Full model, ESM residualised | Change |
| --- | --- | --- | --- |
| LOGO | ρ +0.211, AUC 0.667 | ρ +0.200, AUC 0.650 | −5% ρ |
| position blocks | ρ +0.389, AUC 0.753 | ρ +0.308, AUC 0.700 | −21% ρ |

**An earlier draft of this report read the full-model rows as evidence that
ESM-2 contributes independent signal.** That was the wrong comparison. The
full model retains structure, so its resilience says nothing about ESM. The
ESM-alone rows are the test that answers the question, and they answer it
against ESM-2.

This does not make the ESM-2 feature group useless — it reaches ρ = +0.305 out
of fold and clears its permutation null at the floor, so it is a legitimate
*proxy*. It means ESM-2 is not an independent line of evidence, and a claim
that a language model has learned something about resistance beyond
conservation and structural context is not supported by these data.

### Leakage

Passed. Leakage would show as a permutation null **above** zero. Maximum null
mean observed: +0.079 (position-blocks) and −0.159 (LOGO). No evidence of
information crossing fold boundaries.

---

## What this means for Aim 4

The protocol's question was whether sequence-based models beat a frequency-only
baseline. The answer is yes, decisively, and it is less informative than it
sounds. Three statements are better supported:

**Structural context predicts measured resistance; population frequency does
not.** Distance to the bound ligand plus solvent accessibility reaches ρ = 0.48
and AUC 0.81 out of fold. Frequency reaches ρ = −0.08 and AUC 0.50, and fails
its own permutation test (*p* = 0.19).

**This does not contradict Aim 1.** Frequency answers "is this catalogued
mutation actually a common polymorphism?" — a question about the catalogue.
It does not answer "how resistant is this variant?" The two aims operate on
nearly disjoint subsets: 128 of 644 catalogued substitutions are observed in
circulating virus, but only 14 of 119 phenotyped Tier A ones. Conflating these
would be the single easiest mistake to make in writing this up.

**The frequency signal is one bit.** Across all 241 variants with a numeric
phenotype, continuous frequency gives ρ = −0.248 (*p* = 0.0001) and the binary
observed/not-observed flag gives ρ = −0.243 (*p* = 0.0001). Identical. Knowing
how common an allele is adds nothing over knowing whether it has ever been
seen. Median fold change is 1.5× for observed variants and 5.0× for unobserved.

**Whether more features help depends on the estimator and the metric, which is
itself the finding.** Under ridge, structure alone beats the full staged model
on the primary set under both schemes (+0.29 vs +0.21 leave-one-gene-out;
+0.48 vs +0.39 position blocks). Under random forest the full model wins
(+0.34 vs +0.31 and +0.46 vs +0.41). And on average precision — the metric
pre-declared as primary, scored in Phase 4 — structure-only is best either way
(0.916 against the full model's 0.893).

A result that changes rank with the estimator is a result at the resolution
limit of 119 labels with 32 negatives, and should be reported as one. The
defensible reading is that structure carries the signal and the other three
feature groups add nothing reliable on the primary set; whether they add noise
or a little signal is not resolvable here. Hyperparameters were selected in an
inner grouped loop throughout, so this is not a tuning failure.

The secondary set is consistent with that reading. At n = 235 under position
blocks the staged models pull ahead (+0.492 for `+structure` ridge against
structure-only's +0.386 ridge, though structure-only's forest at +0.488 is
within noise of both), and `properties_only` goes from failing the gate to
passing it (+0.186). The extra feature groups are not noise in themselves —
there were not enough labels in the primary set to fit them.

---

## Limitations that must be reported

**The sample is small and the genes are unbalanced.** 119 primary labels across
six virus–gene combinations, two of which have fewer than ten. Any claim about
cross-gene generalisation rests on six groups.

**Only 14 of 119 primary variants are observed in circulating virus**, so the
frequency features are near-constant on the primary set and the frequency
baseline is weak for structural reasons, not because frequency is a bad idea.

**Tier A labels are engineered genotypes measured in different laboratories**
with different assays and cell lines. `max_fold_change` pools across drugs and
studies. Between-study assay variation is unmodelled and is a plausible
contributor to the residual.

**UL52 has no resolved active site** and its distance feature is not comparable
to the other three genes'.

**HSV-2 structural features are inherited from HSV-1** by alignment, so they
carry the assumption that the two orthologues have the same fold — reasonable
at 70–80% identity, but it is an assumption.

**The conservation-matched test is underpowered** (n ≈ 64) and bounds the
conservation contribution loosely.

**Leave-one-gene-out results must never be quoted against a zero null.** Doing
so inverts the sign of the conclusion for the weaker feature sets.

---

## Bugs found and fixed during this phase

Recorded because every one of these produced a pipeline that ran to completion
and reported success while being wrong — the same failure mode as Phases 1
and 2.

1. **The conservation-matched test silently tested nothing.** It drew negatives
   to match the positives' conservation distribution, but negatives are the
   minority class here (32 vs 87), so every per-bin request exceeded the pool,
   `min(k, len(pool))` returned all of them, and all 200 draws were the
   identical full dataset. The giveaway was a confidence interval of
   `[0.2663, 0.2663]`. Now both classes are subsampled to the smaller per-bin
   count. **This had reported a circularity pass that had not been earned.**
2. **Gate 4 compared against a baseline that scores below zero.** The original
   criterion was "does any model beat frequency-only?". With frequency at
   −0.153, a model at +0.05 would have passed while having no predictive skill.
   Gate 4 now requires skill above chance *and* advantage over frequency.
3. **Scores were referenced to zero when the null is −0.21.** This inverted the
   reading of three baselines, which were reported as "at or below chance" when
   they were on the noise floor. Fixed by permutation testing against the
   scheme's own null.
4. **Hyperparameters were being tuned inside the permutation null.** Tuning
   against a shuffled target is meaningless and makes the test invalid.
   Observed and permuted statistics now come from an identical fixed-parameter
   procedure.
5. **`significant_vs_null` conflated "beats the null" with "predicts".** A
   feature set with ρ = −0.21 was labelled significant. A separate `usable`
   flag now requires a positive observed correlation.
6. **ESM-2 was spilling GPU memory and running ~40× slow.** Batch 32 at
   1000 tokens exceeded 6 GB and fell back to host memory: 13 minutes per
   window with no error raised. fp16 and batch 8 cut a whole protein from
   >13 minutes to 19 seconds. Per-protein checkpointing was added, because the
   original wrote output only at the very end and an OOM at protein seven would
   have discarded everything.
7. **UL52 distances were measuring another subunit's ligand.** Added the
   `site_resolved` flag rather than letting the column look comparable.
8. **The ESM-2 residualisation result was read off the wrong rows.** The audit
   emits both a full-model and an ESM-alone residualisation; an earlier draft
   quoted the full-model pair (−5% to −21%) as evidence that ESM-2 carries
   independent signal. The full model still has the structure features after
   ESM is residualised, so it cannot test that. The ESM-alone pair
   (+0.31 → −0.01) shows the opposite conclusion. No code was wrong here —
   the output was correct and the interpretation was not, which is the harder
   failure to catch.

---

## Files produced

| Path | Contents |
| --- | --- |
| `data/processed/esm2_scores.tsv` | 135,166 masked-marginal scores, all 8 proteins |
| `data/processed/structure_features.tsv` | 6,235 rows: RSA, ligand distance, site flag |
| `data/processed/model_dataset.tsv` | 235 labelled rows with all feature groups |
| `data/processed/model_predictions_primary_tier_A.tsv` | out-of-fold predictions, 44 model/scheme combinations |
| `data/processed/model_predictions_secondary_all_tiers.tsv` | same for the 235-row set |
| `reports/phase3/structure_summary.json` | PDB selection, coverage, contact counts |
| `reports/phase3/dataset_summary.json` | label budget and exclusions |
| `reports/phase3/baselines_summary.json` | all baseline and staged results |
| `reports/phase3/gate4_summary.json` | corrected Gate 4 verdicts |
| `reports/phase3/audit_summary.json` | permutation nulls, matching, residualisation |
| `reports/audit_run_200perm_partial.log` | higher-powered 200-replicate LOGO nulls |

---

## Outstanding before publication

1. ~~HSV-1 UL23 R222C still needs its two source submissions read.~~
   **Resolved — it is an ascertainment artefact.** Patient-level collapse took
   it from 19 carriers to 11 and 7.6% to 4.6%; 8 of the 11 come from one study
   (`max_study_share` 0.73), and leave-one-study-out drops the frequency 3.4×
   to 1.3%. R176Q is worse (study share 1.00, LOSO minimum 0.0) and R41H is
   labelled Sensitive by the published audit. **None of the three candidate
   circulating resistance alleles survives.** The structural result was
   consistent with resistance — R222 sits 3.27 Å from the bound ligand with
   RSA ≈ 0 — which is exactly why the frequency side had to be settled
   separately. See `PHASE2_REPORT.md`.
2. **Assay heterogeneity in the Tier A labels is unmodelled.** A random effect
   per study would be the honest treatment. `statsmodels` is not installed in
   this environment; PMIDs for most Tier A rows are now recoverable from the
   Dähne et al. supplement (item 4), so this is tractable rather than blocked.
3. **The matched-subset test deserves more power.** Relaxing matching from
   quartiles to tertiles would retain more rows at the cost of coarser
   matching; worth reporting both.
4. **The Dähne et al. (2025) supplementary tables are not yet integrated.**
   Two independent uses, in order of value: (a) cross-catalogue agreement as a
   result in its own right — their audit changed the status of N23S and K36E as
   "reported incorrectly" and flagged R89Q and R281Q for unindicated clinical
   contradictions, and this project's population analysis independently puts
   those at 91.2%, 91.2%, 56.1% and 6.3%; (b) ~500 substitutions they label
   Sensitive, against the 32 negatives available here, which would give the
   binary form of Aim 4 real power for the first time. Their coordinate system
   matches this project's exactly (NC_001806 / NC_001798). Provenance and
   licence are recorded in `data/raw/PROVENANCE.md`.
5. **Gate 2's 100% pass rate is narrower than it looks.** It verifies that the
   catalogue is internally consistent with the reference sequence, not that the
   catalogue is faithful to its sources. `L139V` and `A164T` passed this gate
   while being flagged by the published audit.
