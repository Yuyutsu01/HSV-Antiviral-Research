# Pre-declaration record and protocol compliance audit

**Written:** 2026-11-02
**Covers:** the analysis plan declared before data acquisition, what was
actually run, and every point on which the two differ.

The protocol requires that "the primary metric and the primary model are
declared before any evaluation is run, and recorded in the repository with a
timestamp. No selection of either on test performance." This file is that
record. It is written after the analysis, so it is a compliance audit as well
as a declaration, and it is written to be read by someone looking for the places
where the plan was not followed.

---

## What the timestamps do and do not establish

| Event | Time (UTC) | Evidence |
| --- | --- | --- |
| Analysis plan authored, including the metric hierarchy | 2026-11-01T18:31:54Z | session transcript |
| HerpesDRG database retrieved and pinned | 2026-11-01T19:30:05Z | `data/raw/PROVENANCE.md` |
| GenBank query executed, HSV-1 | 2026-11-01T20:02:36Z | `data/raw/accessions/HSV1.json` |
| GenBank query executed, HSV-2 | 2026-11-01T20:07:23Z | `data/raw/accessions/HSV2.json` |
| First model evaluation | 2026-11-02 | `reports/phase3/baselines_summary.json` |

**The ordering is real: the plan predates data acquisition by about an hour and
model evaluation by a day.** The metric hierarchy, the 2× fold-change
threshold, the evidence tiers, the gate criteria and the commitment to report a
negative Aim 4 result were all fixed before any outcome was visible.

**This is not a public pre-registration, and must not be described as one.** It
was not deposited with a third party before the analysis. The ordering is
evidenced only by file timestamps and the session transcript, both of which are
under the authors' control. A reader who does not trust the authors has no
independent confirmation. Stating this plainly is the point: the protections a
pre-registration provides against outcome-dependent analysis choices are
*partially* in place here, not fully.

---

## Declared primary metric and model

**Primary metric: average precision**, always reported beside the positive-class
prevalence, which is the no-skill floor. The protocol's reasoning was that
"positives are rare relative to the candidate pool", so AUROC would be
dominated by easy true negatives.

**That premise turned out to be false, in the opposite direction.** On the
phenotyped set positives are the *majority*: prevalence is 0.731 on the primary
Tier A set and 0.668 on the secondary set. The imbalance runs the other way
than anticipated, which means average precision is *inflated* by the high
positive rate rather than deflated by a rare one. An AP of 0.916 against a
0.731 floor is a lift of +0.185, not an 0.916-sized achievement.

The metric was kept anyway, because changing the primary metric after seeing
which one flatters the result is exactly what pre-declaration exists to
prevent. The correct response is to report AP beside its floor, as declared,
and to say that the declared justification does not hold. Both are done.

**Primary model: regularised logistic regression (elastic net), with
gradient-boosted trees as Stage 2, justified only by capacity.**

**Implemented: ridge regression and random forests.** This is a deviation, in
two parts:

- *Logistic to ridge* follows from the protocol's own choice of regression as
  the primary task ("the primary task is regression on log fold-change; binary
  classification at a declared 2× threshold is secondary"). A logistic model
  cannot fit a continuous target. Ridge is the L2 analogue of the declared
  elastic net with the L1 component dropped. Reasonable, but it is a change and
  is recorded as one.
- *Random forests in place of gradient-boosted trees* is a substantive
  deviation. The protocol admitted a tree ensemble only as Stage 2 and only if
  justified by capacity. Random forests were run alongside ridge from the
  start, not after ridge was found wanting, and they win in several
  scheme/metric combinations. Both model families are reported everywhere so
  that a reader can apply the declared rule themselves: **under the declared
  rule, ridge is the primary model, and on the primary Tier A set under
  position blocks it is also the better of the two for structure-only
  features** (AP 0.905 ridge against 0.916 forest, intervals overlapping
  heavily), so the deviation does not change the headline.

---

## Compliance, item by item

### Followed as declared

| Item | Where |
| --- | --- |
| Average precision as primary, reported beside prevalence | `reports/phase4/metrics_summary.json` |
| Precision at 95% recall as the clinical operating point | same |
| Matthews correlation coefficient | same |
| Accuracy and F1 **not** reported | — |
| AUROC reported for comparability, no claim resting on it | labelled `auroc_not_primary` in the JSON |
| Bootstrap intervals resampled over **mutations**, not predictions | `bootstrap_ci()` in `p4_01_metrics.py` |
| Breakdown by evidence tier, gene and drug mechanism class | `breakdowns` in the metrics JSON |
| Feature correlation matrix published in full | `data/processed/feature_correlation_matrix.tsv` |
| ESM-2 with conservation regressed out | `reports/phase3/audit_summary.json` |
| Matched-negative gap reported as a headline number | `reports/phase4/circularity_gap.json` |
| Permuted-label runs for leakage | `reports/phase3/audit_summary.json` |
| Paired comparison with multiplicity correction | Holm, in `p4_01_metrics.py` |
| Deep networks and fine-tuned language models excluded; ESM-2 zero-shot only | `p3_01_esm.py` |
| Negative Aim 4 result to be reported as the finding | reported; see below |
| Circularity audit in main text, not a supplement | `PHASE3_REPORT.md`, `MANUSCRIPT.md` |

### Deviations

1. **Primary model substitution**, as set out above.
2. **Phase 3 led with Spearman and AUROC.** The pre-declared primary metric was
   not computed until Phase 4. The Phase 3 report's framing therefore rested on
   a metric the protocol had explicitly demoted. Phase 4 recomputes every
   model on the declared metrics from the saved out-of-fold predictions — the
   same fits, scored as declared — and the headline changes materially: under
   position blocks on the primary set, 8 of 22 models cleared the Phase 3
   criterion, but only **2 of 21** beat the frequency baseline on average
   precision after Holm correction. Both are structure-only.
3. **Multiplicity was uncorrected in Phase 3.** Added in Phase 4. This is the
   single change that most reduced the apparent result, and the protocol
   anticipated it in terms: "an uncorrected pairwise win on one split is not
   evidence of superiority".
4. **The matched-negative control was not size-matched in Phase 3.** The audit
   compared a conservation-matched subset of 64 rows against the full 119, which
   confounds matching with sample size. Phase 4 pairs each matched draw with a
   random draw of identical class counts. The gap then comes out at between
   −0.015 and +0.003 across all eight model/scheme combinations: **essentially
   zero**. Phase 3's hedged conclusion, that matching "bounds the conservation
   contribution loosely rather than eliminating it", was an artefact of the
   unmatched sample size.
5. **Brier score and calibration curve not computed.** The declared primary
   task is regression, so the models emit fold-change estimates rather than
   probabilities and a Brier score is undefined for them. Fitting a calibrator
   on the same out-of-fold predictions used for evaluation would leak. The
   regression calibration — slope and intercept of observed on predicted — is
   reported instead; the slope is 0.71 to 0.80 throughout, meaning the models
   systematically under-disperse, as expected at this sample size.
6. **SHAP / permutation importance not computed.** The protocol made this
   conditional. With only two structural features carrying the result, the
   correlation matrix and the staged feature-group comparison already answer
   what importance ranking would.

### Temporal hold-out: run, with one choice made after the fact

The protocol calls this "the single most persuasive validation in the plan". It
was initially blocked: publication year is recoverable from the DOI string for
only 35 of 89 distinct studies, and the recoverable ones cluster by publisher,
so a partial split would have confounded era with journal. Resolving the DOIs
through Crossref returned a year for **89 of 89**, spanning 1987-2025 across 10
publishers, so all 119 primary and 235 secondary variants carry a
first-characterisation year and the confound is gone. Years are cached in
`data/raw/doi_years.json` with their retrieval date.

**The cut-off year was not pre-declared, and this is the one place in the
project where an analysis parameter was chosen after seeing data.** The
mitigation is that it was not chosen on performance: all five cut-offs from
2005 to 2015 are reported, and the nominated primary (2012) is the one that
splits the primary set most evenly, which is an outcome-independent criterion.
The result does not depend on the choice — structure-only is the best feature
group at every cut-off on the primary set and at four of five on the
secondary — so a reader who distrusts the nomination can read the sensitivity
table instead and reach the same conclusion.

Only DOIs were sent to Crossref. No `mailto` parameter and no address in the
User-Agent: that would place a personal identifier in third-party request logs
for no analytical gain.

### Pre-declared and NOT done

**Positive-unlabelled formulation (Stage 3) — not attempted, deliberately.**

It was declared as the extension where the clinical value lies, and that
judgement was right: a clinical laboratory meeting a substitution in a patient
isolate is far more likely to be looking at an unphenotyped variant than a
characterised one. Of the variants observed in circulating virus, **610 have no
measured fold change**, 518 of them not catalogued at all, against 235 with a
phenotype. Extending scoring to that pool is the obvious next step and it is
not taken here, for three reasons that are worth recording in full because a
reader will otherwise assume it was simply overlooked.

**1. The standard PU assumption is violated by this project's own subject
matter.** Most positive-unlabelled methods assume positives are selected
completely at random from all positives. Here a variant is catalogued *because*
someone found it in a resistant isolate. Ascertainment bias is precisely what
this project set out to measure — 59 of 120 observed variants are
study-dominated — so applying a method that assumes it away, inside a paper
about it, would be incoherent. Methods that relax the assumption exist and need
more positives than the 87 available.

**2. The class prior cannot be pinned down at this sample size.** PU
performance depends almost entirely on the estimated fraction of the unlabelled
pool that is positive. With 87 positives that estimate would carry an interval
wide enough to make the resulting ranking unstable, and an unstable ranking
presented as a triage list is worse than no list.

**3. The one feature group that works does not cover the pool.** Section 5.3 of
the manuscript shows the structural signal is confined to genes with a
ligand-resolved active site: AP 0.988 for UL23, 0.567 for UL5, 0.411 for UL52.
The VUS pool is concentrated in the wrong genes — UL30 and UL52 carry 418 of
the 610, UL23 only 81 — and just **282 of the 610 have `site_resolved = 1`**. A
score over UL52 variants of unknown significance would be close to noise while
looking authoritative. Combined with calibration slopes of 0.71–0.80 and
degeneracy at the 2× threshold, a probability attached to a patient's variant
would not mean what it says.

**What is delivered instead**, with no modelling assumptions at all: the 40
`unphenotyped_but_common` variants in the reclassification table. These are
common in circulation and have never been phenotyped, which makes them a
phenotyping-priority list derived from observation rather than from a model. If
the PU extension is attempted later, the defensible scope is a ranked list
restricted to the 282 site-resolved VUS, published as a phenotyping priority
and not as clinical interpretation.

---

## The pre-declared negative result

The protocol committed in advance to reporting a negative Aim 4 result as the
finding if the models failed to beat the frequency baseline. The models did beat
it. But the protocol's framing assumed a frequency baseline with something to
say, and **only 13 of 119 phenotyped Tier A variants are observed in
circulating virus at all**, so the baseline is near-constant on the primary set.

The commitment is therefore honoured in a different place than expected: the
negative result to report is not "models fail to beat frequency" but
**"frequency does not predict phenotype, and beating it establishes very
little"**, together with **"at the pre-declared 2× operating point the best
model makes no usable decision"** — Matthews correlation 0.15 on the primary
set, with 9 of 22 models classifying every variant identically at that
threshold. Good ranking and an unusable operating point are compatible, and the
protocol's choice to require precision at 95% recall and MCC alongside a
ranking metric is what exposed it.
