# Phase 2 report — primary analysis

**Date:** 2026-11-02
**Reproduce with:** `p2_01_reclassification.py` → `p2_02_selection.py` → `p2_03_hp_baseline.py`
**Dataset:** 465 isolates, 1,484 QC-passed gene sequences, frozen accession lists from Phase 1

---

## Verdict

**Gate 3 passed. All three primary aims produced results, two positive and one
negative.** The reclassification table is built and shows the effect the project
was designed to detect. The selection contrast supports H2 in HSV-1 with a
clean margin. H3 is not supported, and the helicase-primase baseline turned out
to be more interesting than the hypothesis attached to it.

---

## Gate 3

Threshold: at least one labelled variant above 1% in naive isolates.
**Observed: 62.** The paper does not become a confirmation study.

---

## Aim 1 — the reclassification table

644 catalogued substitutions joined to their frequency in circulating virus.
120 are observed at all; 524 are absent from 465 isolates.

| Frequency band | n |
| --- | --- |
| Majority allele (>50%) | 11 |
| Common polymorphism (10–50%) | 15 |
| Uncommon polymorphism (1–10%) | 36 |
| Rare (<1%) | 11 |
| Singleton (one isolate) | 47 |
| Not observed | 524 |

Every variant carries a leave-one-study-out range, recomputed by dropping each
contributing study in turn. This matters more than expected: it marks four of
the eight reclassification candidates as resting on a single submission.

### Ten reference artefacts

Positions where the reference strain carries the minority allele, so the
catalogued "mutation" is what nearly all circulating virus has.

| Variant | Frequency | LOSO range | Studies | Fold change |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 A265T | 0.980 (246/251) | 0.977–0.980 | 19 | — |
| HSV-2 UL30 L60P | 0.929 | 0.922–0.938 | 6 | — |
| HSV-2 UL30 P15S | 0.923 | 0.922–0.932 | 6 | — |
| HSV-1 UL23 N23S | 0.916 | 0.906–0.926 | 17 | 1.0 |
| HSV-1 UL23 K36E | 0.916 | 0.906–0.921 | 17 | 1.0 |
| HSV-1 UL30 S33G | 0.862 | 0.844–0.866 | 19 | — |
| HSV-2 UL30 A9T | 0.858 | 0.857–0.870 | 6 | — |
| HSV-1 UL30 V905M | 0.810 | 0.787–0.839 | 17 | — |
| HSV-1 UL30 T1208A | 0.798 | 0.777–0.849 | 18 | — |
| HSV-1 UL23 R89Q | 0.590 | 0.538–0.606 | 14 | — |

All ten are robust to study dropout, and all are drawn from 6 to 19 independent
studies. N23S and K36E are the cleanest cases in the whole table: measured
fold change exactly 1.0 — experimentally shown not to confer resistance — and
present in 92% of circulating HSV-1.

### One majority allele with measured resistance

**HSV-2 UL23 G39E, frequency 0.792, fold change 2.4×, Tier A, four studies,
robust.** Engineered, phenotyped, mildly resistant, and carried by four in five
circulating HSV-2 isolates. This is the single most interesting entry in the
table, because it cannot be dismissed as a reference artefact: the phenotype is
real and was measured in a clean background. It argues that a 2-fold threshold
is simply too low to be clinically meaningful for aciclovir in HSV-2.

### Eight candidates for reclassification as polymorphism

Catalogued as resistance-associated, measured fold change at or below 2×,
present above 1% in circulating virus.

| Variant | Tier | Frequency | Fold change | Study-robust |
| --- | --- | --- | --- | --- |
| HSV-1 UL30 P920S | A | 0.257 | 1.05 | yes |
| HSV-1 UL30 D672N | C | 0.229 | 1.3 | yes |
| HSV-1 UL23 V348I | C | 0.084 | 0.8 | yes |
| HSV-1 UL23 E43A | B | 0.024 | 1.3 | no |
| HSV-1 UL23 I326V | B | 0.016 | 1.3 | yes |
| HSV-2 UL30 T1200A | C | 0.016 | 1.0 | no |
| HSV-2 UL23 G150D | A | 0.013 | 2.0 | no |
| HSV-1 UL30 M1226I | A | 0.012 | 1.5 | no |

The top three are solid: ten, ten and six independent studies, robust to
dropout. The bottom four rest on one or two studies and should be reported as
provisional.

### Three circulating resistance alleles — the clinically important finding

These run the opposite way: catalogued, *measured* resistance above 2×, and
present in circulating virus above 1%.

| Variant | Tier | Frequency | LOSO range | Fold change | Robust |
| --- | --- | --- | --- | --- | --- |
| HSV-1 UL23 R222C | B | 0.076 (19/251) | 0.018–0.088 | **98×** | yes |
| HSV-1 UL23 R41H | A | 0.048 (12/251) | 0.040–0.054 | 2.3× | yes |
| HSV-1 UL23 R176Q | C | 0.016 (4/251) | 0.000–0.019 | 200× | no |

**R222C deserves direct follow-up.** A 98-fold aciclovir resistance
substitution in 7.6% of sequenced HSV-1 isolates is either an important
surveillance finding or an artefact of which cohorts get sequenced. It comes
from only two studies and the LOSO floor drops to 1.8%, so one of those studies
dominates. The phenotype is Tier B — lab-selected, so passenger mutations are
possible. Resolving this means reading the two source submissions directly.

R176Q at 200-fold is flagged not robust: it rests entirely on one study and the
LOSO floor is zero.

---

## Aim 2 — selection contrast, H2 supported in HSV-1

pN/pS computed codon-wise by Nei-Gojobori with Jukes-Cantor correction,
bootstrapped over isolates (1,000 resamples, seed 20261102).

| Gene | n | pN/pS | 95% CI | Reliable |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 thymidine kinase | 251 | **0.382** | 0.352–0.422 | yes |
| HSV-1 UL5 helicase | 256 | 0.177 | 0.158–0.201 | yes |
| HSV-1 UL30 polymerase | 253 | **0.175** | 0.161–0.190 | yes |
| HSV-1 UL52 primase | 253 | 0.134 | 0.124–0.147 | yes |
| HSV-2 UL30 | 183 | 0.507 | 0.420–0.632 | yes |
| HSV-2 UL52 | 156 | 0.549 | 0.466–0.647 | no |
| HSV-2 UL5 | 173 | 0.156 | 0.123–0.193 | no |
| HSV-2 UL23 | 159 | 2.488 | 1.768–3.918 | no |

**The Aim 2 contrast, HSV-1: UL23 0.382 against UL30 0.175, difference +0.208
with bootstrap interval +0.171 to +0.249, and UL23 higher in 100% of 1,000
resamples.** Confidence intervals do not overlap.

H2 predicted that the dispensable thymidine kinase would be under weaker
purifying selection than the essential polymerase. It is, by more than a
factor of two. This connects directly to the catalogue: UL23 holds the most
recorded resistance mutations *and* tolerates the most natural variation, which
is why distinguishing the two matters.

### A secondary finding worth reporting

Within the helicase-primase complex, **UL52 primase (0.134) is more constrained
than UL30 polymerase (0.175), while UL5 helicase (0.177) is indistinguishable
from the polymerase.** The two subunits of one complex sit at opposite ends of
the constraint range measured here. That was not predicted and is worth a
figure.

---

## Aim 3 — pre-therapeutic helicase-primase baseline

**No HSV-1 isolate carries a catalogued helicase-primase resistance
substitution.**

| | Isolates | Carrying a catalogued substitution |
| --- | --- | --- |
| HSV-1 UL5 | 256 | **0 (0.00%)** |
| HSV-1 UL52 | 253 | **0 (0.00%)** |
| HSV-2 UL5 | 173 | 2 (1.16%) |
| HSV-2 UL52 | 156 | 0 (0.00%) |

This is the baseline the aim was designed to capture, and it is a clean number:
*prior to widespread clinical use of pritelivir and amenamevir, zero of 509
HSV-1 helicase-primase gene sequences from circulating isolates carried any
catalogued resistance substitution.* Any future surveillance study has a
comparator that cannot be reconstructed once selection pressure is applied.

The contrast with UL23 and UL30 is the point: in the nucleoside-analogue
targets, 26 catalogued substitutions sit above 10% frequency. In the
helicase-primase targets, none sit above zero.

### H3 is not supported

H3 predicted that catalogued resistance positions would be more conserved than
the gene as a whole. Permutation test over position sets, 10,000 permutations:

| | Variation at catalogued positions | Gene-wide | Ratio | p |
| --- | --- | --- | --- | --- |
| HSV-1 UL5 | 0.00000 | 0.00575 | 0.00 | 0.54 |
| HSV-1 UL52 | 0.00066 | 0.00406 | 0.16 | 0.32 |
| HSV-2 UL5 | 0.00289 | 0.00085 | 3.39 | 0.039 |

The HSV-1 ratios point the predicted way but are nowhere near significance,
because only 7 and 6 positions are catalogued and the genes carry little
variation to begin with. HSV-2 UL5 is nominally significant in the *opposite*
direction on five positions, driven by the two carriers above, and should not
be read as evidence of anything.

**Reported as a negative result.** The test is underpowered by roughly an order
of magnitude, and no amount of additional public sequence will fix that while
the catalogue holds thirteen helicase-primase positions.

---

## Two methodological problems found and fixed

Both produced plausible-looking numbers that were wrong, and both are worth a
methods paragraph because anyone repeating this work will hit them.

### 1. pN/pS must be measured against the population consensus, not the reference

The first run gave HSV-2 UL30 pN/pS = 1.34 — above neutrality for an essential
polymerase, which is not credible.

Cause: in a pairwise-to-reference design, every position where the reference
strain carries a minority allele registers as a difference in nearly every
isolate. HSV-2 UL30 has three such positions, each nonsynonymous and each
present in about 93% of isolates. Between them they contributed roughly 510
nonsynonymous differences across 183 isolates.

These are not polymorphism — they are a property of the chosen reference. The
analysis now builds a consensus codon at each reference position and measures
against that, keeping the reference only for coordinates and reading frame.
HSV-2 UL30 fell from 1.34 to 0.507, and HSV-1 results tightened.

Reference-minority positions per gene: HSV-1 UL23 4, UL30 3, UL5 4, UL52 2;
HSV-2 UL23 1, UL30 3, UL5 0, UL52 2.

### 2. HSV-2 has enough genomes but not enough diversity

HSV-2 UL23 still returns pN/pS = 2.49 after the consensus fix. Investigation
rather than reporting: **the gene carries only 31 synonymous differences in
total across 159 isolates**, and 80% of its nonsynonymous differences come from
five positions — 78, 140, 39, 27 and 29, all common amino-acid polymorphisms.

The ratio is real arithmetic on a denominator too small to interpret. A
reliability rule is now applied before any gene enters the contrast: at least
50 pooled synonymous differences, and the top five sites contributing under 60%
of nonsynonymous differences. HSV-1 passes on all four genes; HSV-2 passes only
on UL30, so the HSV-2 contrast is not computed.

**This qualifies the Phase 0 conclusion that HSV-2 is a full comparative arm.**
It is, for allele frequencies, where 546 genomes give real confidence
intervals. It is not, for selection analysis, where what matters is synonymous
diversity rather than sequence count. Genome count and diversity are different
quantities and this dataset separates them clearly.

---

## Files produced

| Path | Contents |
| --- | --- |
| `data/processed/reclassification_table.tsv` | 644 labels with frequency, LOSO range, band, flag |
| `data/processed/selection_by_gene.tsv` | pN/pS per gene with CI and reliability |
| `data/processed/selection_by_site.tsv` | Per-position synonymous and nonsynonymous counts |
| `data/processed/hp_baseline.tsv` | Helicase-primase catalogued positions and frequencies |
| `reports/phase2/*.json` | Machine-readable summaries |

---

## Entering Phase 3

The modelling arm has what it needs: 103 Tier A HSV-1 positives, 33
experimentally validated negatives, a frequency-only baseline that is now
quantified, and a conservation measure from the selection analysis.

One item carries forward. **HSV-1 UL23 R222C needs its two source submissions
read before publication** — a 98-fold resistance substitution at 7.6% frequency
is either a notable surveillance result or an ascertainment artefact, and the
difference is visible only in the source records.
