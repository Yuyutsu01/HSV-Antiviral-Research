# Phase 2 report — primary analysis

**Date:** 2026-11-02
**Reproduce with:** `p2_01_reclassification.py` → `p2_02_selection.py` → `p2_03_hp_baseline.py`
**Dataset:** 599 isolates, 1,974 QC-passed gene sequences, frozen accession lists from Phase 1

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
128 are observed at all; 516 are absent from 599 isolates.

| Frequency band | n |
| --- | --- |
| Majority allele (>50%) | 11 |
| Common polymorphism (10–50%) | 15 |
| Uncommon polymorphism (1–10%) | 36 |
| Rare (<1%) | 21 |
| Singleton (one isolate) | 45 |
| Not observed | 516 |

Every variant carries a leave-one-study-out range, recomputed by dropping each
contributing study in turn. This matters more than expected: it marks three of
the seven reclassification candidates as resting on one or two submissions.

### Eleven reference artefacts

Positions where the reference strain carries the minority allele, so the
catalogued "mutation" is what nearly all circulating virus has.

| Variant | Frequency | LOSO range | Studies | Max study share | Fold change |
| --- | --- | --- | --- | --- | --- |
| HSV-1 UL23 A265T | 0.979 (234/239) | 0.975–0.996 | 34 | 0.15 | — |
| HSV-2 UL30 L60P | 0.921 (233/253) | 0.913–0.925 | 12 | 0.42 | — |
| HSV-2 UL30 P15S | 0.913 (231/253) | 0.906–0.918 | 12 | 0.42 | — |
| HSV-1 UL23 N23S | 0.912 (218/239) | 0.902–0.928 | 31 | 0.14 | 1.0 |
| HSV-1 UL23 K36E | 0.912 (218/239) | 0.906–0.928 | 31 | 0.14 | 1.0 |
| HSV-1 UL30 S33G | 0.850 (255/300) | 0.842–0.862 | 38 | 0.13 | — |
| HSV-2 UL30 A9T | 0.838 (212/253) | 0.832–0.853 | 12 | 0.42 | — |
| HSV-1 UL30 T1208A | 0.827 (248/300) | 0.812–0.849 | 36 | 0.14 | — |
| HSV-1 UL30 V905M | 0.807 (242/300) | 0.796–0.830 | 33 | 0.14 | — |
| HSV-2 UL23 G39E | 0.794 (127/160) | 0.786–0.816 | 13 | 0.47 | 2.4 |
| HSV-1 UL23 R89Q | 0.561 (134/239) | 0.531–0.575 | 24 | 0.13 | — |

All eleven are robust to study dropout, drawn from 12 to 38 independent studies,
and **none is study-dominated** — the largest single study contributes at most
15% of carriers for every HSV-1 entry. That matters because it is precisely the
property the rare-but-resistant alleles below turn out to lack. These
reclassification claims rest on broad multi-study support; those do not.

N23S and K36E are the cleanest cases in the whole table: measured fold change
exactly 1.0 — experimentally shown not to confer resistance — and present in 91%
of circulating HSV-1, across 31 studies.

### Independent corroboration of the reference artefacts

A published audit of the Sauerbrei dataset, from which HerpesDRG partly
descends, **changed the resistance status of N23S and K36E on the grounds that
it had been "reported incorrectly"** in the source literature — three
publications for N23S, one for K36E. It also flags R89Q and R281Q as having
contradictions between publications that the catalogue had not indicated.

Those are the same variants this table identifies, reached by a completely
different route: that audit read the papers, this analysis counted alleles in
599 genomes. Two independent methods converging on the same specific
mis-annotations is the strongest external support the project's premise has.

R89Q is a particularly good illustration — 56.1% population frequency across 24
studies here, and independently flagged for undisclosed contradictory phenotype
reports there.

### One majority allele with measured resistance

**HSV-2 UL23 G39E, frequency 0.794, fold change 2.4×, Tier A, thirteen carrying submissions,
robust.** Engineered, phenotyped, mildly resistant, and carried by four in five
circulating HSV-2 isolates. This is the single most interesting entry in the
table, because it cannot be dismissed as a reference artefact: the phenotype is
real and was measured in a clean background. It argues that a 2-fold threshold
is simply too low to be clinically meaningful for aciclovir in HSV-2.

### Seven candidates for reclassification as polymorphism

Catalogued as resistance-associated, measured fold change at or below 2×,
present above 1% in circulating virus.

| Variant | Tier | Frequency | Fold change | Studies | Max study share | Study-robust |
| --- | --- | --- | --- | --- | --- | --- |
| HSV-1 UL30 P920S | A | 0.233 | 1.05 | 22 | 0.13 | yes |
| HSV-1 UL30 D672N | C | 0.203 | 1.3 | 22 | 0.16 | yes |
| HSV-1 UL23 V348I | C | 0.092 | 0.8 | 13 | 0.18 | yes |
| HSV-1 UL23 E43A | B | 0.025 | 1.3 | 3 | 0.67 | no |
| HSV-1 UL23 I326V | B | 0.017 | 1.3 | 4 | 0.25 | yes |
| HSV-2 UL23 G150D | A | 0.013 | 2.0 | 2 | 0.50 | yes |
| HSV-2 UL30 T1200A | C | 0.012 | 1.0 | 1 | 1.00 | no |

The top three are solid: 22, 22 and 13 independent submissions, robust to
dropout. The bottom four rest on one to four and should be reported as
provisional — **T1200A worst of all, since every carrier comes from a single
submission.** The flag that selects this list tests frequency and fold change
but applies no ascertainment filter of its own, so a `max_study_share` of 1.00
does not keep a variant off it. Read the list with that column, not without it.

**This list had eight entries on the smaller isolate set, and HSV-1 UL30 M1226I
left it on a boundary rather than on evidence.** With 300 UL30 sequences its
three carriers give exactly 0.0100, and the band requires *above* 1%, so it
moved to `rare` and `consistent_rare`. Nothing about the variant changed. A
threshold that a single additional sequence can flip is a property of the
threshold, and M1226I should be read as sitting on the line in either
direction.

### Three candidate circulating resistance alleles — and none survives

These were reported as running the opposite way: catalogued, *measured*
resistance above 2×, and present in circulating virus above 1%. On
investigation **all three fail**, for three different reasons. This section
previously presented them as a clinically important finding; that was wrong.

| Variant | Tier | Frequency | LOSO min | Max study share | Fold change | Verdict |
| --- | --- | --- | --- | --- | --- | --- |
| HSV-1 UL23 R222C | B | 0.046 (11/239) | **0.013** | **0.73** | 98× | study-dominated |
| HSV-1 UL23 R41H | A | 0.050 (12/239) | 0.038 | 0.33 | 2.3× | **not resistant** |
| HSV-1 UL23 R176Q | C | 0.017 (4/239) | **0.000** | **1.00** | 200× | one study only |

**R222C is an ascertainment artefact.** Its apparent frequency fell from 7.6% to
4.6% once same-patient serial isolates were collapsed — 15 of its original 19
carriers came from a single submission, and their strain names (`10-6-1`,
`10-6-2`, `10-7-1`, `10-11-3`) encode one patient sampled repeatedly. Of the 11
carriers that remain, **8 come from one study**, and dropping that study takes
the estimate to 1.34%, a 3.4-fold fall. It clears the 1% robustness threshold
only barely, which is why the binary flag still reads "robust" and why that flag
should never be quoted on its own.

Its *phenotype* is not in doubt: an independent curation records R222C as
resistant with confirmed clinical susceptibility across three publications. The
98-fold figure is real. What fails is the claim that it circulates at
appreciable frequency.

**R176Q is worse** — 100% of its carriers come from one study and its frequency
falls to zero on dropout. The 200-fold phenotype has no independent population
support whatsoever.

**R41H is not a resistance allele at all.** It was the only one of the three to
survive the ascertainment checks — 33% max study share, holding at 3.8% under
dropout — but an independent curation labels it **Sensitive** for aciclovir,
penciclovir and brivudine across three publications, one with an explicit
clinical-susceptibility contradiction. Its 2.3× HerpesDRG fold change sits just
over the 2× threshold this analysis declared, which is exactly the kind of
borderline call that threshold was always going to produce.

**The honest conclusion is that this dataset provides no evidence of
high-level resistance alleles circulating at appreciable frequency in sequenced
HSV-1 genomes.** That is a weaker claim than the one this report originally
made, and a more defensible one.

### Half the observed catalogue is study-dominated

**60 of the 128 observed labelled variants more than halve in frequency when a
single study is dropped.** This is reported as a first-class limitation rather
than a footnote, because it determines which frequency estimates can carry a
claim.

Two new columns make it checkable per variant: `max_study_share` (the largest
single study's share of carriers) and `study_dominated` (whether dropout more
than halves the estimate). The split is informative rather than uniformly
discouraging — the eleven reference artefacts are all non-dominated with max
shares around 0.15, while the three failed resistance candidates sit at 0.73,
1.00 and 0.33. High-frequency alleles are supported by many studies; rare ones
usually are not. Frequency estimates below about 10% in this dataset should be
treated as provisional unless the share column says otherwise.

An implementation flaw inflated apparent robustness before this was added:
records with no PubMed ID or BioProject received an empty study key, and empty
keys were excluded from the dropout loop — so those isolates could never be
dropped, and a variant carried mostly by one study plus a few unattributed
records looked robust when it was not. The key now falls back to the accession
prefix, which keeps a submission block together and makes every isolate
droppable. R222C's "supported by 2 studies" was a product of that flaw.

---

## Aim 2 — selection contrast, H2 supported in HSV-1

pN/pS computed codon-wise by Nei-Gojobori with Jukes-Cantor correction,
bootstrapped over isolates (1,000 resamples, seed 20261102).

| Gene | n | pN/pS | 95% CI | Reliable |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 thymidine kinase | 239 | **0.371** | 0.338–0.407 | yes |
| HSV-1 UL5 helicase | 251 | 0.181 | 0.159–0.205 | yes |
| HSV-1 UL30 polymerase | 300 | **0.173** | 0.161–0.185 | yes |
| HSV-1 UL52 primase | 299 | 0.140 | 0.128–0.154 | yes |
| HSV-2 UL30 | 253 | 0.491 | 0.409–0.612 | yes |
| HSV-2 UL52 | 225 | 0.542 | 0.480–0.621 | no |
| HSV-2 UL5 | 247 | 0.167 | 0.140–0.196 | no |
| HSV-2 UL23 | 160 | 2.456 | 1.744–3.674 | no |

**The Aim 2 contrast, HSV-1: UL23 0.371 against UL30 0.173, difference +0.1980
with bootstrap interval +0.1664 to +0.2361, and UL23 higher in 100% of 1,000
resamples.** Confidence intervals do not overlap.

This result has now survived both of the corrections that moved other numbers
substantially. The patient-level de-duplication that overturned the Aim 1
resistance-allele claims left it alone (+0.208 [+0.171, +0.249] on 465 isolates,
+0.2005 [+0.1652, +0.2389] on 451), and recovering 148 further isolates moved it
by 0.0025 (+0.1980 [+0.1664, +0.2361] on 599). Selection measured across a whole
gene is insensitive both to a handful of serial isolates and to a third more
sequence, which is a useful contrast with the single-variant frequency
estimates.

H2 predicted that the dispensable thymidine kinase would be under weaker
purifying selection than the essential polymerase. It is, by more than a
factor of two. This connects directly to the catalogue: UL23 holds the most
recorded resistance mutations *and* tolerates the most natural variation, which
is why distinguishing the two matters.

### A secondary finding worth reporting

Within the helicase-primase complex, **UL52 primase (0.140) is more constrained
than UL30 polymerase (0.173), while UL5 helicase (0.181) is indistinguishable
from the polymerase.** The two subunits of one complex sit at opposite ends of
the constraint range measured here. That was not predicted and is worth a
figure.

---

## Aim 3 — pre-therapeutic helicase-primase baseline

**No HSV-1 isolate carries a catalogued helicase-primase resistance
substitution.**

| | Isolates | Carrying a catalogued substitution |
| --- | --- | --- |
| HSV-1 UL5 | 251 | **0 (0.00%)** |
| HSV-1 UL52 | 299 | **0 (0.00%)** |
| HSV-2 UL5 | 247 | 2 (0.81%) |
| HSV-2 UL52 | 225 | 0 (0.00%) |

This is the baseline the aim was designed to capture, and it is a clean number:
*prior to widespread clinical use of pritelivir and amenamevir, zero of 550
HSV-1 helicase-primase gene sequences from circulating isolates carried any
catalogued resistance substitution.* Any future surveillance study has a
comparator that cannot be reconstructed once selection pressure is applied.

The contrast with UL23 and UL30 is the point: in the nucleoside-analogue
targets, 26 catalogued substitutions sit above 10% frequency. In the
helicase-primase targets, exactly one is observed at all — the HSV-2 UL5
substitution below, at 0.8% — and none reaches 1%.

### H3 is not supported

H3 predicted that catalogued resistance positions would be more conserved than
the gene as a whole. Permutation test over position sets, 10,000 permutations:

| | Variation at catalogued positions | Gene-wide | Ratio | p |
| --- | --- | --- | --- | --- |
| HSV-1 UL5 | 0.00000 | 0.00562 | 0.00 | 0.54 |
| HSV-1 UL52 | 0.00056 | 0.00406 | 0.14 | 0.37 |
| HSV-2 UL5 | 0.00202 | 0.00092 | 2.21 | 0.052 |

The HSV-1 ratios point the predicted way but are nowhere near significance,
because only 8 and 6 positions are catalogued and the genes carry little
variation to begin with. HSV-2 UL5 runs in the *opposite* direction on four
positions, driven by the two carriers above, and at p = 0.052 is not
significant either.

**On the smaller isolate set that HSV-2 UL5 test read p = 0.039 at a ratio of
3.39, and adding 74 sequences moved it across 0.05.** The direction did not
change; the significance did. A four-position test behaving that way is the
expected behaviour, not a surprise, and it is the reason the earlier report's
instruction not to interpret it was the right call even while the number looked
interesting.

**Reported as a negative result.** The test is underpowered by roughly an order
of magnitude, and no amount of additional public sequence will fix that while
the catalogue holds fourteen HSV-1 helicase-primase positions.

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
present in 84% to 92% of isolates, so between them they contributed a
nonsynonymous difference in almost every isolate in the set.

These are not polymorphism — they are a property of the chosen reference. The
analysis now builds a consensus codon at each reference position and measures
against that, keeping the reference only for coordinates and reading frame.
HSV-2 UL30 fell from 1.34 to 0.491, and HSV-1 results tightened.

Reference-minority positions per gene: HSV-1 UL23 4, UL30 3, UL5 4, UL52 2;
HSV-2 UL23 1, UL30 3, UL5 0, UL52 2.

### 2. HSV-2 has enough genomes but not enough diversity

HSV-2 UL23 still returns pN/pS = 2.46 after the consensus fix. Investigation
rather than reporting: **the gene carries only 32 synonymous differences in
total across 160 isolates**, and 79% of its nonsynonymous differences come from
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

**The R222C question is now resolved: it is an ascertainment artefact.** It was
carried forward from this report's first version as needing its source
submissions read. Collapsing same-patient serial isolates and repairing the
leave-one-study-out key answered it without the literature step — 8 of 11
remaining carriers come from one study and the estimate falls 3.4-fold on
dropout. See the revised Aim 1 section above, where none of the three candidate
circulating resistance alleles survives.
