# Most catalogued HSV antiviral resistance mutations have never been checked against the virus that circulates

**Draft manuscript — Phase 4**
**Date:** 2026-11-02
**Status:** complete draft for internal review. Two pre-declared validations
are outstanding; see Limitations.

---

## Abstract

Catalogues of antiviral resistance mutations in herpes simplex virus are
assembled from case reports and phenotyping studies. They record that a
substitution was present in a resistant isolate, but not how common that
substitution is in virus from people who were never treated. Without that
denominator a natural polymorphism and a resistance determinant are
indistinguishable.

We assembled allele frequencies for every catalogued HSV mutation in thymidine
kinase (UL23), DNA polymerase (UL30), helicase (UL5) and primase (UL52) from
451 publicly deposited genomes, and used them to recontextualise the catalogue.

Of 644 catalogued substitutions, **120 occur in circulating virus and 524 do
not**. Eleven are the *majority* allele, meaning the reference sequence carries
the minority residue; a further 15 exceed 5%. Eight are recommended for
reclassification as polymorphisms. Forty are common in circulation but have
never been phenotyped, and are the clearest target for future assay work.
Three initially appeared to be circulating resistance alleles; none survives
scrutiny, two as single-submission artefacts and one on a 2.3-fold phenotype.

Thymidine kinase is under markedly weaker purifying selection than DNA
polymerase (pN/pS 0.372 against 0.172; difference +0.2005, 95% CI +0.1652 to
+0.2389, 100% of bootstrap resamples), consistent with its dispensability for
replication in vitro and with the greater tolerance of its sequence for
substitution.

No isolate among 242 HSV-1 UL5 and 239 HSV-1 UL52 sequences carries a
catalogued helicase–primase resistance substitution, establishing a
pre-therapeutic baseline of effectively zero before pritelivir and amenamevir
enter wide use.

Finally, we asked whether sequence- and structure-based features predict
measured fold change better than population frequency. They do, but the finding
that matters is narrower: distance to the bound ligand and solvent accessibility
reach an average precision of 0.916 against a no-skill floor of 0.731, while
population frequency reaches 0.723 — **below its own floor**. The structural
signal is confined to thymidine kinase, where a nucleoside pocket is resolved
(AP 0.988); it disappears for helicase and primase, where it is not (0.567 and
0.411). At the pre-declared 2-fold operating point the best model's Matthews
correlation is 0.15, and 9 of 22 models classify every variant identically. A
model that ranks well and decides uselessly is the honest summary.

---

## 1. Why the denominator is missing

A resistance catalogue is a record of ascertainment. An entry means someone
sequenced a resistant isolate, found a substitution, and reported it. Nothing in
that process establishes how often the substitution occurs anyway.

The consequence is not hypothetical. Eleven substitutions catalogued as
resistance-associated in HSV are the commonest residue at their position in
circulating virus. For those, the catalogue has recorded a property of the
reference sequence rather than of the virus.

The gap persists because the two kinds of data sit in different places.
Resistance phenotypes live in the primary literature and in curated databases;
population allele frequencies require assembling and aligning public genomes.
Nobody is being careless. The denominator is simply a different task, and this
is that task.

## 2. Methods

### 2.1 Sources, frozen

HerpesDRG, pinned at commit `22fc186` (retrieved 2026-11-01T19:30:05Z, MD5
`d78ad69a…`), active entries only (`status == "A"`). GenBank nucleotide records
for taxa 10298 (HSV-1) and 10310 (HSV-2) restricted to 140–200 kb, queried
2026-11-01T20:02:36Z and 20:07:23Z. All 1,143 returned accessions are listed in
`data/release/accession_list.tsv` **with their fate**, so that an absent
accession can be distinguished from an excluded one.

### 2.2 From 1,143 accessions to 451 isolates

| Fate | n |
| --- | --- |
| Analysed | **451** |
| No target-gene CDS extracted | 411 |
| Excluded by a downstream coverage or divergence filter | 183 |
| Failed CDS quality control | 36 |
| Laboratory-adapted strain | 30 |
| Serial isolate from an already-represented patient | 18 |
| Non-human host | 14 |

De-duplication operates at three levels: accession, strain and patient. The
patient level matters more than it sounds — see §6.1.

Per-gene coverage: HSV-1 237–242 sequences, HSV-2 156–183.

### 2.3 Evidence tiers

Tier is taken from `strain_generation`, keeping the strongest available per
mutation: **A, engineered** (recombinant BAC or marker transfer; 133
substitutions), **B, lab-selected** (89), **C, clinical isolate** (478). Tier A
is the primary set for anything involving phenotype, because only there is the
genotype–phenotype link experimentally isolated.

### 2.4 Alignment, and a constraint worth declaring

No multiple-sequence-alignment tool was available. Alignment is therefore
reference-anchored and pairwise (Biopython `PairwiseAligner`, BLOSUM62, free
terminal gaps). For coding sequences of 70–80% identity to a known reference
this is adequate, and it makes every position directly comparable to reference
coordinates without a column-mapping step. It is nonetheless a limitation, and
indels are handled less well than a codon-aware MSA would handle them.

### 2.5 Frequencies and ascertainment

Frequencies are counts over per-position coverage, with 95% Wilson score
intervals. Every frequency additionally carries:

- `n_studies_carrying` — distinct submissions contributing a carrier;
- `loso_min_freq` — the lowest frequency under leave-one-study-out;
- `max_study_share` — the largest share of carriers from one submission;
- `study_dominated` — whether dropping one study halves the frequency.

These are not supplementary. A 5% allele drawn entirely from one submission and
a 5% allele seen across five are different claims, and the table refuses to
state one without the other.

### 2.6 Selection

pN/pS by Nei–Gojobori with Jukes–Cantor correction, measured **against the
population consensus rather than the reference**. Using the reference would
count every position where the reference is the minority residue as
polymorphism, inflating the estimate by a term that has nothing to do with
selection. A gene is reported as reliable only when synonymous sites are
sufficient; HSV-2 UL23 (Sd = 31) is not, and the HSV-2 contrast is withheld
rather than reported with a wide interval.

### 2.7 Features and models for Aim 4

Five feature groups: frequency, conservation (per-site non-synonymous rate),
amino-acid properties, structure, and ESM-2.

**Gene identity is excluded as a feature.** It predicts label prevalence
directly — UL30 contributes 45 labels at a different base rate from UL52's 6 —
so including it would let a model recover prevalence instead of biology. Gene
is a grouping key only.

**The active site is defined by the ligands present in each structure, not by a
list of catalytic residues.** A hand-written list would import the published
claim about where resistance mutations sit, which is what Aim 4 is testing.
Structures were auto-selected by alignment coverage: 1vtk (UL23), 8exx (UL30),
9vlq (UL5), 9nnp (UL52). The geometry validated itself without any residue
being named: the HSV-1 thymidine pocket residues came out buried and in contact
range (E83 2.79 Å, Y101 2.67 Å, Q125 2.85 Å, R163 2.45 Å, R222 3.27 Å, all with
relative accessibility ≈ 0).

**UL52 has no resolved site of its own** — zero residues within 5 Å of any
ligand, median distance 42 Å, its nearest ligand belonging to a neighbouring
subunit. Those rows carry `site_resolved = 0`. This turns out to matter a great
deal (§5.3).

ESM-2 `esm2_t33_650M_UR50D` masked marginals, zero-shot. Deep networks and
fine-tuned language models were excluded in advance: at a few hundred labelled
positives they cannot be fitted or validated.

Two cross-validation schemes, both grouped: leave-one-gene-out, and
leave-one-block-out over 50-residue position blocks. Hyperparameters are tuned
in an inner grouped loop. The primary metric is average precision, reported
beside the positive-class prevalence; the 2-fold threshold and the metric
hierarchy were fixed before data acquisition (`PRE_DECLARATION.md`).

## 3. Aim 1 — what the catalogue looks like with a denominator

Of 644 catalogued substitutions, 120 occur in circulating virus.

| Band | n | |
| --- | --- | --- |
| Not observed | 524 | No population evidence either way |
| Singleton | 47 | One carrier |
| Rare (<1%) | 11 | |
| Uncommon (1–5%) | 36 | |
| **Common (5–50%)** | **15** | |
| **Majority (>50%)** | **11** | The reference carries the minority residue |

Interpretive flags:

| Flag | n | What it means |
| --- | --- | --- |
| `reference_artifact` | 10 | Majority allele, no measured resistance |
| `reclassify_as_polymorphism` | 8 | Catalogued as resistance-associated, common, fold change below 2× |
| `unphenotyped_but_common` | 40 | Common in circulation, never phenotyped |
| `consistent_rare` | 53 | Rare, as a resistance determinant should be |
| `consistent_with_resistance` | 5 | Rare and phenotypically resistant |
| `circulating_resistance_allele` | 3 | Common **and** resistant — all three fail on inspection |
| `majority_allele_with_measured_resistance` | 1 | |

**The eight `reclassify_as_polymorphism` entries are the paper's concrete
deliverable**, and the 40 `unphenotyped_but_common` entries are its clearest
call to action: these are common enough that a laboratory encountering them will
have to interpret them, and no phenotype exists.

### 3.1 The three apparent circulating resistance alleles

| Variant | Freq | Carriers | Studies | `max_study_share` | `loso_min` | FC | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| UL23 R222C | 0.046 | 11 | 3 | 0.73 | 0.0135 | 98× | Ascertainment artefact |
| UL23 R176Q | 0.017 | 4 | 1 | 1.00 | 0.0 | 200× | Ascertainment artefact |
| UL23 R41H | 0.051 | 12 | 5 | 0.33 | 0.0385 | 2.3× | Sound frequency, marginal phenotype |

R222C was the single most consequential row in the table: 98-fold resistance at
an apparent 7.6% frequency would have been a surveillance finding. It is not.
Patient-level collapse reduced it from 19 carriers to 11 and 7.6% to 4.6%, and 8
of the 11 remaining carriers come from one submission; dropping that submission
takes the frequency to 1.4%.

R41H is the interesting case, because it is **not** an artefact — 12 carriers
across five submissions, robust to study dropout. What disqualifies it is the
phenotype: 2.3-fold is marginal against the declared 2-fold threshold, and an
independent published audit labels it Sensitive. Frequency and phenotype each
need their own scrutiny, and `study_dominated = 0` alone would have retained it.

**59 of the 120 observed variants are study-dominated.** Any use of these
frequencies that ignores that column will overstate about half of them. The 11
reference artefacts, by contrast, are all non-dominated at study shares around
0.15 — they are real.

## 4. Aims 2 and 3

### 4.1 Thymidine kinase is under weaker constraint than polymerase

| Gene | pN/pS | 95% CI | Sd |
| --- | --- | --- | --- |
| HSV-1 UL23 | **0.372** | 0.339–0.406 | 609 |
| HSV-1 UL30 | **0.172** | 0.157–0.187 | 1251 |
| HSV-1 UL5 | 0.176 | 0.155–0.198 | 935 |
| HSV-1 UL52 | 0.134 | 0.122–0.147 | 1647 |

Difference UL23 − UL30: **+0.2005 [+0.1652, +0.2389]**, with UL23 higher in
**100%** of 1,000 bootstrap resamples.

This is the expected direction — thymidine kinase is dispensable for replication
in culture while polymerase is not — and it supplies a population-genetic
reason why TK accumulates substitutions that then appear in resistance
catalogues. It is also robust: the patient-level de-duplication fix that moved
R222C by 40% moved this contrast by 0.004, because gene-wide rates are
insensitive to serial isolates in a way single-variant frequencies are not.

HSV-2 is withheld: UL23 has only 31 synonymous sites and a pN/pS point estimate
of 2.49 with an interval spanning 1.77–3.92, which is an artefact of the
denominator rather than evidence of diversifying selection.

### 4.2 The helicase–primase baseline is effectively zero

| Gene | Isolates | Carrying a catalogued substitution |
| --- | --- | --- |
| HSV-1 UL5 | 242 | **0** |
| HSV-1 UL52 | 239 | **0** |
| HSV-2 UL5 | 173 | 2 (1.2%) |
| HSV-2 UL52 | 156 | **0** |

Pritelivir and amenamevir are new enough that this is a genuine pre-therapeutic
baseline, and it is the cleanest result in the project: there is essentially no
pre-existing helicase–primase resistance in circulating HSV. Any future rise is
measurable against this.

One caveat against over-reading it: the catalogue lists only 8 HSV-1 UL5 and 6
UL52 positions, so the test has little to find. A permutation test of whether
catalogued positions are more conserved than random ones is non-significant for
HSV-1 UL5 (p = 0.53) and UL52 (p = 0.32). For HSV-2 UL5 catalogued positions
are 3.4× *more* variable than gene-wide (p = 0.039) — the opposite of
constraint — but that rests on four positions and should not be interpreted.

## 5. Aim 4 — models, and what beating a frequency baseline is worth

### 5.1 The protocol's premise was wrong, in a direction worth recording

Average precision was pre-declared as the primary metric on the grounds that
"positives are rare relative to the candidate pool", which would make AUROC
dominated by easy negatives. **On the phenotyped set positives are the
majority**: prevalence 0.731 on Tier A. The imbalance runs the other way, so AP
is inflated by a high positive rate rather than deflated by a rare one.

The metric was kept, because changing the primary metric after seeing which one
flatters the result is what pre-declaration exists to prevent. The consequence
is that **AP must be read against its floor**, and the lift is the number that
means anything.

### 5.2 Results, primary Tier A, position-blocks CV

No-skill floor (prevalence): **0.731**.

| Feature set | AP | Lift | P@95% recall | MCC @2× |
| --- | --- | --- | --- | --- |
| **Structure only (forest)** | **0.916** [0.853, 0.962] | **+0.185** | 0.790 | 0.15 |
| Structure only (ridge) | 0.905 [0.834, 0.957] | +0.174 | 0.798 | 0.15 |
| ESM-2 only | 0.848 | +0.116 | 0.728 | degenerate |
| **Frequency only** | **0.723** | **−0.008** | 0.768 | −0.06 |
| Conservation only | 0.623 | −0.108 | 0.728 | 0.07 |

**Frequency scores below its own no-skill floor.** So does conservation, more
heavily. "The models beat the frequency baseline" is true and nearly vacuous:
the baseline has no skill to beat. The reason is structural — only **13 of 119**
phenotyped Tier A variants are observed in circulating virus at all, so the
frequency features are near-constant on this set.

Under multiplicity correction across the 21-model grid, **only 2 models beat
frequency on average precision** (Holm-adjusted p < 0.05), and both are
structure-only. Phase 3, reporting Spearman without correction, had counted 8.
The protocol anticipated this in terms: "an uncorrected pairwise win on one
split is not evidence of superiority."

### 5.3 Where the structural signal actually lives

| Gene | n | Prevalence | Structure AP | Frequency AP |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 | 31 | 0.84 | **0.988** | 0.774 |
| HSV-2 UL23 | 12 | 0.67 | **1.000** | 0.647 |
| HSV-1 UL30 | 45 | 0.62 | 0.823 | 0.622 |
| HSV-1 UL5 | 7 | 0.57 | **0.567** | 0.525 |
| HSV-1 UL52 | 6 | 0.50 | **0.411** | 0.533 |

**The model works where the active site is resolved and fails where it is not.**
Thymidine kinase, with a well-resolved nucleoside pocket, is essentially solved.
Polymerase is intermediate. Helicase and primase — the two genes carrying
`site_resolved = 0` or near it — show no skill, and primase scores *below* its
prevalence floor and below frequency.

By drug mechanism class the same pattern appears: nucleoside analogues AP 0.916,
pyrophosphate analogue 0.875, nucleotide analogue 0.927, but
**helicase–primase inhibitors 0.582 against a 0.57 floor**.

This is mechanistically coherent rather than a curiosity: the feature is
distance to a bound ligand, and where no ligand defines the site the feature
measures nothing. It bounds the claim tightly, and it means the headline AP of
0.916 should never be quoted as a property of HSV resistance prediction in
general.

### 5.4 Good ranking, unusable decisions

At the pre-declared 2-fold operating point the best model's Matthews correlation
is **0.15**, and **9 of 22 models predict a single class for every variant** —
their MCC is undefined, not missing. Calibration slopes run 0.71–0.80, so the
models systematically under-disperse.

A model can rank variants well and still make no usable decision at the
threshold anyone would deploy. The protocol's insistence on reporting precision
at 95% recall and MCC *alongside* a ranking metric is what exposed this; had we
reported AUROC alone, as is common in this literature, it would have been
invisible.

### 5.5 External validation: cross-virus transfer

Training on HSV-2 (30 labels) and testing on HSV-1 (89) is the honest direction,
because HSV-2 structural features are inherited from HSV-1 geometry by
alignment — so the HSV-1 → HSV-2 direction tests a model against features
derived from its own training virus.

| Direction | Test n | Floor | Structure AP | Lift |
| --- | --- | --- | --- | --- |
| HSV-2 → HSV-1 | 89 | 0.685 | **0.904** | **+0.218** |
| HSV-1 → HSV-2 | 30 | 0.867 | 1.000 | +0.133 |

The honest direction holds up: trained on 30 labels from the other virus,
structure alone transfers at a lift larger than it achieves in
cross-validation. The full model does not — it collapses to AP 0.631, below the
floor, which is what overfitting 30 training rows across five feature groups
looks like.

## 6. The circularity audit

Reported here rather than in a supplement, because it is the project's main
methodological contribution.

### 6.1 Leave-one-gene-out has a negative null

Permuting the target **within** cross-validation groups — preserving each gene's
mean label, destroying only the feature–label link — gives a null centred at
**−0.15 to −0.28**, not zero.

The mechanism: the model predicts approximately the training mean for a
held-out gene, and removing a gene shifts that mean *away* from the gene's own
mean. Where genes differ in average resistance, and they do, the scheme
manufactures anti-correlation from noise.

Scored against zero, frequency-only (−0.153), conservation-only (−0.210) and
BLOSUM-only (−0.254) look "worse than chance". They are on the noise floor.
Position-blocks has a null near zero (+0.017, +0.016, −0.028, +0.002) and
therefore carries the headline results. **Any leave-one-gene-out figure quoted
against a zero null inverts the conclusion for the weaker feature sets.**

### 6.2 The matched-negative gap is zero

Resistance mutations are catalogued at conserved, buried, active-site
positions — and so is every other functionally important residue. A model that
has only learned "this position is conserved" would look predictive.

The protocol asks for identical models evaluated against random and against
conservation-matched negatives, with the difference reported as a headline
number. Each conservation-matched draw is paired here with a random draw of
**identical class counts**, so the only thing differing between arms is whether
the conservation distributions were equalised.

| Scheme | Model | Random AP | Matched AP | **Gap** |
| --- | --- | --- | --- | --- |
| Position blocks | Structure | 0.760 | 0.757 | **+0.003** |
| Position blocks | Structure (forest) | 0.758 | 0.770 | **−0.013** |
| Position blocks | Full | 0.690 | 0.698 | −0.008 |
| Position blocks | Full (forest) | 0.742 | 0.756 | −0.014 |
| Leave-one-gene-out | Structure | 0.733 | 0.732 | +0.001 |
| Leave-one-gene-out | Full (forest) | 0.690 | 0.690 | −0.000 |

**The gap is zero across all eight combinations**, range −0.015 to +0.003. The
structural signal is not a conservation proxy.

An earlier version of this analysis compared a matched subset of 64 rows against
the full 119 and concluded, hedged, that matching "bounds the conservation
contribution loosely rather than eliminating it". That hedge was an artefact of
comparing different sample sizes. With the size-matched control the answer is
clean.

### 6.3 ESM-2 does not survive residualisation

| Scheme | ESM-2 alone | Residualised against conservation + structure |
| --- | --- | --- |
| Position blocks | ρ +0.305 | **ρ +0.022** |
| Leave-one-gene-out | ρ +0.138 | **ρ −0.063** |

Conservation and structure linearly explain only 32% of the ESM-2 masked-marginal
score, so on a variance accounting two thirds of ESM-2 is something else. **That
two thirds carries no predictive signal.** Removing the explainable 32% destroys
essentially all of ESM-2's ability to rank fold change.

ESM-2 is therefore a legitimate *proxy* — it reaches AP 0.848 and clears its
permutation null at the floor — but it is not an independent line of evidence,
and the claim that a protein language model has learned something about
resistance beyond conservation and structural context is not supported here.

### 6.4 The feature correlation matrix explains the rest

Published in full in `data/processed/feature_correlation_matrix.tsv`. Every pair
at |ρ| ≥ 0.5:

| | | ρ |
| --- | --- | --- |
| `freq_count` | `freq_frequency` | +1.000 |
| `freq_count` | `freq_observed` | +0.998 |
| `freq_frequency` | `freq_observed` | +0.998 |
| **`cons_site_variability`** | **`freq_frequency`** | **+0.907** |
| `cons_site_variability` | `freq_count` | +0.907 |
| `cons_site_variability` | `freq_observed` | +0.905 |
| `esm_masked_marginal` | `esm_ref_logprob` | −0.889 |

Two things follow. **Conservation and frequency are very nearly the same
feature** (ρ = +0.907) — which is why they fail together and in the same
direction, and why "add conservation to frequency" never helped. And **no
structural feature appears in this list at all**: structure is uncorrelated
with conservation and frequency, which is precisely why it is the only group
that works and why §6.2 comes out at zero.

The three frequency features being mutually collinear at ≈1.0 also confirms
Aim 1's finding that the frequency signal is **one bit**: across all 241
variants with a numeric phenotype, continuous frequency gives ρ = −0.240 and
the binary observed/not-observed flag gives ρ = −0.236, both p = 0.0002.
Knowing how common an allele is adds nothing over knowing whether it has ever
been seen. Median fold change is 1.5× for observed variants and 5.0× for
unobserved.

### 6.5 Leakage

No evidence. Leakage would push a permutation null **above** zero; the maximum
null mean observed is +0.093 (position blocks) and −0.152 (leave-one-gene-out).

## 7. Discussion

The useful result of this project is not a model. It is a table.

Aim 1's reclassification table changes what is known about 644 catalogued
mutations by supplying the one thing the catalogue structurally cannot contain:
how often each substitution occurs in virus that was never under drug pressure.
Eleven entries turn out to describe the reference sequence rather than the
virus. Eight should be reclassified. Forty are common enough that a laboratory
will meet them and has nothing to go on.

Aim 4 is worth reporting mainly for what it rules out. Population frequency does
not predict phenotype — it scores below its own no-skill floor — and that is not
a failure of Aim 1 but a statement about disjoint populations: 120 of 644
catalogued mutations are seen in circulating virus, but only 13 of 119
phenotyped ones. The two aims barely overlap, and conflating them is the easiest
mistake available with these results.

What does predict phenotype is structural context, and only where the structure
resolves a site. That is a narrower claim than "machine learning predicts HSV
resistance", and it comes with a built-in test for where it applies: if no
ligand defines the active site, the feature measures nothing and the model
should not be used. Helicase and primase are exactly that case, and they are
also the genes where new drugs are arriving.

The methodological findings may outlast the biological ones. A grouped
cross-validation scheme can manufacture negative correlation from noise, so that
scoring against zero inverts conclusions. A conservation-matching control
compared against a different sample size produces a hedge that the data do not
require. A protein language model can clear every significance test and still
contribute nothing beyond features already in the model. And a model can reach
an average precision of 0.916 while making no usable decision at the threshold
it was declared against. Each of these was caught only because the protocol
pre-committed to a test that was inconvenient.

## 8. Limitations

**Only 13 of 119 primary phenotyped variants are observed in circulating
virus.** The frequency baseline is weak for structural reasons, not because
frequency is a poor idea, and Aim 4's comparison is correspondingly weak.

**119 labels across six virus–gene combinations, two with fewer than ten.** Any
cross-gene generalisation claim rests on six groups.

**Tier A labels pool across laboratories, assays and cell lines.**
`max_fold_change` is a maximum over drugs and studies, so it is an upper bound.
Between-study assay variation is unmodelled; a per-study random effect is the
honest treatment and is not implemented.

**UL52 has no resolved active site** and its distance feature is not comparable
to the other genes'. The breakdown in §5.3 shows the consequence directly.

**HSV-2 structural features are inherited from HSV-1** by alignment. Reasonable
at 70–80% identity, but an assumption, and it makes the HSV-1 → HSV-2 transfer
direction weaker evidence than it appears.

**Pairwise rather than codon-aware multiple alignment**, for want of an MSA
tool. Indels are handled less well than they should be.

**Gate 2's 100% pass rate is narrower than it looks.** It verifies that the
catalogue is internally consistent with the reference sequence, not that the
catalogue is faithful to its sources. Two substitutions passed this gate while
being flagged as misreported by an independent published audit.

**Two pre-declared validations were not run.** The temporal hold-out — training
on mutations characterised before a cut-off year and testing on those after,
which the protocol called "the single most persuasive validation in the plan" —
requires publication years that are recoverable offline for only 35 of 89
studies, and the recoverable ones cluster by publisher. A split on
partially-available year would confound era with journal. The positive-unlabelled
extension to variants of unknown significance was also not attempted.

**This is not a public pre-registration.** The analysis plan demonstrably
predates data acquisition by file timestamp, but it was not deposited with a
third party, and the evidence is under the authors' control.

## 9. Availability

| Artefact | Path |
| --- | --- |
| Frozen accession list, with fates | `data/release/accession_list.tsv` |
| Master variant table | `data/release/master_variant_table.tsv` |
| **Reclassification table** (deposited resource) | `data/release/reclassification_table.tsv` |
| Data dictionary | `data/release/DATA_DICTIONARY.md` |
| Feature correlation matrix | `data/processed/feature_correlation_matrix.tsv` |
| Pre-declaration and compliance audit | `reports/phase4/PRE_DECLARATION.md` |
| Pipeline | `scripts/`, driven by `scripts/run_all.py` |

Checksums for the release tables are in `reports/phase4/release_summary.json`.
