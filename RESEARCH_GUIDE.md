# A complete guide to this research

*What the study asks, how it was built, what it found, and why each choice was
made the way it was.*

This document is the explainer. The formal write-up is
[`reports/phase4/MANUSCRIPT.md`](reports/phase4/MANUSCRIPT.md); this is the
version that shows its reasoning — the design decisions, the things that went
wrong, the alternatives rejected and why. If you read one file to understand
the project, read this one.

**Contents**

1. [The short version](#1-the-short-version)
2. [The problem](#2-the-problem-a-catalogue-with-no-denominator)
3. [The four questions](#3-the-four-questions)
4. [The data](#4-the-data-and-how-1143-records-became-599-isolates)
5. [Methods, and why each one](#5-methods-and-why-each-one)
6. [How the experiment runs](#6-how-the-experiment-runs)
7. [The gates](#7-the-gates)
8. [Results](#8-results)
9. [Benchmarks](#9-benchmarks)
10. [The circularity audit](#10-the-circularity-audit)
11. [External validation](#11-external-validation)
12. [What we deliberately did not do](#12-what-we-deliberately-did-not-do)
13. [Limitations](#13-limitations)
14. [How the numbers are kept honest](#14-how-the-numbers-are-kept-honest)
15. [Reproducing it](#15-reproducing-it)
16. [Where everything lives](#16-where-everything-lives)
17. [Glossary](#17-glossary)

---

## 1. The short version

Herpes simplex virus has catalogues of mutations said to cause antiviral
resistance. Those catalogues are built from case reports: someone sequenced a
virus from a patient whose treatment failed, found a mutation, and published
it. What nobody had systematically checked is how often those same mutations
occur in ordinary virus from people who were never treated.

That number is the denominator, and without it a resistance mutation and a
harmless natural variant look identical on paper.

We computed it for all 644 catalogued substitutions in the four genes that
antiviral drugs target, using 599 publicly deposited HSV genomes. The headline:

| Finding | Number |
| --- | --- |
| Catalogued substitutions examined | 644 |
| Never seen in circulating virus | 516 |
| Seen in circulating virus | 128 |
| **Are the *majority* allele** — the reference genome carries the rarer residue | **11** |
| **Recommended for reclassification as polymorphisms** | **7** |
| **Common in circulation but never phenotyped** — the priority list for lab work | **41** |
| Genuinely common *and* strongly resistant | **0** |

That last row is the reassuring one. Among the 241 substitutions with a
measured resistance value, nothing above 5% frequency exceeds 2.4-fold
resistance. There is no common, strongly resistant allele hiding in circulating
HSV — but you could not have known that without building the denominator.

---

## 2. The problem: a catalogue with no denominator

Imagine a list of "mutations associated with resistance". Every entry means
the same thing: *this substitution was present in a virus that resisted a
drug*. The entry is evidence of an association. It is not evidence of a cause,
and it says nothing about how common the substitution is anyway.

Suppose a substitution is carried by 60% of all HSV worldwide. A patient whose
treatment fails has a 60% chance of carrying it for reasons that have nothing
to do with their treatment failing. Sequence enough failures and it will appear
in a resistance catalogue — not because it causes resistance, but because it is
*everywhere*.

This is not hypothetical. **Eleven catalogued HSV "resistance mutations" are the
commonest residue at their position in circulating virus.** For those, the
catalogue has recorded a property of the reference sequence rather than a
property of the virus. The reference genome happens to carry the minority
residue, so everything else looks like a mutation.

**Why the gap persisted.** The two halves of the question live in different
places. Resistance phenotypes are in the primary literature and curated
databases. Population allele frequencies require downloading, aligning and
de-duplicating hundreds of public genomes. Nobody was being careless — the
denominator is simply a separate piece of work, and this project is that work.

---

## 3. The four questions

| Aim | Question | Result |
| --- | --- | --- |
| **1** | What does the catalogue look like once every entry carries a population frequency? | 7 reclassifications, 41 phenotyping priorities, 11 reference artefacts |
| **2** | Is thymidine kinase under weaker evolutionary constraint than DNA polymerase? | Yes, decisively |
| **3** | How much helicase–primase resistance exists *before* the new drugs are used? | Essentially none |
| **4** | Can sequence and structure predict resistance better than frequency alone? | Yes, but only where a structure resolves the drug-binding site |

Aims 1–3 are descriptive measurements. Aim 4 is a modelling question, and it is
the one where it is easiest to fool yourself — which is why roughly a third of
the project is spent trying to.

---

## 4. The data, and how 1,143 records became 599 isolates

### Sources, frozen at a point in time

| Source | What | Pinned how |
| --- | --- | --- |
| **HerpesDRG** | The resistance mutation catalogue | Commit `22fc186`, MD5 recorded, active entries only (`status == "A"`) |
| **GenBank** | HSV-1 (taxon 10298) and HSV-2 (10310) genomes, 140–200 kb | Frozen query, timestamped to the second |
| **PDB** | Protein structures: `1vtk`, `8exx`, `9vlq`, `9nnp` | Auto-selected by alignment coverage |
| **Crossref** | Publication years, for the temporal validation | Cached with retrieval date |

Everything is pinned because a result you cannot re-derive is not a result. The
GenBank query returned 1,143 accessions, and **every one is listed in
`data/release/accession_list.tsv` with its fate** — so a reader can tell an
accession that was never retrieved from one that was retrieved and excluded.
That distinction is usually invisible in published work.

### The funnel

| Fate | n | Why |
| --- | --- | --- |
| **Analysed** | **599** | |
| Duplicate of an isolate already counted | 187 | Identical sequence or strain |
| No target-gene CDS extracted | 168 | Annotation did not yield a usable gene |
| Failed CDS quality control | 127 | Length, frame, ambiguity or divergence |
| Laboratory-adapted strain | 30 | Not circulating virus |
| Serial isolate from an already-represented patient | 18 | Would double-count one person |
| Non-human host | 14 | Not circulating human virus |

### De-duplication, which matters more than it sounds

Three levels, removing 281, 5 and 18 records respectively:

1. **Identical sequence** across all four target genes
2. **Same strain name**
3. **Same patient** — serial isolates taken from one person over time

The third is small in count and large in consequence. A patient with a
persistent infection may be sampled a dozen times; depositing all twelve is
good practice, but counting them as twelve independent observations inflates
whatever that patient happened to carry.

> **This single step changed the project's most dramatic finding.** One
> substitution, UL23 R222C, appeared at 7.6% frequency with 98-fold measured
> resistance — which would have been a genuine surveillance alarm. Collapsing
> serial isolates to one per patient took it from 19 carriers to 11, and 7.6%
> to 4.6%. It is not a circulating resistance allele. See §8.1.

---

## 5. Methods, and why each one

This section is the heart of the explainer. Every choice below had a plausible
alternative, and the reason for rejecting it matters as much as the choice.

### 5.1 Pairwise alignment, not a multiple-sequence alignment

**What we did.** Reference-anchored pairwise alignment with Biopython's
`PairwiseAligner` (BLOSUM62, free terminal gaps).

**Why.** No MSA tool (MAFFT, MUSCLE, Clustal) was available in the environment.
For coding sequences at 70–80% identity to a known reference, pairwise
alignment is adequate, and it has one real advantage: every position maps
directly onto reference coordinates with no column-mapping step, so a catalogue
entry naming "position 222" needs no translation.

**The cost, stated plainly.** Indels are handled less well than a codon-aware
MSA would handle them. This is a genuine limitation and is declared as one
rather than hidden.

### 5.2 pN/pS against the population consensus, not the reference

**What we did.** Measured selection relative to the consensus sequence of the
population itself.

**Why.** Using the reference would count every position where the reference
carries the minority residue as "polymorphism" — inflating the estimate by a
quantity that reflects which genome someone picked as reference decades ago,
not selection. Since this project's *first finding* is that 11 catalogued
positions have a minority-residue reference, using the reference here would
have imported exactly the artefact we set out to measure.

### 5.3 Every frequency carries its ascertainment

A frequency alone is not a claim you can act on. Each one is published with:

| Column | What it tells you |
| --- | --- |
| `n_studies_carrying` | How many independent submissions contributed a carrier |
| `max_study_share` | The largest share of carriers from a single submission |
| `loso_min_freq` | The frequency if you drop the most influential study |
| `study_dominated` | Whether dropping one study halves the frequency |

**Why this is not supplementary.** A 5% allele drawn entirely from one
submission and a 5% allele seen across five submissions are different claims
about the world. **60 of the 128 observed variants are study-dominated** —
roughly half. Any downstream use of these frequencies that ignores that column
will overstate about half of them.

### 5.4 The active site is defined by ligands, not by a list

**What we did.** Structural features measure distance to whatever ligand is
actually present in the crystal structure, and solvent accessibility.

**Why not a list of known catalytic residues?** Because Aim 4 asks whether
structural context predicts resistance. Hand-writing a list of "the important
residues" would import the published claim about where resistance mutations sit
— and then the model would confirm it. That is circular.

**It validated itself.** Without any residue being named in advance, the
thymidine-binding pocket came out buried and in contact range: E83 at 2.79 Å,
Y101 2.67 Å, Q125 2.85 Å, R163 2.45 Å, R222 3.27 Å, all with relative solvent
accessibility near zero. The geometry found the pocket on its own.

**And it found its own boundary.** UL52 (primase) has *no* resolved site —
zero residues within 5 Å of any ligand, median distance 42 Å, its nearest
ligand belonging to a neighbouring protein. Those rows are flagged
`site_resolved = 0`. This turns out to explain where the whole model works and
where it fails (§8.4).

### 5.5 Gene identity is excluded as a feature

UL30 contributes 45 labels at one base rate; UL52 contributes 6 at another.
Giving a model the gene name lets it recover label prevalence directly instead
of learning anything biological. Gene is used as a *grouping key* for
cross-validation and never as a predictor.

### 5.6 ESM-2 zero-shot, and nothing deeper

**What we did.** ESM-2 (`esm2_t33_650M_UR50D`) masked marginals, used zero-shot.

**Why not fine-tune, or use a bigger model?** With a few hundred labelled
positives you cannot fit a deep network and you certainly cannot validate one.
This was ruled out *in advance*, in writing, rather than after seeing that it
would not work.

### 5.7 Grouped cross-validation, two schemes

Random train/test splits would be meaningless here: adjacent residues in the
same gene are not independent — they share a structural context, a conservation
column, and often a source study. Two grouped schemes are used, and the unit
being split is always one amino-acid substitution.

**Leave-one-gene-out.** The group key is `virus:gene`, so the 119 Tier A
variants fall into six groups: HSV-1 UL30 (45), HSV-1 UL23 (31), HSV-2 UL30
(18), HSV-2 UL23 (12), HSV-1 UL5 (7) and HSV-1 UL52 (6). Six folds; each trains
on five genes and tests on the sixth. Nothing from the held-out gene is in
training.

**Leave-one-block-out.** The group key is `virus:gene` plus the position
divided by 50, giving 44 blocks of 1 to 13 variants. Five folds, blocks
assigned whole, so roughly 95 variants train and 24 test. The same gene appears
on both sides, but never the same 50-residue neighbourhood.

**Hyperparameters never see the test fold.** Inside each training fold, an
inner grouped split over that fold's own groups picks between ridge penalties
and forest settings, scored by Spearman correlation. The winner is refit on the
whole training fold and predicts the test fold once. Every fold's predictions
go into one out-of-fold vector of 119, and that vector is what gets scored —
not any within-fold number.

Running both schemes is not belt-and-braces. It is what revealed that one of
them is broken (§10.1).

### 5.8 Where the mathematics is written out

This guide explains the reasoning; the manuscript carries the formulae. Every
estimator and test is defined in full there, so nothing depends on a convention
the reader has to guess:

| What | Manuscript section |
| --- | --- |
| How the fold change is formed, and the drug it may come from | §2.3 |
| Wilson score interval, and why not Wald | §2.5 |
| Nei–Gojobori site counting, pathway averaging, Jukes–Cantor | §2.6 |
| Relative solvent accessibility, ESM-2 masked marginal | §2.7 |
| The two cross-validation schemes, stated completely | §2.7 |
| Average precision, P@95, MCC, calibration, Spearman | §2.8 |
| Bootstrap intervals, paired Δ-AP, two-sided *p* | §2.8 |
| Holm adjustment, permutation null, matching, residualisation | §2.8 |

Two of those deserve a note here, because they are the ones most often left
implicit. The pathway-averaging rule in §2.6 is what makes Nei–Gojobori a
counting method rather than a model fit, and it is not neutral for codons
differing at more than one position. And the average precision used throughout
is the step-wise sum over positives, not the area under an interpolated
precision–recall curve; at this sample size the two disagree in the third
decimal.

The manuscript's §11 is a reference list. Its HSV phenotype citations are
labelled *via HerpesDRG*, because each was taken from the frozen catalogue's
own `ref_doi` field for the variant under discussion rather than from an
independent literature search. Citing them asserts what the catalogue
attributes to each work — which is the claim the paper makes — and not that the
primary measurements were re-examined.

### 5.9 Average precision, always beside its floor

**The primary metric was declared before any data was acquired**, along with
the 2-fold resistance threshold and the whole metric hierarchy
(`reports/phase4/PRE_DECLARATION.md`).

**The declared justification turned out to be wrong**, and this is recorded
rather than quietly fixed. The protocol chose average precision because
"positives are rare relative to the candidate pool", which would make AUROC
dominated by easy negatives. In fact, on the phenotyped set **positives are the
majority** — prevalence 0.731. The imbalance runs the other way.

**The metric was kept anyway.** Changing your primary metric after seeing which
one flatters your result is precisely what pre-declaration exists to prevent.
The correct response is to report AP beside its no-skill floor and say the
justification does not hold. An AP of 0.916 against a floor of 0.731 is a lift
of **+0.185**, not a 0.916-sized achievement.

---

## 6. How the experiment runs

The whole study is 21 scripted stages, driven by a single runner:

```bash
python scripts/run_all.py
```

| Phase | What happens |
| --- | --- |
| **0** | Pin the catalogue, assign evidence tiers, Gate 1 |
| **1** | Fetch genomes, verify coordinates (Gate 2), extract genes, QC, de-duplicate, compute frequencies |
| **2** | Reclassification table, selection analysis, helicase–primase baseline (Gate 3) |
| **3** | ESM-2 scores, structural features, assemble the modelling set, baselines, circularity audit, Gate 4 |
| **4** | Pre-declared metrics, matched-negative control, release tables, cross-virus and temporal validation, figures, verification gate |

**Determinism.** Everything runs from `SEED = 20261102`. Two stages are not
pure — the GenBank fetch needs the network and the ESM-2 scoring wants a GPU —
and both are skipped when their outputs already exist. Since those outputs are
committed, **every number downstream of them reproduces with no network
connection and no graphics card**.

The runner stops at the first failure rather than carrying on with stale
inputs.

---

## 7. The gates

Four pre-declared checkpoints, each with a threshold fixed before it was
evaluated. A gate exists so that the project can fail honestly at a known
point rather than drift toward whatever result the data happens to allow.

| Gate | Asks | Threshold | Outcome |
| --- | --- | --- | --- |
| **1** | Are there enough high-quality labels to attempt this? | ≥ 50 Tier A HSV-1 mutations | **PASS** — 103 |
| **2** | Does the catalogue agree with the reference sequence? | Coordinate verification | **PASS** — 100% agreement |
| **3** | Does any labelled variant exceed 1% in untreated isolates? | At least one | **PASS** — 62 |
| **4** | Do the models beat the frequency baseline? | Two-part criterion | **PASS** — in all four scheme/set combinations |

> **Gate 2's 100% is narrower than it looks**, and the project says so. It
> verifies that the catalogue is internally consistent with the reference
> sequence — not that the catalogue is faithful to its original sources. Two
> substitutions passed this gate while being queried by an external review of
> this work as possibly misreported relative to their primary sources. Both are
> settled in the manuscript's §3.1, against the catalogue's own per-drug
> records.

---

## 8. Results

### 8.1 Aim 1 — the catalogue with a denominator

Of 644 catalogued substitutions, **128 occur in circulating virus**.

| Band | n | Meaning |
| --- | --- | --- |
| Not observed | 516 | No population evidence either way |
| Singleton | 45 | Exactly one carrier |
| Rare (<1%) | 21 | |
| Uncommon (1–10%) | 36 | |
| **Common (10–50%)** | **15** | |
| **Majority (>50%)** | **11** | The reference carries the minority residue |

Each entry also receives an interpretive flag:

| Flag | n | What it means |
| --- | --- | --- |
| `consistent_rare` | 60 | Rare, as a resistance determinant should be |
| **`unphenotyped_but_common`** | **41** | Above 1% in circulation, never phenotyped |
| `reference_artifact` | 10 | Majority allele, no measured resistance |
| **`reclassify_as_polymorphism`** | **7** | Catalogued as resistance-associated, above 1%, fold change ≤ 2× |
| `consistent_with_resistance` | 6 | Rare *and* phenotypically resistant |
| `circulating_resistance_allele` | 3 | Common and resistant — all three fail on inspection |
| `majority_allele_with_measured_resistance` | 1 | HSV-2 UL23 G39E, 79% frequency, 2.4× |

**The 7 reclassifications are the concrete deliverable. The 41 unphenotyped
common variants are the clearest call to action** — each is carried by more
than 1% of circulating virus, so a clinical laboratory will encounter them, and
there is no phenotype to interpret them with.

#### The three that looked like circulating resistance

| Variant | Freq | Carriers | Studies | `max_study_share` | `loso_min` | Fold change | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| UL23 R222C | 0.046 | 11 | 3 | 0.73 | 0.013 | 98× | Ascertainment artefact |
| UL23 R176Q | 0.017 | 4 | 1 | 1.00 | 0.0 | 200× | Ascertainment artefact |
| UL23 R41H | 0.051 | 12 | 5 | 0.33 | 0.038 | 2.3× | Sound frequency, phenotype does not hold up |

**R222C** is the one that would have mattered: 98-fold resistance at an
apparent 7.6% frequency. Patient-level collapse took it to 11 carriers and
4.6%, and 8 of those 11 come from a single submission — dropping it leaves
1.4%.

The obvious objection is that de-duplication threw away real carriers, so
every record in the frozen query carrying R222C was re-examined. There are 38;
11 are counted. Of the other 27, **ten are serial isolates from four patients
already represented** (five from one patient, with strain names like `10-6-1`,
`10-7-1`, `10-11-3`) and **seventeen are sequence-identical across all four
target genes** to a record already kept. None is independent. Even counting all
38 would leave 63% of them from the one dominant submission.

**R41H** is the interesting case because it is *not* an artefact — 12 carriers
across five submissions, robust to study dropout. What disqualifies it is the
phenotype, and the catalogue's own records say so once they are read drug by
drug. The 2.3-fold that put it in the table is a **foscarnet** value, and
foscarnet inhibits the DNA polymerase directly without needing thymidine kinase
to activate it — so a UL23 substitution has no mechanism by which to shift it.
For aciclovir, the drug UL23 does activate, three independent studies give 0.8
to 1.0-fold. This is the cost of taking the fold change as the maximum over all
drugs tested, which the manuscript declares in §2.3 and works through in §3.1.
**Frequency and phenotype each need their own scrutiny**; the ascertainment flag
alone would have retained it.

#### The empty region

Among the 241 substitutions with a numeric phenotype, 30 are observed in
circulating virus and 211 are not. The two groups differ: median fold change
**1.5× for observed**, **5.0× for unobserved**.

Above 5% frequency, the strongest phenotype in the entire catalogue is
**2.4×**. Two substitutions clear the 2-fold threshold at that frequency, both
marginally, and **none reaches 5×**. The informative corner of the plot — common
*and* strongly resistant — is empty.

### 8.2 Aim 2 — selection

| Gene | pN/pS | 95% CI | Synonymous sites |
| --- | --- | --- | --- |
| **HSV-1 UL23** (thymidine kinase) | **0.371** | 0.338–0.407 | 620 |
| **HSV-1 UL30** (DNA polymerase) | **0.173** | 0.161–0.185 | 1,536 |
| HSV-1 UL5 (helicase) | 0.181 | 0.159–0.205 | 955 |
| HSV-1 UL52 (primase) | 0.140 | 0.128–0.154 | 2,014 |

Difference UL23 − UL30: **+0.1980 [+0.1664, +0.2361]**, with UL23 higher in
**100% of 1,000 bootstrap resamples**.

Thymidine kinase is under markedly weaker purifying selection than polymerase.
This is the expected direction — TK is dispensable for replication in culture,
polymerase is not — and it supplies a population-genetic reason why TK
accumulates the substitutions that then populate resistance catalogues.

**HSV-2 UL23 is withheld rather than reported.** It has only 32 synonymous
sites, giving a point estimate of 2.46 with an interval spanning 1.74–3.67.
That is an artefact of a tiny denominator, not evidence of diversifying
selection, and publishing it with a wide interval would invite exactly the
wrong reading.

### 8.3 Aim 3 — the pre-therapeutic baseline

| Gene | Isolates | Carrying a catalogued substitution |
| --- | --- | --- |
| HSV-1 UL5 | 251 | **0** |
| HSV-1 UL52 | 299 | **0** |
| HSV-2 UL5 | 247 | 2 (0.8%) |
| HSV-2 UL52 | 225 | **0** |

Pritelivir and amenamevir — helicase–primase inhibitors — are new enough that
this is a genuine *pre-therapeutic* baseline. There is essentially no
pre-existing helicase–primase resistance in circulating HSV, and any future
rise is now measurable against a recorded starting point.

**The honest caveat.** The catalogue lists only 8 HSV-1 UL5 and 6 UL52
positions, so the test has little to find. A permutation test of whether
catalogued positions are more conserved than random ones is non-significant for
both (p = 0.54, p = 0.37). For HSV-2 UL5 the catalogued positions come out 2.2×
*more* variable than gene-wide — the opposite of constraint — but at p = 0.052
on four positions, which is a number to report and not to interpret.

### 8.4 Aim 4 — where prediction works, and where it stops

The headline is not "machine learning predicts HSV resistance". It is narrower
and more useful:

> **Structural context predicts measured resistance — but only in genes where a
> crystal structure actually resolves the drug-binding site.**

| Gene | n | Prevalence (floor) | Structure AP | Frequency AP |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 | 31 | 0.84 | **0.988** | 0.774 |
| HSV-2 UL23 | 12 | 0.67 | **1.000** | 0.647 |
| HSV-1 UL30 | 45 | 0.62 | 0.823 | 0.622 |
| HSV-1 UL5 | 7 | 0.57 | 0.567 | 0.525 |
| HSV-1 UL52 | 6 | 0.50 | **0.411** | 0.533 |
| HSV-2 UL30 | 18 | 1.00 | undefined | undefined |

Thymidine kinase, with a well-resolved nucleoside pocket, is essentially
solved. Polymerase is intermediate. Helicase and primase — the genes flagged
`site_resolved = 0` — show no skill at all, and primase scores *below* its own
prevalence floor.

*(HSV-2 UL30 has 18 variants, all resistant. With only one class present,
average precision is undefined — it is shown rather than dropped so the rows
sum to the 119 used elsewhere.)*

The same pattern appears by drug class: nucleoside analogues 0.916,
pyrophosphate analogue 0.875, nucleotide analogue 0.927 — but **helicase–primase
inhibitors 0.582 against a 0.57 floor**.

This is mechanistically coherent rather than a curiosity. The feature is
*distance to a bound ligand*. Where no ligand defines the site, the feature
measures nothing. **The claim comes with a built-in test for where it applies**
— and helicase and primase, where it does not apply, are exactly the genes
where the new drugs are arriving.

#### Good ranking, unusable decisions

At the pre-declared 2-fold operating point, the best model's Matthews
correlation is **0.15**, and **8 of 22 models predict a single class for every
variant**.

Calibration is not the problem — the three surviving models have slopes of
0.77, 0.95 and 1.20. The failing feature sets fail by predicting almost nothing:
frequency's calibration slope is 0.06, conservation's 0.17. That is a
near-constant prediction dressed as a model.

> **A model can rank variants well and still make no usable decision at the
> threshold anyone would deploy.** Requiring precision-at-95%-recall and MCC
> *alongside* a ranking metric is what exposed this. Under the AUROC-only
> reporting common in this literature, it would have been invisible.

---

## 9. Benchmarks

Primary Tier A set (119 variants), position-blocks cross-validation.
**No-skill floor (prevalence): 0.731.** Lift is the only number that means
anything.

| Feature set | AP | **Lift** | P@95% recall | MCC @2× |
| --- | --- | --- | --- | --- |
| **Structure only (forest)** | **0.916** [0.853, 0.962] | **+0.185** | 0.790 | 0.15 |
| Structure only (ridge) | 0.905 [0.834, 0.957] | +0.174 | 0.798 | 0.15 |
| ESM-2 only | 0.848 | +0.116 | 0.728 | degenerate |
| BLOSUM62 only (forest) | 0.802 | +0.071 | 0.722 | −0.10 |
| Amino-acid properties only | 0.751 | +0.020 | 0.735 | 0.02 |
| **Frequency only** | **0.723** | **−0.008** | 0.762 | −0.08 |
| Conservation only | 0.717 | −0.014 | 0.741 | 0.07 |
| BLOSUM62 only (ridge) | 0.715 | −0.016 | 0.735 | degenerate |

**Read the bottom of that table first.** Frequency scores *below its own
no-skill floor*, and so does conservation. "The models beat the frequency
baseline" is therefore true and nearly vacuous — the baseline has no skill to
beat.

The reason is structural: **only 14 of the 119 phenotyped Tier A variants are
observed in circulating virus at all**, so the frequency features are
near-constant on this set. Aims 1 and 4 barely overlap, and conflating them is
the easiest available mistake with these results.

**After multiplicity correction** across the 21 comparisons (every model except
the frequency baseline, which is the thing being compared against), **only 3
models beat frequency** at Holm-adjusted p < 0.05: the two structure-only
variants and the full model. Seven clear p < 0.05 uncorrected — so the
correction is the difference between a reportable result and a grid of
near-misses.

---

## 10. The circularity audit

This is the project's main methodological contribution, and it is in the main
text rather than a supplement on purpose. The question throughout: *is the
model learning biology, or learning an artefact of how the data was assembled?*

### 10.1 One cross-validation scheme is broken

Permuting the target *within* cross-validation groups — preserving each gene's
mean label, destroying only the feature–label link — should give a null
centred at zero.

For leave-one-gene-out it gives **−0.16 to −0.29**.

**The mechanism.** The model predicts roughly the training mean for a held-out
gene. Removing that gene shifts the training mean *away* from the gene's own
mean. Where genes differ in average resistance — and they do — the scheme
manufactures anti-correlation out of pure noise.

**The consequence.** Scored against zero, frequency-only (−0.150),
conservation-only (−0.228) and BLOSUM-only (−0.254) look "worse than chance".
They are not. They are sitting on the noise floor.

Position-blocks is better behaved and carries the headline results — but only
for feature sets that have signal. Its null sits between −0.014 and +0.079 for
structure, ESM-2 and the full model, but is **−0.167 for frequency and −0.115
for conservation**. So no scheme here licenses scoring those two against zero.

> **Any grouped-CV figure quoted against a zero null inverts the conclusion for
> the weaker feature sets** — and the weaker feature sets are exactly the ones a
> reader wants to compare against.

### 10.2 Is the model just detecting "conserved"?

Resistance mutations sit at conserved, buried, active-site positions. So does
every other functionally important residue. A model that learned only "this
position is conserved" would look predictive and mean nothing.

The test: run identical models against random negatives and against
conservation-matched negatives, and report the difference.

**The gap is zero across all eight model/scheme combinations** — range −0.016
to +0.003. The structural signal is not a conservation proxy.

**Why the control's design matters.** Conservation matching necessarily
discards rows: 64 of the 119 survive it. Comparing that matched arm against the
*full* 119 differs in sample size as well as in conservation, and average
precision is sensitive enough to *n* at this scale that such a comparison
returns a gap of ambiguous origin. Pairing each matched draw with a random draw
of **identical size and class balance** removes that confound. Only then does a
gap of zero actually mean what it appears to mean.

### 10.3 ESM-2 does not survive residualisation

| Scheme | ESM-2 alone | Residualised against conservation + structure |
| --- | --- | --- |
| Position blocks | ρ +0.305 | **ρ −0.011** |
| Leave-one-gene-out | ρ +0.138 | **ρ −0.113** |

Conservation and structure linearly explain only 34% of the ESM-2 score — so on
a variance accounting, two thirds of ESM-2 is something else. **That two thirds
carries no predictive signal.** Removing the explainable third does not merely
weaken ESM-2's ranking; it eliminates it.

ESM-2 is a legitimate *proxy* — it reaches AP 0.848 and clears its permutation
null. But it is not an independent line of evidence, and the claim that a
protein language model has learned something about resistance beyond
conservation and structural context **is not supported here**.

### 10.4 The correlation matrix explains the rest

Published in full. Every pair at |ρ| ≥ 0.5 is listed; two facts matter:

- **Conservation and frequency are very nearly the same feature** (ρ = +0.888).
  That is why they fail together, in the same direction, and why adding
  conservation to frequency never helped.
- **No structural feature appears in that list at all.** Structure is
  uncorrelated with conservation and frequency — which is precisely why it is
  the only group that works, and why §10.2 comes out at zero.

The frequency signal is also **one bit**. Across all 241 variants with a
numeric phenotype, continuous frequency gives ρ = −0.248 and the binary
*seen/not-seen* flag gives ρ = −0.243. Knowing *how common* an allele is adds
essentially nothing over knowing *whether it has ever been seen*.

### 10.5 Leakage

No evidence. Leakage would push a permutation null *above* zero; the maximum
null mean observed is +0.079.

---

## 11. External validation

Two tests the models were never tuned against.

### Cross-virus transfer

| Direction | Test n | Floor | Structure AP | Lift |
| --- | --- | --- | --- | --- |
| **HSV-2 → HSV-1** | 89 | 0.685 | **0.904** | **+0.218** |
| HSV-1 → HSV-2 | 30 | 0.867 | 1.000 | +0.133 |

HSV-2 → HSV-1 is the honest direction, because HSV-2 structural features are
inherited from HSV-1 geometry by alignment — so the reverse direction would
test a model against features derived from its own training virus.

Trained on just 30 labels from the other virus, structure alone transfers at a
*larger* lift than it achieves in cross-validation. The full model does not —
it collapses to AP 0.632, below the floor, which is what overfitting 30
training rows across five feature groups looks like.

### The temporal hold-out

The strongest declared test: train on mutations characterised before a cut-off
year, test on those characterised after. It approximates genuine prospective
validation — *knowing only what was published by a given date, can the model
identify what was discovered next?*

All 89 distinct study DOIs resolved to publication years via Crossref, spanning
1987–2025 across 10 publishers, so the split is not confounded with journal.

| Cut-off | Train → test | Structure | Frequency | Conservation | ESM-2 | Full |
| --- | --- | --- | --- | --- | --- | --- |
| 2005 | 41 → 78 | **+0.195** | −0.077 | +0.085 | +0.118 | +0.110 |
| 2008 | 47 → 72 | **+0.215** | +0.019 | +0.100 | +0.052 | +0.126 |
| 2010 | 52 → 67 | **+0.218** | +0.020 | +0.062 | −0.011 | +0.123 |
| **2012** | **59 → 60** | **+0.229** | +0.030 | +0.042 | +0.050 | +0.090 |
| 2015 | 70 → 49 | **+0.170** | −0.006 | +0.162 | +0.101 | +0.066 |

**Structure wins at every cut-off, on both label sets, with no alternative
reaching +0.163 at any split.** Nothing else is stable — conservation ranges
from +0.042 to +0.162, ESM-2 from −0.011 to +0.118, and none holds a consistent
rank against the others. Frequency straddles zero throughout.

This is the most persuasive single result in the project: five independent
prospective splits, no tuning against any of them, same winner each time.

**The cut-off year was not pre-declared** — the one place in the project where
an analysis parameter was chosen after seeing data. The mitigation is that it
was not chosen on *performance*: all five cut-offs are reported, and 2012 is
nominated because it splits the set most evenly, which is outcome-independent.

> **A finding about the literature, not the models.** Resistance prevalence
> falls sharply from train to test at every cut-off — 0.86 → 0.60 at 2012.
> Early resistance work characterised almost exclusively resistant variants;
> later work reports many more negatives. The test sets are systematically
> *harder* than the training sets, which makes a stable positive lift more
> impressive — and means a model trained on older literature is calibrated to a
> base rate that no longer holds.

---

## 12. What we deliberately did not do

Recording what was skipped, and why, is part of the method.

### Positive-unlabelled learning — declared, then declined

This was pre-declared as the extension where the clinical value lies, and that
judgement was right: a laboratory meeting a substitution in a patient isolate is
far more likely to be looking at an *unphenotyped* variant than a characterised
one. There are **717 observed variants with no measured fold change**, 619 of
them not catalogued at all.

It was not attempted, for three reasons:

1. **The standard PU assumption is violated by this project's own subject
   matter.** PU methods assume positives are selected at random from all
   positives. Here, a variant is catalogued *because* someone found it in a
   resistant isolate. Ascertainment bias is exactly what this project
   measures — applying a method that assumes it away, inside a paper about it,
   would be incoherent.
2. **The class prior cannot be pinned down** at 87 positives, and an unstable
   ranking presented as a triage list is worse than no list.
3. **The one feature group that works does not cover the pool.** The structural
   signal is confined to genes with a resolved site; the unphenotyped pool is
   concentrated in UL30 and UL52, and only 317 of the 717 have a resolved site.
   A score over UL52 variants would be close to noise while looking
   authoritative.

**What is delivered instead**, with no modelling assumptions at all: the 41
`unphenotyped_but_common` variants — a phenotyping priority list derived from
observation rather than from a model.

### Also skipped

| Not done | Why |
| --- | --- |
| Brier score / calibration curve | The declared primary task is regression, so models emit fold-change estimates, not probabilities. Regression calibration is reported instead. |
| SHAP / permutation importance | Made conditional in the protocol. With two structural features carrying the result, the correlation matrix and staged comparison already answer what it would. |
| Fine-tuned language models | Ruled out in advance — unfittable and unvalidatable at this label count. |

---

## 13. Limitations

Stated as the paper states them, without softening.

- **Only 14 of 119 phenotyped Tier A variants are observed in circulating
  virus.** The frequency baseline is weak for structural reasons, and Aim 4's
  central comparison is correspondingly weak.
- **119 labels across six virus–gene combinations, two with fewer than ten.**
  Any cross-gene generalisation rests on six groups.
- **Labels pool across laboratories, assays and cell lines.** `max_fold_change`
  is a maximum over drugs and studies — an upper bound. Between-study assay
  variation is unmodelled; a per-study random effect is the honest treatment
  and is not implemented.
- **UL52 has no resolved active site**, so its distance feature is not
  comparable to the other genes'.
- **HSV-2 structural features are inherited from HSV-1** by alignment —
  reasonable at 70–80% identity, but an assumption.
- **Pairwise rather than codon-aware multiple alignment.**
- **Gate 2's 100% pass rate is narrower than it looks** — internal consistency,
  not fidelity to sources.
- **The temporal cut-off year was not pre-declared.**
- **This is not a public pre-registration.** The analysis plan demonstrably
  predates data acquisition by file timestamp, but it was not deposited with a
  third party, and the evidence is under the authors' control. The protections
  a real pre-registration offers are *partially* in place here, not fully.

> **No variant score in this work should be used to interpret a patient
> isolate.**

---

## 14. How the numbers are kept honest

Every figure quoted in the manuscript, README, pre-declaration and data
dictionary is re-read from its source summary and compared, by
`scripts/p4_05_verify.py`. It currently runs **259 checks** and the pipeline
fails if any disagrees.

This is not decoration. The gate has caught real errors that survived multiple
readings:

| Error | How it hid |
| --- | --- |
| Tier counts of 133/89/478 summing to 700 in a section about 644 substitutions | Tier A is 133 under *both* partitions — only the sum gave it away |
| Frequency bands labelled 1–5% and 5–50% when the code cuts at 1% and 10% | The *counts* were right, so every count check passed |
| 187 isolates described as a "coverage or divergence filter" | They are de-duplication removals — the label said the data were bad when they were redundant |
| A gene breakdown summing to 101 against a 119 used everywhere else | The missing group has only one class, so it was silently dropped |

The checks now cover partitions summing to their totals, breakdowns summing to
their sets, and the counts behind each named claim — the class of error that a
pipeline comparing JSON to prose structurally cannot catch on its own.

---

## 15. Reproducing it

```bash
python scripts/run_all.py --list      # show every stage and what it needs
```

```bash
python scripts/run_all.py             # run everything
```

```bash
python scripts/p4_05_verify.py        # just the verification gate
```

The network and GPU stages are skipped when their outputs exist, and those
outputs are committed — so the full analysis reproduces offline on a laptop.
Expect roughly 100 minutes for the three slow modelling stages.

---

## 16. Where everything lives

| Path | What |
| --- | --- |
| `reports/phase4/MANUSCRIPT.md` | The formal write-up |
| `reports/phase4/PRE_DECLARATION.md` | Analysis plan and compliance audit, including every deviation |
| **`data/release/reclassification_table.tsv`** | **The deposited resource — the paper's main deliverable** |
| `data/release/accession_list.tsv` | All 1,143 accessions with their fate |
| `data/release/master_variant_table.tsv` | Every variant with frequency, flags and features |
| `data/release/DATA_DICTIONARY.md` | Column definitions |
| `data/processed/feature_correlation_matrix.tsv` | Published in full |
| `reports/phase4/figures/` | The three manuscript figures, plus every plotted value as JSON |
| `reports/phase0/` … `phase4/` | Per-phase reports, including bugs found and fixed |
| `scripts/` | The 21 pipeline stages and the verifier |

---

## 17. Glossary

| Term | Meaning |
| --- | --- |
| **Allele frequency** | How often a particular residue occurs at a position across a population |
| **Fold change** | How much more drug is needed to inhibit a mutant than the wild type. 2× is this project's declared resistance threshold |
| **pN/pS** | Ratio of amino-acid-changing to silent mutations. Below 1 means purifying selection — the gene is under pressure to stay the same |
| **Average precision (AP)** | A ranking metric. Only meaningful against its no-skill floor, which equals the positive-class prevalence |
| **Lift** | AP minus the floor. The number that actually means something |
| **Prevalence / floor** | The fraction of positives. A model that ranks randomly achieves exactly this AP |
| **MCC** | Matthews correlation coefficient — a balanced measure of decisions at a threshold, as opposed to ranking |
| **Tier A/B/C** | Evidence quality: engineered (strongest), lab-selected, clinical isolate |
| **Study-dominated** | More than half the carriers come from a single submission — a frequency you should not trust alone |
| **LOSO** | Leave-one-study-out: recompute the frequency with the most influential study removed |
| **Permutation null** | Shuffle the labels and re-score. Tells you what "no signal" actually looks like under your scheme — which is not always zero |
| **Residualisation** | Remove the part of feature A that feature B explains, then ask whether what remains still predicts |
| **Holm correction** | Adjusts p-values for testing many models at once |
| **VUS** | Variant of unknown significance — seen in a patient, no phenotype available |
