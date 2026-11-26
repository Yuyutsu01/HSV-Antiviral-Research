# Most catalogued HSV antiviral resistance mutations have never been checked against the virus that circulates

*Authors and affiliations to be completed.*

**Correspondence:** *to be completed*

**Keywords:** herpes simplex virus; antiviral resistance; allele frequency;
thymidine kinase; DNA polymerase; helicase–primase; aciclovir; pritelivir

**Running title:** Population frequencies of catalogued HSV resistance mutations

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
599 publicly deposited genomes, and used them to recontextualise the catalogue.

Of 644 catalogued substitutions, **128 occur in circulating virus and 516 do
not**. Eleven are the *majority* allele, meaning the reference sequence carries
the minority residue; a further 15 exceed 10%. Seven are recommended for
reclassification as polymorphisms. Forty-one are common in circulation but have
never been phenotyped, and are the clearest target for future assay work.
Three initially appeared to be circulating resistance alleles; none survives
scrutiny: one is carried by a single submission, one is dominated by one
submission among three, and one rests on a 2.3-fold phenotype.
Among the 241 substitutions with a measured phenotype, none that exceeds 5%
frequency exceeds 2.4-fold resistance: there is no common strongly resistant
allele in circulating HSV.

Thymidine kinase is under markedly weaker purifying selection than DNA
polymerase (pN/pS 0.371 against 0.173; difference +0.1980, 95% CI +0.1664 to
+0.2361, 100% of bootstrap resamples), consistent with its dispensability for
replication in vitro and with the greater tolerance of its sequence for
substitution.

No isolate among 251 HSV-1 UL5 and 299 HSV-1 UL52 sequences carries a
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
0.411). It survives a prospective test: across five temporal hold-outs
splitting the literature from 2005 to 2015, structure is the best feature group
at every cut-off, at a lift of +0.170 to +0.229 against a best-of-the-rest that
never exceeds +0.163. At the pre-declared 2-fold operating point the best
model's Matthews correlation is 0.15, and 8 of 22 models classify every variant
identically. A model that ranks well and decides uselessly is the honest
summary.

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

HerpesDRG (Charles et al. 2024), pinned at commit `22fc186` (retrieved 2026-11-01T19:30:05Z, MD5
`d78ad69a…`), active entries only (`status == "A"`). GenBank nucleotide records
for taxa 10298 (HSV-1) and 10310 (HSV-2) restricted to 140–200 kb, queried
2026-11-01T20:02:36Z and 20:07:23Z. All 1,143 returned accessions are listed in
`data/release/accession_list.tsv` **with their fate**, so that an absent
accession can be distinguished from an excluded one.

### 2.2 From 1,143 accessions to 599 isolates

| Fate | n |
| --- | --- |
| Analysed | **599** |
| Duplicate of an isolate already counted | 187 |
| No target-gene CDS extracted | 168 |
| Failed CDS quality control | 127 |
| Laboratory-adapted strain | 30 |
| Serial isolate from an already-represented patient | 18 |
| Non-human host | 14 |

De-duplication operates at three levels: identical sequence across all four
genes, strain name, and patient. The three remove 281, 5 and 18 records; the
first is much the largest and is the reason the analysed set is 599 rather
than the 943 records from which a target gene could be read at all. The
patient level matters more than it sounds — see §3.1.

Per-gene coverage: HSV-1 239–300 sequences, HSV-2 160–253.

### 2.3 Evidence tiers

The catalogue holds 700 active entries for the four genes, of which 644 are
single-residue substitutions; the remaining 56 are 30 nonsense, 24
indel/frameshift and 2 other variants, which this work does not analyse. Every
count in this paper is over the 644 unless stated otherwise.

Tier is taken from `strain_generation`, keeping the strongest available per
mutation. Across the 644 substitutions: **A, engineered** (recombinant BAC or
marker transfer) **133**, **B, lab-selected** **66**, **C, clinical isolate**
**445**. Tier A is the primary set for anything involving phenotype, because
only there is the genotype–phenotype link experimentally isolated.

**Of the 133 Tier A substitutions, 119 reach the modelling set**: 11 carry a
qualitative resistance call with no numeric fold change, and 3 more are
dropped as co-mutations or for want of a reference protein. The 119 is the
primary set throughout §§5 and 6.

**How the fold change is formed, and the limitation it carries.** A catalogue
entry may hold several phenotype records, from different studies, assays and
drugs. The value used here is the maximum over all of them:

$$\mathrm{FC}(v) = \max_{r \,\in\, \mathcal{R}(v)} \; \max_{d \,\in\, \mathcal{D}(r)} \mathrm{fc}(r, d),$$

for records $\mathcal{R}(v)$ of variant $v$ and drugs $\mathcal{D}(r)$ tested in
record $r$. Taking the maximum is the conservative choice for Aim 1, whose
question is whether a *common* allele could be resistant *under any reading of
the catalogue* — understating a phenotype would manufacture the reassuring
answer. It has a cost, declared here because it matters in §3.1: the maximum
may come from a drug whose activation does not involve the mutated gene, in
which case it is a measurement of assay noise rather than of resistance.

### 2.4 Alignment, and a constraint worth declaring

No multiple-sequence-alignment tool was available: neither MAFFT (Katoh and
Standley 2013) nor BLAST (Altschul et al. 1990) could be installed in the
environment this was developed in. Alignment is therefore reference-anchored
and pairwise (Biopython `PairwiseAligner`, Cock et al. 2009), scored with
BLOSUM62 (Henikoff and Henikoff 1992) and free terminal gaps. For coding sequences of 70–80% identity to a known reference
this is adequate, and it makes every position directly comparable to reference
coordinates without a column-mapping step. It is nonetheless a limitation, and
indels are handled less well than a codon-aware MSA would handle them.

### 2.5 Frequencies and ascertainment

A variant's frequency is its carrier count over the number of isolates with
usable coverage at that position, $\hat{p} = k/n$, reported with a 95% Wilson
score interval (Wilson 1927). For $z = 1.96$ the interval is

$$
\frac{\hat{p} + \dfrac{z^{2}}{2n}}{1 + \dfrac{z^{2}}{n}}
\;\pm\;
\frac{z}{1 + \dfrac{z^{2}}{n}}
\sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^{2}}{4n^{2}}},
$$

truncated to $[0, 1]$. The Wilson interval rather than the Wald interval
$\hat{p} \pm z\sqrt{\hat{p}(1-\hat{p})/n}$ because nearly every frequency here
is small and many are zero, where Wald coverage falls far below nominal and the
interval collapses to zero width at $\hat{p} = 0$ (Agresti and Coull 1998).
Implementation: `scripts/p1_04_frequencies.py`.

Every frequency additionally carries:

- `n_studies_carrying` — distinct submissions contributing a carrier;
- `loso_min_freq` — the lowest frequency under leave-one-study-out;
- `max_study_share` — the largest share of carriers from one submission;
- `study_dominated` — whether dropping one study halves the frequency.

These are not supplementary. A 5% allele drawn entirely from one submission and
a 5% allele seen across five are different claims, and the table refuses to
state one without the other.

### 2.6 Selection

pN/pS by the Nei–Gojobori counting method (Nei and Gojobori 1986) with a
Jukes–Cantor correction (Jukes and Cantor 1969), measured **against the
population consensus rather than the reference**. Using the reference would
count every position where the reference is the minority residue as
polymorphism, inflating the estimate by a term that has nothing to do with
selection.

For a codon $c$, the synonymous and non-synonymous site counts are

$$
S(c) = \sum_{i=1}^{3} \frac{1}{3} \sum_{b \neq c_i}
  \mathbf{1}\!\left[\, \tau(c^{i \to b}) = \tau(c) \,\right],
\qquad N(c) = 3 - S(c),
$$

where $c^{i \to b}$ is $c$ with position $i$ replaced by base $b$ and $\tau$ is
the standard genetic code. Each of the three alternative bases at a position
contributes $1/3$ of a site, so the three positions sum to 3.

For a pair of aligned codons differing at $m$ positions, the synonymous and
non-synonymous **differences** are averaged over all $m!$ shortest mutational
pathways between them, each pathway weighted equally:

$$
S_d = \frac{1}{m!} \sum_{\pi \in \Pi_m} s(\pi),
\qquad
N_d = \frac{1}{m!} \sum_{\pi \in \Pi_m} \bigl(m - s(\pi)\bigr),
$$

with $s(\pi)$ the number of synonymous steps along pathway $\pi$. This is what
makes the method a *counting* method rather than a model fit, and it is the
step most often omitted when "Nei–Gojobori" is cited; it is written out here
because the choice is not neutral for codons differing at two or three
positions. Implementation: `scripts/p2_02_selection.py`.

Summing $S$, $N$, $S_d$ and $N_d$ over all compared codons gives the raw
proportions $p_S = S_d / S$ and $p_N = N_d / N$, to which the Jukes–Cantor
correction for multiple hits at the same site is applied:

$$
d = -\frac{3}{4}\ln\!\left(1 - \frac{4p}{3}\right),
\qquad \text{undefined for } p \geq \tfrac{3}{4},
$$

and the reported statistic is the ratio $d_N / d_S$. At the divergences seen
here ($p \ll 0.1$) the correction is small, and it is applied rather than
skipped so that genes of different divergence remain comparable. Confidence
intervals are percentile bootstrap over isolates, 1,000 replicates.

A gene is reported as reliable only when synonymous sites are
sufficient; HSV-2 UL23 (Sd = 32) is not, and the HSV-2 contrast is withheld
rather than reported with a wide interval.

### 2.7 Features and models for Aim 4

Six feature groups: frequency, conservation (per-site non-synonymous rate),
amino-acid properties, BLOSUM62 exchangeability, structure, and ESM-2. Each is
run alone as a baseline and the staged models add them in that order.

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

Structures are taken from the Protein Data Bank (Berman et al. 2000). Solvent
accessibility is computed by the Shrake–Rupley algorithm (Shrake and Rupley
1973) as implemented in Biopython, and reported as a **relative** accessibility

$$\mathrm{RSA}(i) = \frac{\mathrm{SASA}(i)}{\mathrm{ASA}_{\max}(a_i)},$$

where $\mathrm{ASA}_{\max}$ is the theoretical maximum accessible area for the
residue's amino-acid type, taken from Tien et al. (2013). Normalising matters
because raw SASA is not comparable across residue types: a fully exposed
glycine and a buried tryptophan can return similar absolute areas.

**UL52 has no resolved site of its own** — zero residues within 5 Å of any
ligand, median distance 42 Å, its nearest ligand belonging to a neighbouring
subunit. Those rows carry `site_resolved = 0`. This turns out to matter a great
deal (§5.3).

ESM-2 `esm2_t33_650M_UR50D` (Lin et al. 2023), scored zero-shot by the masked
marginal of Meier et al. (2021). The position of interest is replaced by the
mask token, the model is run once over the masked sequence, and the variant is
scored by the log-ratio of the two residues' probabilities at that position:

$$
s(i, a \!\rightarrow\! b) =
  \log P\!\left(x_i = b \mid x_{\setminus i}\right)
  - \log P\!\left(x_i = a \mid x_{\setminus i}\right),
$$

both terms taken from the same softmax over the amino-acid vocabulary, computed
in float32 regardless of the model's inference precision. A negative score
means the model prefers the reference residue. Nothing is fitted: no label is
seen, so this feature cannot leak the target.

Deep networks and fine-tuned language models were excluded in advance: at a few
hundred labelled positives they cannot be fitted or validated.

**Cross-validation, stated completely.** The unit is one amino-acid
substitution; the primary set is the 119 Tier A variants. Nothing is split at
random, because positions within a gene are not independent — they share a
structural context, a conservation column and often a source study, and a
random fold would place neighbouring residues on both sides of the split
(Roberts et al. 2017). Two grouped schemes are run:

| Scheme | Group key | Groups | Outer splitter | Test rows per fold |
| --- | --- | --- | --- | --- |
| Leave-one-gene-out | `virus:gene` | 6 | `LeaveOneGroupOut` | 6 to 45 |
| Position blocks | `virus:gene:floor(pos/50)` | 44 | `GroupKFold(5)` | about 24 |

The six gene groups hold 45 (HSV-1 UL30), 31 (HSV-1 UL23), 18 (HSV-2 UL30), 12
(HSV-2 UL23), 7 (HSV-1 UL5) and 6 (HSV-1 UL52) variants. The 44 position
blocks hold 1 to 13 variants each and are assigned whole to folds, so no
residue within 50 positions of a test residue in the same gene is ever in
training.

Hyperparameters are selected **inside each training fold only**, by an inner
`GroupKFold` with $\min(3, |G_{\text{train}}|)$ splits over that fold's own
groups, scored by Spearman $\rho$. The selected configuration is refit on the
whole training fold and predicts the test fold once; this nesting is what keeps
the reported figure from being a model-selection score rather than a
generalisation score (Varma and Simon 2006). Estimators are ridge regression
(Hoerl and Kennard 1970) over $\alpha \in \{0.1, 1, 10, 100\}$ and a random
forest (Breiman 2001) of 300 trees over `min_samples_leaf` $\in \{1, 5\}$ and
`max_features` $\in \{0.5, 1.0\}$. Every fold's predictions are written back
into a single out-of-fold vector of length 119, and that vector — not any
within-fold score — is what §5 evaluates.

The primary metric is average precision, reported beside the positive-class
prevalence; the 2-fold threshold and the metric hierarchy were fixed before
data acquisition (`PRE_DECLARATION.md`).

### 2.8 Statistical definitions

Every estimator and test used in §§5 and 6 is defined here, so that no reported
number depends on a convention the reader has to guess. Implementation:
`scripts/p4_01_metrics.py` and `scripts/p3_05_audit.py`.

**Average precision.** Variants are ranked by score, descending, with stable
tie-breaking. With $y_k \in \{0,1\}$ the label of the $k$-th ranked variant and
$P = \sum_k y_k$ the number of positives,

$$
\mathrm{AP} = \frac{1}{P} \sum_{k=1}^{n} y_k \cdot
  \frac{\sum_{j \le k} y_j}{k}.
$$

This is the step-wise sum — precision at each positive, averaged over the
positives — not the trapezoidal area under an interpolated precision–recall
curve. The two differ at small $n$, and the choice is stated because it moves
the third decimal of every figure in §5.2. A scoring independent of the labels
has $\mathbb{E}[\mathrm{AP}] = P/n$, the prevalence, which is why AP is never
reported here without it (Davis and Goadrich 2006; Saito and Rehmsmeier 2015).

**Precision at 95% recall.** Precision at the smallest rank $k$ at which
cumulative recall reaches 0.95.

**Matthews correlation** (Matthews 1975), at the pre-declared operating point
$s \ge \log_{10} 2$:

$$
\mathrm{MCC} = \frac{TP \cdot TN - FP \cdot FN}
  {\sqrt{(TP{+}FP)(TP{+}FN)(TN{+}FP)(TN{+}FN)}}.
$$

When the prediction is constant at that threshold the denominator is zero and
MCC is undefined. That is reported as `threshold_degenerate`, not as a missing
value: it is the finding that the model makes no usable decision.

**Rank agreement** is Spearman's $\rho$ (Spearman 1904) between predicted score
and observed $\log_{10}$ fold change — the statistic in which every figure in
§6.1, §6.2 and §6.3 is quoted.

**Calibration** is the ordinary least-squares slope and intercept of observed
$\log_{10}$ fold change on predicted score. Slope 1 with intercept 0 is
perfect; slope below 1 is under-dispersion, that is, shrinkage toward the mean.

**Bootstrap intervals** are percentile intervals resampled **over mutations**
(Efron 1979): 2,000 replicates for every Phase 4 metric, 1,000 for the Phase 2
selection ratios, seeded once per run. Resampling mutations rather than
predictions is deliberate — the latter treats each prediction as an independent
draw and understates the interval, when the sampling unit is the variant.

**Model-versus-baseline comparison** is a paired bootstrap of the AP advantage.
On each of 2,000 resamples $b$ the same resampled mutations score both models,
giving $\Delta_b = \mathrm{AP}_b(\text{model}) - \mathrm{AP}_b(\text{frequency})$;
resamples containing only one class are discarded. Pairing on the same resample
removes between-resample variance from the comparison. The interval is the
2.5th and 97.5th percentiles of $\{\Delta_b\}$, and the two-sided $p$ is the
sign-flip proportion, doubled and capped:

$$
p = \min\!\left(1,\; 2 \min\left(
  \widehat{\Pr}[\Delta \le 0],\;
  \widehat{\Pr}[\Delta \ge 0] \right)\right).
$$

**Multiplicity.** Holm–Bonferroni (Holm 1979) across the 21
model-versus-baseline comparisons. With $p_{(1)} \le \dots \le p_{(m)}$ ordered,
the adjusted value is

$$
\tilde{p}_{(i)} = \max_{j \le i} \; \min\!\left(1,\; (m - j + 1)\, p_{(j)}\right),
$$

the running maximum enforcing monotonicity. Holm rather than
Benjamini–Hochberg because the question is whether *any* model beats the
baseline, which calls for control of the family-wise error rate rather than of
the false discovery rate.

**Permutation nulls.** The target is permuted **within** cross-validation
groups, preserving each group's mean label and destroying only the
feature–label link, and the whole nested cross-validation is re-run. With
$R = 100$ replicates and $\rho_{\text{obs}}$ the observed value, the one-sided
$p$ uses the add-one estimator

$$
p = \frac{1 + \#\{\, \rho^{*}_{r} \ge \rho_{\text{obs}} \,\}}{1 + R},
$$

which cannot return zero: the smallest attainable value is
$1/101 \approx 0.0099$, and no $p$ below that is claimed anywhere in this
paper. A result is called usable only when $p < 0.05$ **and**
$\rho_{\text{obs}} > 0$, because §6.1 shows the null is not centred at zero.

**Conservation-matched negatives.** Conservation is binned at the quartiles of
the pooled distribution, so that neither class defines the grid. Within each
bin both classes are subsampled to the smaller of the two counts, equalising
the conservation distribution by construction; 100 independent draws are taken
and the reported figure is the mean with a percentile interval. Both classes
are subsampled rather than only the negatives because the negatives are the
minority here (32 against 87), so a draw sized to the positives would silently
return the whole pool and test nothing.

**Residualisation.** For each ESM feature $E_j$, an OLS fit on the conservation
and structure features $C$ gives $\hat{E}_j = C\hat{\beta}_j$, and the feature
is replaced by its residual $E_j - \hat{E}_j$. What survives is the part of
ESM-2's signal that conservation and structure cannot account for.

**Seeds.** One seed per stage, recorded in that stage's JSON summary beside the
replicate counts, so every interval and $p$-value in this paper is reproducible
exactly.

## 3. Aim 1 — what the catalogue looks like with a denominator

Of 644 catalogued substitutions, 128 occur in circulating virus (Figure 1A).

| Band | n | |
| --- | --- | --- |
| Not observed | 516 | No population evidence either way |
| Singleton | 45 | One carrier |
| Rare (<1%) | 21 | |
| Uncommon (1–10%) | 36 | |
| **Common (10–50%)** | **15** | |
| **Majority (>50%)** | **11** | The reference carries the minority residue |

Interpretive flags:

| Flag | n | What it means |
| --- | --- | --- |
| `reference_artifact` | 10 | Majority allele, no measured resistance |
| `reclassify_as_polymorphism` | 7 | Catalogued as resistance-associated, above 1%, fold change at or below 2× |
| `unphenotyped_but_common` | 41 | Above 1% in circulation, never phenotyped |
| `consistent_rare` | 60 | Rare, as a resistance determinant should be |
| `consistent_with_resistance` | 6 | Rare and phenotypically resistant |
| `circulating_resistance_allele` | 3 | Between 1% and 50%, fold change above 2× — all three fail on inspection |
| `majority_allele_with_measured_resistance` | 1 | Above 50% and above 2×: HSV-2 UL23 G39E, at 79% and 2.4× |

**The seven `reclassify_as_polymorphism` entries are the paper's concrete
deliverable**, and the 41 `unphenotyped_but_common` entries are its clearest
call to action: each is carried by more than 1% of circulating virus, so a
laboratory will meet them and has no phenotype to interpret them with.

### 3.1 The three apparent circulating resistance alleles

| Variant | Freq | Carriers | Studies | `max_study_share` | `loso_min` | FC | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| UL23 R222C | 0.046 | 11 | 3 | 0.73 | 0.0134 | 98× | Ascertainment artefact |
| UL23 R176Q | 0.017 | 4 | 1 | 1.00 | 0.0 | 200× | Ascertainment artefact |
| UL23 R41H | 0.051 | 12 | 5 | 0.33 | 0.0381 | 2.3× | Sound frequency, phenotype does not hold up |

R222C is the entry that would have mattered most: 98-fold resistance at
an apparent 7.6% frequency would have been a surveillance finding. It is not.
The 98-fold itself is the aciclovir value from a laboratory-selected strain
(Schalkwijk et al. 2024); the catalogue's clinical-isolate records for the same
substitution give 7.1-fold and 9.8-fold (Gaudreau et al. 1998) and 25-fold
(Bestman-Smith et al. 2001), with a further record reporting a qualitative
resistant call (Schubert et al. 2014). The substitution is genuinely
resistance-associated. What fails is the frequency.
Patient-level collapse reduced it from 19 carriers to 11 and 7.6% to 4.6%, and 8
of the 11 remaining carriers come from one submission; dropping that submission
takes the frequency to 1.4%.

The obvious challenge to that verdict is that de-duplication discarded real
carriers, so every record in the frozen query carrying R222C was re-examined:
38 in total, of which 11 are counted. Of the other 27, ten are serial isolates
from four patients already represented — five of them from one patient, with
strain names of the form `10-6-1`, `10-7-1`, `10-11-3` — and seventeen are
sequence-identical across all four target genes to a record already kept.
**None is an independent carrier.** Counting all 38 would still leave 24 of
them, 63%, from the single dominant submission.

R41H is the interesting case, because it is **not** an artefact — 12 carriers
across five submissions, robust to study dropout. What disqualifies it is the
phenotype, and the catalogue's own records say so once they are read drug by
drug. Three independent studies contribute. Chibo et al. (2004) record it as a
polymorphism for both aciclovir and foscarnet. Brunnemann et al. (2016), in a
recombinant BAC background, give aciclovir 1.0-fold, penciclovir 1.1-fold and
brivudine 1.0-fold. Inagaki et al. (2018), also in a recombinant background,
give aciclovir 0.8-fold.

**The 2.3-fold that put R41H in this table is the foscarnet value from
Brunnemann et al. (2016).** Foscarnet is a pyrophosphate analogue that inhibits
the viral DNA polymerase directly and requires no activation by thymidine
kinase, so a UL23 substitution has no mechanism by which to shift it; a
2.3-fold reading on an assay whose other three drugs sit at 1.0-fold is better
read as the assay's own dispersion. For aciclovir — the drug UL23 does
activate — the three studies agree on 0.8 to 1.0-fold.

This is the cost of the maximum-over-drugs rule declared in §2.3, and R41H is
the case where it bites. The rule is kept because Aim 1's question is whether
any common allele could be resistant on any reading of the catalogue, and
loosening it would manufacture the reassuring answer. The verdict stands, with
a firmer basis than "marginal": R41H is a common natural polymorphism of HSV-1
UL23 with no demonstrated resistance phenotype against the drug its gene acts
on. Frequency and phenotype each need their own scrutiny, and
`study_dominated = 0` alone would have retained it.

**60 of the 128 observed variants are study-dominated.** Any use of these
frequencies that ignores that column will overstate about half of them. The 11
reference artefacts, by contrast, are all non-dominated at study shares around
0.15 — they are real.

### 3.2 Frequency and measured resistance are close to disjoint

Plotting the two together bounds the surveillance question directly (Figure 1B).
Of the 241 catalogued substitutions carrying a numeric phenotype, 30 are
observed in circulating virus and 211 are not, and the two groups differ in
measured resistance: median fold change 1.5× for the observed and 5.0× for the
unobserved.

The informative part of the plot is empty. Taking 5% as the cut — half the
band boundary used in the table above, so that the claim is made over a wider
set of substitutions rather than a narrower one — the strongest phenotype
above it is **2.4×**. Two substitutions clear the declared 2-fold threshold at
that frequency, both marginally (HSV-1 UL23 R41H at 2.3× and HSV-2 UL23 G39E
at 2.4×, the latter a majority allele at 79%), and **none reaches 5×**. Both
dissolve on inspection in the same way: §3.1 shows R41H's 2.3× to be a
foscarnet reading on a thymidine-kinase substitution, and G39E's 2.4× is a
penciclovir dye-uptake value from Brunnemann et al. (2018), whose own
plaque-reduction assay on the same recombinant gives 0.5× for penciclovir and
0.7× for aciclovir, while Burrel et al. (2012) report 1.2× by functional
kinetics. Neither of the two substitutions that clear the threshold does so on
the drug its gene activates.

The nearest approach to a counterexample is UL23 R222C, at 98-fold and 4.6%,
just below the cut — and §3.1 is the reason it does not qualify on frequency
either: 8 of its 11 carriers come from one submission, and dropping that
submission takes it to 1.4%. Subject to that, there is no substitution in
circulating HSV that is both reasonably common and strongly resistant by
measurement. The absence is the reassuring reading of Aim 1, and it is only
visible once the denominator exists.

## 4. Aims 2 and 3

### 4.1 Thymidine kinase is under weaker constraint than polymerase

| Gene | pN/pS | 95% CI | Sd |
| --- | --- | --- | --- |
| HSV-1 UL23 | **0.371** | 0.338–0.407 | 620 |
| HSV-1 UL30 | **0.173** | 0.161–0.185 | 1536 |
| HSV-1 UL5 | 0.181 | 0.159–0.205 | 955 |
| HSV-1 UL52 | 0.140 | 0.128–0.154 | 2014 |

Difference UL23 − UL30: **+0.1980 [+0.1664, +0.2361]**, with UL23 higher in
**100%** of 1,000 bootstrap resamples. The intervals do not approach each other
(Figure 2).

This is the expected direction — thymidine kinase is dispensable for replication
in culture while polymerase is not — and it supplies a population-genetic
reason why TK accumulates substitutions that then appear in resistance
catalogues. The contrast is also insensitive to the two choices of sequence set
that move individual allele frequencies most. Collapsing serial isolates to one
per patient, which moves R222C's frequency by 40%, moves this contrast by
0.004. Restricting extraction to the 451 isolates whose GenBank records carry
an explicit `gene` qualifier, discarding a quarter of the sequence and the
HSV-2 records disproportionately, moves it by 0.0025. Both lie well inside the
interval of ±0.03. Gene-wide rates are insensitive to sampling in a way
single-variant frequencies are not, which is the general reason to prefer them
for a selection contrast.

HSV-2 is withheld: UL23 has only 32 synonymous sites and a pN/pS point estimate
of 2.46 with an interval spanning 1.74–3.67, which is an artefact of the
denominator rather than evidence of diversifying selection.

### 4.2 The helicase–primase baseline is effectively zero

| Gene | Isolates | Carrying a catalogued substitution |
| --- | --- | --- |
| HSV-1 UL5 | 251 | **0** |
| HSV-1 UL52 | 299 | **0** |
| HSV-2 UL5 | 247 | 2 (0.8%) |
| HSV-2 UL52 | 225 | **0** |

Pritelivir and amenamevir are new enough that this is a genuine pre-therapeutic
baseline, and it is the cleanest result in the project: there is essentially no
pre-existing helicase–primase resistance in circulating HSV. Any future rise is
measurable against this.

One caveat against over-reading it: the catalogue lists only 8 HSV-1 UL5 and 6
UL52 positions, so the test has little to find. A permutation test of whether
catalogued positions are more conserved than random ones is non-significant for
HSV-1 UL5 (p = 0.54) and UL52 (p = 0.37). For HSV-2 UL5 catalogued positions
come out 2.2× *more* variable than gene-wide — the opposite of constraint — but
at p = 0.052 that is not significant either, and it rests on four positions. A
test with four informative positions moves across a conventional threshold on
the addition or removal of a few dozen sequences, so the ratio is reported and
not interpreted.

## 5. Aim 4 — models, and what beating a frequency baseline is worth

### 5.1 Average precision has to be read against its floor

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
| BLOSUM62 only (forest) | 0.802 | +0.071 | 0.722 | −0.10 |
| Amino-acid properties only | 0.751 | +0.020 | 0.735 | 0.02 |
| **Frequency only** | **0.723** | **−0.008** | 0.762 | −0.08 |
| Conservation only | 0.717 | −0.014 | 0.741 | 0.07 |
| BLOSUM62 only (ridge) | 0.715 | −0.016 | 0.735 | degenerate |

**Frequency scores below its own no-skill floor, and so does conservation.**
"The models beat the frequency baseline" is true and nearly vacuous: the
baseline has no skill to beat. The reason is structural — only **14 of 119**
phenotyped Tier A variants are observed in circulating virus at all, so the
frequency features are near-constant on this set.

Under multiplicity correction across the 21 comparisons — every model except
the frequency baseline itself, which is the thing being compared against —
**only 3 models beat frequency on average precision** (Holm-adjusted p < 0.05): the two
structure-only variants, and the full model with ESM-2 added. Seven clear
p < 0.05 uncorrected, so the correction is what separates a reportable result
from a grid of near-misses. The protocol anticipated this in terms: "an uncorrected
pairwise win on one split is not evidence of superiority."

### 5.3 Where the structural signal actually lives

Figure 3A gives the same breakdown against each gene's prevalence floor.

| Gene | n | Prevalence | Structure AP | Frequency AP |
| --- | --- | --- | --- | --- |
| HSV-1 UL23 | 31 | 0.84 | **0.988** | 0.774 |
| HSV-2 UL23 | 12 | 0.67 | **1.000** | 0.647 |
| HSV-1 UL30 | 45 | 0.62 | 0.823 | 0.622 |
| HSV-1 UL5 | 7 | 0.57 | **0.567** | 0.525 |
| HSV-1 UL52 | 6 | 0.50 | **0.411** | 0.533 |
| HSV-2 UL30 | 18 | 1.00 | undefined | undefined |

The six rows sum to the 119 primary variants. HSV-2 UL30 contributes 18, all
of them resistant, so its prevalence is 1.00 and average precision is
undefined on it — a ranking metric needs both classes present. It is shown
rather than dropped because omitting it would leave the breakdown summing to
101 against a 119 used everywhere else.

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
is **0.15**, and **8 of 22 models predict a single class for every variant** —
their MCC is undefined, not missing. Calibration is not the problem: the three
models that survive correction have slopes of 0.77, 0.95 and 1.20, close enough
to one. The problem is that the feature sets that fail do so by predicting
almost nothing at all — frequency's calibration slope is 0.06 and
conservation's 0.17, which is a near-constant prediction dressed as a model.

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
cross-validation. The full model does not — it collapses to AP 0.632, below the
floor, which is what overfitting 30 training rows across five feature groups
looks like.

### 5.6 External validation: the temporal hold-out

The protocol's strongest declared test: train on mutations characterised before
a cut-off year, test on those characterised after. It approximates prospective
validation — knowing only what was published by a given date, can the model
identify what was subsequently discovered?

A mutation's year is the earliest among its references. All 89 distinct study
DOIs resolved to a publication year through Crossref, spanning 1987–2025 across
10 publishers, so **every one of the 235 labelled variants carries a year**
and the split is not confounded with journal. The secondary set of 235 is all
tiers together and contains the 119 Tier A primary set, rather than standing
beside it.

The cut-off was not pre-declared. Rather than choose one, all five are reported;
2012 is nominated primary because it splits the primary set most evenly, a
criterion independent of the outcome.

Average-precision **lift over the test-set prevalence floor**, primary Tier A
(Figure 3B):

| Cut-off | Train → test | Structure | Frequency | Conservation | ESM-2 | Full |
| --- | --- | --- | --- | --- | --- | --- |
| 2005 | 41 → 78 | **+0.195** | −0.077 | +0.085 | +0.118 | +0.110 |
| 2008 | 47 → 72 | **+0.215** | +0.019 | +0.100 | +0.052 | +0.126 |
| 2010 | 52 → 67 | **+0.218** | +0.020 | +0.062 | −0.011 | +0.123 |
| **2012** | **59 → 60** | **+0.229** | +0.030 | +0.042 | +0.050 | +0.090 |
| 2015 | 70 → 49 | **+0.170** | −0.006 | +0.162 | +0.101 | +0.066 |

Secondary set, 235 variants: structure lift **+0.245 to +0.349**, frequency
**+0.001 to +0.025**.

**Structure-only is the best model at every cut-off on both sets**, with a lift
that never falls below +0.13 and no alternative reaching +0.163 at any split.
Nothing else is stable: conservation ranges from +0.042 to +0.162 across the
five splits, ESM-2 from −0.011 to +0.118, and the full model from +0.066 to
+0.126 — all positive, none within reach of structure, and none holding a
consistent rank against each other. Frequency straddles zero throughout.

This is the most persuasive single result in the project. Five independent
prospective splits, no tuning against any of them, and the same feature group
wins each time by a margin no alternative approaches.

**One confound the split exposes, which is a finding about the literature
rather than about the models.** Resistance prevalence falls sharply from train
to test at every cut-off: 0.86 → 0.60 on the primary set at 2012, and 0.83 →
0.52 on the secondary. Early resistance work characterised almost exclusively
resistant variants; later work reports many more negatives. So the test sets are
systematically harder than the training sets, which makes a stable positive lift
more impressive than it would otherwise be — and it means a model trained on
pre-2012 literature is calibrated to a base rate that no longer holds.

## 6. The circularity audit

Reported here rather than in a supplement, because it is the project's main
methodological contribution.

### 6.1 Leave-one-gene-out has a negative null

Permuting the target **within** cross-validation groups — preserving each gene's
mean label, destroying only the feature–label link (§2.8) — gives a null
centred at a Spearman $\rho$ of **−0.16 to −0.29**, not zero. Every figure in
this section is a Spearman rank correlation between predicted score and
observed $\log_{10}$ fold change.

The mechanism: the model predicts approximately the training mean for a
held-out gene, and removing a gene shifts that mean *away* from the gene's own
mean. Where genes differ in average resistance, and they do, the scheme
manufactures anti-correlation from noise.

Scored against zero, frequency-only (−0.150), conservation-only (−0.228) and
BLOSUM-only (−0.254) look "worse than chance". They are on the noise floor.

Position-blocks is the better-behaved scheme and carries the headline results,
but only for the feature groups that have signal: its null sits between −0.014
and +0.079 for structure, ESM-2 and the full model. For frequency (−0.167) and
conservation (−0.115) the null is substantially negative under *both* schemes,
so no scheme here licenses scoring those two against zero. **Any grouped-CV
figure quoted against a zero null inverts the conclusion for the weaker feature
sets**, and the weaker feature sets are exactly the ones a reader would want to
compare against.

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
| Position blocks | Full | 0.693 | 0.706 | −0.013 |
| Position blocks | Full (forest) | 0.739 | 0.754 | −0.015 |
| Leave-one-gene-out | Structure | 0.733 | 0.732 | +0.001 |
| Leave-one-gene-out | Structure (forest) | 0.754 | 0.755 | −0.001 |
| Leave-one-gene-out | Full | 0.621 | 0.624 | −0.003 |
| Leave-one-gene-out | Full (forest) | 0.693 | 0.693 | +0.001 |

**The gap is zero across all eight combinations**, range −0.016 to +0.003. The
structural signal is not a conservation proxy.

The size matching is what makes that statement possible. Conservation matching
necessarily discards rows — 64 of the 119 survive it here — so a matched arm
compared against the full set differs in sample size as well as in conservation,
and average precision is sensitive enough to n at this scale that such a
comparison returns a gap of ambiguous origin. Pairing each matched draw with a
random draw of identical size and class balance removes that confound, and only
then does a gap of zero mean that conservation was not carrying the signal.

### 6.3 ESM-2 does not survive residualisation

| Scheme | ESM-2 alone | Residualised against conservation + structure |
| --- | --- | --- |
| Position blocks | ρ +0.305 | **ρ −0.011** |
| Leave-one-gene-out | ρ +0.138 | **ρ −0.113** |

Conservation and structure linearly explain only 34% of the ESM-2 masked-marginal
score, so on a variance accounting two thirds of ESM-2 is something else. **That
two thirds carries no predictive signal.** Removing the explainable 34% does not
merely weaken ESM-2's ranking of fold change, it eliminates it: the residual
correlates with the target at −0.011 and −0.113, either side of zero, under the
two schemes.

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
| **`cons_site_variability`** | **`freq_frequency`** | **+0.888** |
| `cons_site_variability` | `freq_count` | +0.888 |
| `cons_site_variability` | `freq_observed` | +0.886 |
| `esm_masked_marginal` | `esm_ref_logprob` | −0.889 |

Two things follow. **Conservation and frequency are very nearly the same
feature** (ρ = +0.888) — which is why they fail together and in the same
direction, and why "add conservation to frequency" never helped. And **no
structural feature appears in this list at all**: structure is uncorrelated
with conservation and frequency, which is precisely why it is the only group
that works and why §6.2 comes out at zero.

The three frequency features being mutually collinear at ≈1.0 also confirms
Aim 1's finding that the frequency signal is **one bit**: across all 241
variants with a numeric phenotype, continuous frequency gives ρ = −0.248 and
the binary observed/not-observed flag gives ρ = −0.243, both p = 0.0001.
Knowing how common an allele is adds nothing over knowing whether it has ever
been seen. Median fold change is 1.5× for observed variants and 5.0× for
unobserved.

### 6.5 Leakage

No evidence. Leakage would push a permutation null **above** zero; the maximum
null mean observed is +0.079 (position blocks) and −0.159 (leave-one-gene-out).

## 7. Discussion

The useful result of this project is not a model. It is a table.

Aim 1's reclassification table changes what is known about 644 catalogued
mutations by supplying the one thing the catalogue structurally cannot contain:
how often each substitution occurs in virus that was never under drug pressure.
Eleven entries turn out to describe the reference sequence rather than the
virus. Seven should be reclassified. Forty-one are common enough that a
laboratory will meet them and has nothing to go on.

Aim 4 is worth reporting mainly for what it rules out. Population frequency does
not predict phenotype — it scores below its own no-skill floor — and that is not
a failure of Aim 1 but a statement about disjoint populations: 128 of 644
catalogued mutations are seen in circulating virus, but only 14 of 119
phenotyped ones. The two aims barely overlap, and conflating them is the easiest
mistake available with these results.

What does predict phenotype is structural context, and only where the structure
resolves a site. That is a narrower claim than "machine learning predicts HSV
resistance", and it comes with a built-in test for where it applies: if no
ligand defines the active site, the feature measures nothing and the model
should not be used. Helicase and primase are exactly that case, and they are
also the genes where new drugs are arriving.

Four methodological results here generalise beyond HSV, and may outlast the
biological ones. A grouped cross-validation scheme can manufacture negative
correlation from noise, so that a figure scored against a zero null inverts the
conclusion for any feature set without signal (§6.1). A conservation-matched
negative control is uninterpretable unless its random arm is matched on sample
size as well as on class balance (§6.2). A protein language model can clear
every significance test applied to it and still contribute nothing once
conservation and structural context are removed (§6.3). And a model can reach
an average precision of 0.916 while making no usable decision at the threshold
it was declared against (§5.4). Each was visible only because the metric
hierarchy, the operating point and the controls were fixed before any result was
seen, and each would have been invisible under the AUROC-only reporting common
in this literature.

## 8. Limitations

**Only 14 of 119 primary phenotyped variants are observed in circulating
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
being queried by an external review of this work as possibly misreported
relative to their primary sources; §3.1 resolves both against the catalogue's
own per-drug records, and the review was not itself a published audit.

**The temporal hold-out's cut-off year was not pre-declared.** All five
cut-offs are reported for that reason, and the nominated primary is chosen on
split balance rather than performance. Reference years come from Crossref and
are cached with their retrieval date; a study's year is the earliest among its
references.

**One pre-declared extension was not attempted, and the reason bounds what
this work can be used for.** The positive-unlabelled formulation would extend
scoring to variants of unknown significance — **717 variants observed in
circulating virus with no measured fold change**, 619 of them not catalogued at
all, against the 235 that carry a phenotype. That pool, not the labelled set,
is what a clinical laboratory actually encounters.

It is not attempted because three conditions are not met. The standard PU
assumption that positives are selected at random from all positives is
contradicted by this project's own central finding about ascertainment. The
class prior, on which PU performance almost entirely depends, cannot be
estimated stably from 87 positives. And the only feature group that works does
not cover the pool: 517 of the 717 lie in UL30 and UL52, and only 317 have a
ligand-resolved active site, so a score would be least trustworthy exactly
where it would be most used.

The 41 `unphenotyped_but_common` variants in the reclassification table are the
assumption-free substitute — a phenotyping-priority list derived from
observation rather than from a model. **No variant score in this work should be
used to interpret a patient isolate.**

**This is not a public pre-registration.** The analysis plan demonstrably
predates data acquisition by file timestamp, but it was not deposited with a
third party, and the evidence is under the authors' control.

## 9. Figures

**Figure 1. The catalogue against the population it is applied to.**
(**A**) Frequency distribution of all 644 catalogued substitutions in UL23,
UL30, UL5 and UL52 across 599 isolates, binned as in §3. 516 are absent from
circulating virus; 11 are the majority allele, meaning the reference sequence
carries the minority residue at that position. (**B**) Maximum measured fold
change against population frequency, both on log axes, for the 241
substitutions carrying a numeric phenotype. Substitutions never observed in
circulating virus are placed in the jittered strip at the left, outside the
frequency axis, because zero has no position on a log scale. The dashed
horizontal line is the pre-declared 2-fold threshold; the dotted vertical line
is 5% frequency. The shaded region is common (>5%) and resistant above 5×, and
it contains no substitution. Labelled points are those discussed in §3.1 and
§3.2.

**Figure 2. Purifying selection by gene.** pN/pS by Nei–Gojobori with
Jukes–Cantor correction, measured against the population consensus rather than
the reference, with 95% bootstrap intervals over 1,000 resamples. Genes with
insufficient synonymous sites are withheld rather than plotted with a wide
interval; HSV-2 UL23 (Sd = 32) is the only exclusion. The annotated contrast is
the pre-declared Aim 2 comparison.

**Figure 3. The structural signal, and its boundaries.** (**A**) Average
precision by virus and gene for the structure-only model and the frequency
baseline, each against that gene's positive-class prevalence, which is the
no-skill floor. Daggers mark the genes with no ligand-resolved active site,
where the distance feature has nothing to measure. (**B**) Average-precision
lift over the test-set floor for each feature group across the five temporal
hold-outs, training on substitutions first characterised before the cut-off year
and testing on those characterised after. Structure is the best feature group at
every cut-off; no alternative exceeds +0.163 at any split.

Figures are generated by `scripts/p4_08_figures.py` from the committed summary
files alone, as vector PDF and 300 dpi PNG, and every plotted value is written
to `reports/phase4/figures/figure_data.json` so that each is checked against its
source by the verification gate.

## 10. Availability

| Artefact | Path |
| --- | --- |
| Frozen accession list, with fates | `data/release/accession_list.tsv` |
| Master variant table | `data/release/master_variant_table.tsv` |
| **Reclassification table** (deposited resource) | `data/release/reclassification_table.tsv` |
| Data dictionary | `data/release/DATA_DICTIONARY.md` |
| Feature correlation matrix | `data/processed/feature_correlation_matrix.tsv` |
| Pre-declaration and compliance audit | `reports/phase4/PRE_DECLARATION.md` |
| Figures, and every value plotted in them | `reports/phase4/figures/` |
| Pipeline | `scripts/`, driven by `scripts/run_all.py` |

Checksums for the release tables are in `reports/phase4/release_summary.json`.

**Software.** Python 3.11.8 with NumPy (Harris et al. 2020), SciPy (Virtanen et
al. 2020), scikit-learn (Pedregosa et al. 2011), Biopython (Cock et al. 2009),
matplotlib (Hunter 2007), PyTorch and Hugging Face `transformers`. Pinned
versions are in the README.

## 11. References

**How these citations were obtained, and what they do and do not assert.**
Methodological and software references are to the works whose procedures this
pipeline implements. The HSV phenotype references are a different matter and are
labelled as such: each was taken from the `ref_doi` field of the frozen
HerpesDRG snapshot for the specific variant being discussed, not from an
independent literature search, and its bibliographic details were then resolved
against Crossref. Citing them here asserts that the catalogue attributes the
quoted value to that work — which is the claim §3.1 and §3.2 actually make, and
all they need. It does not assert that the present authors have examined the
primary measurements. Readers wishing to go behind the catalogue should treat
the DOI as the starting point, not as a warrant.

Agresti A, Coull BA (1998) Approximate is better than "exact" for interval
estimation of binomial proportions. *The American Statistician* 52:119–126.
doi:10.1080/00031305.1998.10480550

Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) Basic local alignment
search tool. *Journal of Molecular Biology* 215:403–410.
doi:10.1016/S0022-2836(05)80360-2

Berman HM, et al. (2000) The Protein Data Bank. *Nucleic Acids Research*
28:235–242. doi:10.1093/nar/28.1.235

Bestman-Smith J, et al. (2001) Highly reliable heterologous system for
evaluating resistance of clinical herpes simplex virus isolates to nucleoside
analogues. *Journal of Virology* 75:3105–3110.
doi:10.1128/JVI.75.7.3105-3110.2001 — *phenotype reference, via HerpesDRG*

Breiman L (2001) Random forests. *Machine Learning* 45:5–32.
doi:10.1023/A:1010933404324

Brunnemann AK, et al. (2016) Recombinant herpes simplex virus type 1 strains
with targeted mutations relevant for aciclovir susceptibility. *Scientific
Reports* 6:29903. doi:10.1038/srep29903 — *phenotype reference, via HerpesDRG*

Brunnemann AK, et al. (2018) Relevance of non-synonymous thymidine kinase
mutations for antiviral resistance of recombinant herpes simplex virus type 2
strains. *Antiviral Research* 152:53–57. doi:10.1016/j.antiviral.2018.02.004 —
*phenotype reference, via HerpesDRG*

Burrel S, et al. (2012) Impact of novel mutations of herpes simplex virus 1 and
2 thymidine kinases on acyclovir phosphorylation activity. *Antiviral Research*
96:386–390. doi:10.1016/j.antiviral.2012.09.016 — *phenotype reference, via
HerpesDRG*

Charles OJ, et al. (2024) HerpesDRG: a comprehensive resource for human
herpesvirus antiviral drug resistance genotyping. *BMC Bioinformatics* 25.
doi:10.1186/s12859-024-05885-5

Chibo D, et al. (2004) Molecular analysis of clinical isolates of acyclovir
resistant herpes simplex virus. *Antiviral Research* 61:83–91.
doi:10.1016/j.antiviral.2003.08.018 — *phenotype reference, via HerpesDRG*

Cock PJA, et al. (2009) Biopython: freely available Python tools for
computational molecular biology and bioinformatics. *Bioinformatics*
25:1422–1423. doi:10.1093/bioinformatics/btp163

Davis J, Goadrich M (2006) The relationship between precision–recall and ROC
curves. *Proceedings of the 23rd International Conference on Machine Learning*
233–240. doi:10.1145/1143844.1143874

Efron B (1979) Bootstrap methods: another look at the jackknife. *The Annals of
Statistics* 7:1–26. doi:10.1214/aos/1176344552

Gaudreau A, et al. (1998) Phenotypic and genotypic characterization of
acyclovir-resistant herpes simplex viruses from immunocompromised patients.
*Journal of Infectious Diseases* 178:297–303. doi:10.1086/515626 — *phenotype
reference, via HerpesDRG*

Harris CR, et al. (2020) Array programming with NumPy. *Nature* 585:357–362.
doi:10.1038/s41586-020-2649-2

Henikoff S, Henikoff JG (1992) Amino acid substitution matrices from protein
blocks. *Proceedings of the National Academy of Sciences* 89:10915–10919.
doi:10.1073/pnas.89.22.10915

Hoerl AE, Kennard RW (1970) Ridge regression: biased estimation for
nonorthogonal problems. *Technometrics* 12:55–67.
doi:10.1080/00401706.1970.10488634

Holm S (1979) A simple sequentially rejective multiple test procedure.
*Scandinavian Journal of Statistics* 6:65–70.

Hunter JD (2007) Matplotlib: a 2D graphics environment. *Computing in Science &
Engineering* 9:90–95. doi:10.1109/MCSE.2007.55

Inagaki T, et al. (2018) Acyclovir sensitivity and neurovirulence of herpes
simplex virus type 1 with amino acid substitutions in the viral thymidine kinase
gene, which were detected in the patients with intractable herpes simplex
encephalitis previously reported. *Japanese Journal of Infectious Diseases*
71:343–349. doi:10.7883/yoken.JJID.2018.176 — *phenotype reference, via
HerpesDRG*

Jukes TH, Cantor CR (1969) Evolution of protein molecules. In: *Mammalian
Protein Metabolism*, 21–132. doi:10.1016/B978-1-4832-3211-9.50009-7

Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version
7: improvements in performance and usability. *Molecular Biology and Evolution*
30:772–780. doi:10.1093/molbev/mst010

Lin Z, et al. (2023) Evolutionary-scale prediction of atomic-level protein
structure with a language model. *Science* 379:1123–1130.
doi:10.1126/science.ade2574

Matthews BW (1975) Comparison of the predicted and observed secondary structure
of T4 phage lysozyme. *Biochimica et Biophysica Acta* 405:442–451.
doi:10.1016/0005-2795(75)90109-9

Meier J, Rao R, Verkuil R, Liu J, Sercu T, Rives A (2021) Language models enable
zero-shot prediction of the effects of mutations on protein function. *Advances
in Neural Information Processing Systems* 34. Preprint:
doi:10.1101/2021.07.09.450648

Nei M, Gojobori T (1986) Simple methods for estimating the numbers of synonymous
and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution*
3:418–426. doi:10.1093/oxfordjournals.molbev.a040410

Pedregosa F, et al. (2011) Scikit-learn: machine learning in Python. *Journal of
Machine Learning Research* 12:2825–2830.

Roberts DR, et al. (2017) Cross-validation strategies for data with temporal,
spatial, hierarchical, or phylogenetic structure. *Ecography* 40:913–929.
doi:10.1111/ecog.02881

Saito T, Rehmsmeier M (2015) The precision–recall plot is more informative than
the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLOS
ONE* 10:e0118432. doi:10.1371/journal.pone.0118432

Schalkwijk HH, et al. (2024) Combination of ganciclovir and trifluridine
prevents drug-resistance emergence in HSV-1. *Antimicrobial Agents and
Chemotherapy* 68. doi:10.1128/aac.00110-24 — *phenotype reference, via
HerpesDRG*

Schubert A, et al. (2014) Single nucleotide polymorphisms of thymidine kinase
and DNA polymerase genes in clinical herpes simplex virus type 1 isolates
associated with different resistance phenotypes. *Antiviral Research* 107:16–22.
doi:10.1016/j.antiviral.2014.03.015 — *phenotype reference, via HerpesDRG*

Shrake A, Rupley JA (1973) Environment and exposure to solvent of protein atoms.
Lysozyme and insulin. *Journal of Molecular Biology* 79:351–371.
doi:10.1016/0022-2836(73)90011-9

Spearman C (1904) The proof and measurement of association between two things.
*The American Journal of Psychology* 15:72–101. doi:10.2307/1412159

Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed
solvent accessibilities of residues in proteins. *PLoS ONE* 8:e80635.
doi:10.1371/journal.pone.0080635

Varma S, Simon R (2006) Bias in error estimation when using cross-validation for
model selection. *BMC Bioinformatics* 7:91. doi:10.1186/1471-2105-7-91

Virtanen P, et al. (2020) SciPy 1.0: fundamental algorithms for scientific
computing in Python. *Nature Methods* 17:261–272.
doi:10.1038/s41592-019-0686-2

Wilson EB (1927) Probable inference, the law of succession, and statistical
inference. *Journal of the American Statistical Association* 22:209–212.
doi:10.1080/01621459.1927.10502953
