# Data dictionary — release tables

Three files, all tab-separated, UTF-8, with a header row. SHA-256 checksums and
row counts are in `reports/phase4/release_summary.json`.

| File | Rows | What it is |
| --- | --- | --- |
| `accession_list.tsv` | 1,143 | Every GenBank accession the frozen query returned, with its fate |
| `master_variant_table.tsv` | 1,162 | One row per unique variant, catalogued or observed |
| `reclassification_table.tsv` | 644 | The deposited resource: every catalogued HSV mutation with its population frequency |

Reference coordinates are **NC_001806** (HSV-1) and **NC_001798** (HSV-2)
throughout. Positions are 1-based amino-acid positions in the reference protein.

---

## `accession_list.tsv`

The frozen accession list. It records what was **excluded** as well as what was
kept, because a list of only the 451 analysed isolates is not reproducible: a
reader could not tell whether an absent accession was missed by the query or
dropped by a filter.

| Column | Type | Meaning |
| --- | --- | --- |
| `accession` | text | GenBank nucleotide accession |
| `virus` | HSV1 / HSV2 | Assigned by the taxonomic query, not by sequence |
| `fate` | enum | See below |
| `fate_detail` | text | Reason specifics: QC flags, the retained duplicate, the host, the strain name |
| `n_genes_analysed` | int | Of the four target genes, how many contributed variant calls |
| `genes_analysed` | text | Semicolon-separated gene names |
| `study_key` | text | Submission grouping used for leave-one-study-out; empty where unrecoverable |
| `pubmed` | text | PubMed ID where the record carries one |
| `collection_year`, `country`, `host` | text | As annotated in GenBank, not normalised |
| `is_lab_strain` | 0/1 | Laboratory-strain heuristic |
| `strain` | text | Strain or isolate name as annotated |
| `query_term`, `query_date` | text | The frozen query and when it ran |

### `fate` values

| Value | n | Meaning |
| --- | --- | --- |
| `analysed` | 451 | Contributed variant calls to at least one gene |
| `excluded_no_cds_extracted` | 411 | No target-gene CDS could be extracted; usually an unannotated or partial record |
| `excluded_downstream_filter` | 183 | Had a passing CDS but did not reach the coverage table — divergence and per-gene coverage filters |
| `excluded_cds_qc` | 36 | CDS extracted but failed QC; `fate_detail` names the flags |
| `excluded_laboratory_strain` | 30 | Laboratory-adapted strain |
| `excluded_patient_duplicate` | 18 | Serial isolate from a patient already represented; `fate_detail` names the retained record |
| `excluded_non_human_host` | 14 | Non-human host annotation |

**`fate` is reconstructed, not logged at the moment of exclusion.** It is
derived by replaying the filters against the intermediate tables, so the
category is reliable but a record failing two filters is attributed to whichever
is checked first. The order is: analysed → patient duplicate → laboratory strain
→ non-human host → CDS QC → no CDS → downstream. `excluded_downstream_filter` is
the least specific category and should be read as "passed CDS QC, excluded
later" rather than as a single cause.

---

## `master_variant_table.tsv`

One row per unique `(virus, gene, position, ref_aa, alt_aa)`. This is the union
of the catalogue and the observed population, **not their intersection** — that
disjunction is the project's central finding, so either `in_catalogue` or
`observed` may be 0.

| 1,162 rows | |
| --- | --- |
| in the catalogue | 644 |
| observed in circulating virus | 638 |
| **both** | **120** |
| catalogue only | 524 |
| observed only | 518 |
| with a measured phenotype (in the model dataset) | 235 |
| with structural features | 896 |
| with ESM-2 scores | 1,161 |

### Columns

**Identity:** `virus`, `gene`, `position`, `ref_aa`, `alt_aa`, `aa_change`.

**Catalogue side** — empty when `in_catalogue = 0`:

| Column | Meaning |
| --- | --- |
| `in_catalogue` | 0/1, present in HerpesDRG with `status == "A"` |
| `evidence_tier` | `A_engineered` (133), `B_lab_selected` (66), `C_clinical_isolate` (445). Strongest tier available for that mutation |
| `n_refs` | Independent references |
| `drugs_tested` | Semicolon-separated |
| `max_fold_change` | Maximum EC50 fold change across drugs and studies. Pooled, so it is an upper bound, not a per-drug value |
| `resistant_at_2x` | 0/1 at the declared threshold; empty where no numeric fold change exists |

**Population side** — empty when `observed = 0`:

| Column | Meaning |
| --- | --- |
| `observed` | 0/1, seen in at least one analysed isolate |
| `count` | Isolates carrying the allele, after patient-level collapse |
| `coverage` | Isolates with a called residue at that position |
| `frequency` | `count / coverage` |
| `ci_low`, `ci_high` | Wilson score interval, 95% |

**Ascertainment** — the columns that decide whether a frequency is real:

| Column | Meaning |
| --- | --- |
| `n_studies_carrying` | Distinct submissions carrying the allele |
| `loso_min_freq` | Lowest frequency under leave-one-study-out |
| `max_study_share` | Largest share of carriers from a single study. 1.00 means every carrier came from one submission |
| `study_dominated` | 1 if dropping one study cuts the frequency by more than half |
| `frequency_band`, `flag` | See the reclassification table below |

**Features** — present wherever computable, independently of phenotype:

| Column | Meaning |
| --- | --- |
| `cons_site_variability` | Non-synonymous fraction at the site, measured against the **population consensus**, not the reference |
| `cons_site_coverage` | Isolates contributing |
| `struct_rsa` | Relative solvent accessibility, Shrake–Rupley over Tien et al. maxima |
| `struct_dist_to_ligand` | Ångström to the nearest bound ligand heavy atom |
| `struct_site_resolved` | 0 if the structure has no ligand-defined site for that gene. **Always 0 for UL52** — its nearest ligand belongs to a neighbouring subunit, median distance 42 Å, so `struct_dist_to_ligand` is not an active-site measure there |
| `struct_pdb_id` | Source entry: 1vtk (UL23), 8exx (UL30), 9vlq (UL5), 9nnp (UL52) |
| `struct_inherited` | 1 for HSV-2 rows, whose geometry is inherited from HSV-1 structures by pairwise alignment because no HSV-2 structures exist for these genes |
| `esm2_masked_marginal` | ESM-2 650M masked-marginal score, `log P(alt) − log P(ref)` at the masked position |
| `esm2_ref_logprob` | Reference-residue log probability |
| `in_model_dataset` | 1 if the variant carried a numeric fold change and entered Aim 4 |

---

## `reclassification_table.tsv`

The deposited resource. Every catalogued HSV mutation in the four target genes,
with its measured population frequency and the ascertainment checks needed to
decide whether that frequency means anything. Same columns as the catalogue and
population sides above, plus:

### `frequency_band`

| Value | n | Definition |
| --- | --- | --- |
| `not_observed` | 524 | Absent from all analysed isolates |
| `singleton` | 47 | Exactly one carrier |
| `rare` | 11 | Below 1% |
| `uncommon_polymorphism` | 36 | 1% to 5% |
| `common_polymorphism` | 15 | 5% to 50% |
| `majority_allele` | 11 | Above 50% — the reference carries the minority residue |

### `flag`

The interpretive verdict. **These are hypotheses for review, not clinical
determinations.**

| Value | n | Meaning |
| --- | --- | --- |
| `not_observed` | 524 | No population evidence either way |
| `consistent_rare` | 53 | Rare in circulation, as a resistance determinant should be |
| `unphenotyped_but_common` | 40 | Common in circulation, no measured fold change. The group most in need of phenotyping |
| `reference_artifact` | 10 | Majority allele with no measured resistance: the reference sequence, not the virus, is the outlier |
| `reclassify_as_polymorphism` | 8 | Catalogued as resistance-associated, common in circulation, fold change below threshold |
| `consistent_with_resistance` | 5 | Rare and phenotypically resistant |
| `circulating_resistance_allele` | 3 | Common **and** phenotypically resistant. **All three fail ascertainment checks** — see below |
| `majority_allele_with_measured_resistance` | 1 | Majority allele with a measured fold change above threshold |

### Read `circulating_resistance_allele` with the ascertainment columns

These three are the most consequential rows in the table, and none of them
supports a claim of circulating resistance — but **they fail for two different
reasons, and only two of the three fail on ascertainment**:

| Variant | Freq | Carriers | Studies | `max_study_share` | `loso_min_freq` | Fold change | Why it does not stand |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HSV-1 UL23 R222C | 0.046 | 11 | 3 | 0.73 | 0.0135 | 98× | 8 of 11 carriers from one submission; frequency falls 3.4× on dropout |
| HSV-1 UL23 R176Q | 0.017 | 4 | 1 | 1.00 | 0.0 | 200× | Every carrier from a single submission; drops to zero on dropout |
| HSV-1 UL23 R41H | 0.051 | 12 | 5 | 0.33 | 0.0385 | 2.3× | **Ascertainment is sound.** It fails on phenotype: 2.3× is marginal against the 2× threshold, and an independent published audit labels it Sensitive |

R222C and R176Q are ascertainment artefacts, and `study_dominated = 1` catches
both. R41H is not an artefact — it is a genuine 5% polymorphism seen across five
submissions — and `study_dominated = 0` correctly passes it. What disqualifies
it is that its "resistance" rests on a fold change barely above the declared
threshold, which is a question about the catalogue's phenotype, not about
sequence ascertainment.

So a flag of `circulating_resistance_allele` is an invitation to check both
sides: `max_study_share` and `loso_min_freq` for whether the frequency is real,
and `max_fold_change` against the 2× threshold for whether the phenotype is.
Filtering on `study_dominated = 0` alone would have retained R41H.
