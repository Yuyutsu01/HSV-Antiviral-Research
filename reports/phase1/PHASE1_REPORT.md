# Phase 1 report — dataset construction

**Date:** 2026-11-02
**Accession lists frozen:** HSV-1 2026-11-01T20:02:36Z (597), HSV-2 2026-11-01T20:07:23Z (546)
**Reproduce with:** `p1_01_fetch.py` → `p1_02_gate2.py` → `p1_03_extract.py` → `p1_04_frequencies.py`

---

## Verdict

**Gate 2 passed at 100%. The frequency table is built, and it already shows the
effect the project was designed to measure.**

Of 644 labelled substitutions, 130 are observed in circulating virus. Eleven sit
at frequencies above 50%, fifteen more between 10% and 50%. Meanwhile the
mutations with large experimentally measured fold changes are rare, as they
should be. That contrast is the paper.

---

## Gate 2 — coordinate verification

| Gene | Match | Total | Rate |
| --- | --- | --- | --- |
| HSV-1 UL23 | 225 | 225 | 100% |
| HSV-1 UL30 | 277 | 277 | 100% |
| HSV-1 UL5 | 13 | 13 | 100% |
| HSV-1 UL52 | 7 | 7 | 100% |
| HSV-2 UL23 | 40 | 40 | 100% |
| HSV-2 UL30 | 76 | 76 | 100% |
| HSV-2 UL5 | 5 | 5 | 100% |
| **Total** | **643** | **643** | **100%** |

Threshold was 95%. Every reference residue named in a HerpesDRG label matches
the residue at that position in the NCBI reference protein. Reference lengths
resolved correctly (HSV-1 thymidine kinase 376 aa, polymerase 1235 aa), which
confirms the gene matcher selected the intended CDS rather than a paralog.

The single HSV-1 UL29 label has no reference protein in the target set and is
excluded from the denominator rather than counted as a failure.

---

## Dataset construction

| Stage | HSV-1 | HSV-2 | Total |
| --- | --- | --- | --- |
| Accessions in frozen list | 597 | 546 | 1,143 |
| Records returning any CDS | — | — | 975 |
| Records with at least one target gene | — | — | 707 |
| **Isolates retained after filtering** | — | — | **465** |

### Where the records went

**168 records return no coding sequences at all.** These are unannotated
submissions, and several are patent filings (`FV537006`–`FV537013`, `CS176883`)
rather than clinical or epidemiological isolates. Correctly excluded.

**A further 244 records carry CDS annotation but not for the target genes.**
Full annotation is present in roughly 710–730 records; the remainder annotate
only part of the genome.

This is worth stating plainly in the paper: **the effective population
denominator is about 710 genomes, not the 1,143 that a length-filtered NCBI
query returns.** Anyone quoting raw genome counts as a denominator is
overstating their sample by a third.

### De-duplication and exclusion

| Level | Removed |
| --- | --- |
| Identical target-gene sequence | 197 |
| Laboratory-adapted strain | 27 |
| Non-human host | 13 |
| Serial sample, same strain name | 5 |
| Within-study near-duplicate | 0 |

Non-human-host exclusion was added after the Phase 2 selection analysis
surfaced it (see below). Frequencies are claims about virus circulating in
humans, so an isolate recorded from another host does not belong in the
denominator. Records with a blank host field are kept, since most unannotated
submissions are human.

Laboratory strains caught include strain 17 (11 records), KOS (6), McKrae, F,
HG52 and MS. These are passaged stocks circulated between laboratories, not
independent natural isolates, and leaving them in would have inflated the
apparent frequency of whatever alleles they happen to carry.

The within-study level removed nothing, because identical sequences had already
been collapsed at level 1. Study linkage was resolved for 845 of 1,143 records
via PubMed ID or BioProject.

### Sequences per gene after QC

| | UL23 | UL30 | UL5 | UL52 |
| --- | --- | --- | --- | --- |
| HSV-1 | 251 | 253 | 256 | 253 |
| HSV-2 | 159 | 183 | 173 | 156 |

QC rejects sequences for ambiguous residues above 2%, length deviation beyond
10%, frame errors, and amino-acid divergence from the reference above 5%.

The divergence filter was added in Phase 2 after two HSV-2 records
(`KU310662`, `KU310667`) were found to carry 118 and 109 amino-acid differences
in UL23 — 31% divergence — while sitting inside the 10% length tolerance.
Within-species divergence in these genes is well under 2%, so 5% is generous
and removes only genuine outliers.

---

## Allele frequencies

1,048 variant rows across 7,118 reference positions. Variable positions range
from 3.6% of residues (HSV-2 UL5) to 27.1% (HSV-2 UL23).

**70% of variant rows are singletons** — observed in exactly one isolate. At
these sample sizes a singleton is indistinguishable from a sequencing or
assembly error, so headline claims apply a minimum count of two. The singleton
set is retained in the table and flagged, not deleted.

---

## The central result, arriving early

### Nineteen positions where the reference carries the minority allele

HSV-1 UL23 4 · HSV-1 UL30 3 · HSV-1 UL5 4 · HSV-1 UL52 2 · HSV-2 UL23 1 ·
HSV-2 UL30 3 · HSV-2 UL52 2.

At these positions the "mutation" recorded in the catalogue is what almost all
circulating virus carries, and the reference strain is the outlier. Eleven
labelled variants fall here. The most extreme is **HSV-1 UL23 A265T at 98.1%**
(254 of 259 isolates). An annotation that applies to 98% of circulating virus
carries no clinical information.

### Labelled variants by population frequency

| Band | Labelled variants |
| --- | --- |
| Above 50% — reference is the minority allele | 11 |
| 10% to 50% — common natural polymorphism | 15 |
| 1% to 10% | 35 |
| Below 1% | 69 |

The 10–50% band is the more interesting of the two upper bands. These are not
reference artefacts; they are genuinely polymorphic positions where a
substitution catalogued as resistance-associated is carried by a substantial
minority of ordinary isolates.

### Tier A mutations — the cases with experimental phenotypes

These carry engineered, phenotyped fold changes, so frequency and phenotype can
be read together.

| Variant | Frequency | Max fold change | Reading |
| --- | --- | --- | --- |
| HSV-2 UL23 G39E | 0.792 | 2.4× | Majority allele, mild measured resistance. The strongest single reclassification candidate |
| HSV-1 UL30 P920S | 0.252 | 1.05× | Common polymorphism, no resistance. Correctly a negative |
| HSV-1 UL23 R41H | 0.046 | 2.3× | Borderline on both axes |
| HSV-1 UL30 L802F | 0.008 | 7.4× | Genuine resistance, appropriately rare |
| HSV-1 UL30 T821M | 0.004 | 6.25× | Genuine resistance, appropriately rare |
| **HSV-1 UL30 S724N** | **0.004** | **72×** | High-level resistance, one isolate in 258 |

**The pattern holds in the right direction.** Mutations with large measured
fold changes are rare in circulating virus; the catalogued entries that are
common have small or absent measured effects. That is exactly what should
happen if the catalogue is a mixture of genuine resistance determinants and
misclassified natural variation — and it is the first direct evidence for the
project's premise.

---

## Two limitations that must be reported

### Treatment-exposure stratification is not possible from metadata

The protocol specified reporting frequencies stratified by a treatment-naive
proxy. **Zero of 1,143 records could be classified as treatment-associated.**
GenBank source features for HSV genomes essentially never record antiviral
exposure: strain is populated for 738 records, host for 978, geographic
location for 709, collection date for 490 — but treatment history for none.

The stratification therefore cannot be done, and claiming it would be false.
What partly substitutes is the decision taken in Phase 0 to exclude targeted
single-gene amplicons from the denominator. Those submissions are the ones
enriched for treatment failure, and they are not in this dataset. Complete
genomes are predominantly submitted from epidemiological and phylogenetic
studies rather than from resistance testing.

The honest wording for the paper: frequencies are estimates **in sequenced
complete genomes**, which are not a random sample of infections and whose
treatment history is unknown.

### Reference-anchored pairwise alignment, not a multiple alignment

No MAFFT, MUSCLE or BLAST on the analysis machine. Each isolate protein is
aligned pairwise to the reference with BLOSUM62, affine gaps, and free terminal
gaps so partial CDS are not penalised. Within a species the sequences are
98–99% identical, so this recovers the same columns an MSA would while mapping
every residue directly onto a reference coordinate — which is the property that
matters here, because the labels are indexed against that same reference.

**Phase 2 still needs a genuine codon-aware multiple alignment** for the dN/dS
selection analysis in Aim 2. That is a real dependency, not an optional
refinement.

---

## Bugs found and fixed during this phase

Recorded because each produced a pipeline that ran cleanly and silently
produced wrong output.

1. **Metadata parsing returned nothing.** Truncated GenBank records print
   `source          1` with no `a..b` range, so a regex expecting a range
   matched nothing. All source qualifiers came back empty while the script
   reported success. Caught because de-duplication removed zero laboratory
   strains, which is not credible for HSV.
2. **UL8 could be captured as UL5.** The protein-name pattern `dna helicase`
   matched "DNA helicase/primase complex-associated protein", which is UL8.
   Patterns now require the full subunit name, and the ambiguous bare
   "helicase-primase subunit" is left unmatched rather than guessed.
3. **"Variable positions" counted covered positions.** A summary statistic
   reported every covered position as variable. The underlying table was
   correct; only the summary was wrong.

---

## Files produced

| Path | Contents |
| --- | --- |
| `data/raw/accessions/*.json` | Frozen accession lists with query terms and dates |
| `data/raw/reference/*.faa` | Reference proteins for the four target genes |
| `data/interim/isolate_metadata.tsv` | Per-accession strain, host, location, date, study |
| `data/interim/cds_qc.tsv` | Per-sequence QC verdict and flags |
| `data/processed/genes/*.faa` | QC-passed, de-duplicated protein sequences |
| `data/processed/gate2_coordinate_check.tsv` | Per-label coordinate verdict |
| `data/processed/allele_frequencies.tsv` | 1,048 variants with Wilson intervals |
| `reports/phase1/*.json` | Machine-readable summaries |

---

## Phase 2 entry conditions

Gate 3 asks whether any labelled variant exceeds 1% in naive isolates. **It does
— 61 of them.** The paper does not become a confirmation study. Phase 2 proceeds
to the reclassification table, the UL23/UL30 selection contrast, and the
helicase-primase baseline.
