# Phase 0 report — scoping and feasibility

**Date:** 2026-11-02
**Data source:** `ojcharles/herpesdrg-db`, pinned at commit `22fc18676a1ea49cd6d7cd1c3280db210a7580e6`
**Reproduce with:** `python scripts/phase0_audit.py`

---

## Verdict

**All Phase 0 gates pass. The project proceeds as specified, with three changes to the plan.**

The HSV label budget is substantially larger than the HerpesDRG publication's summary
("1341 mutations across 17 genes, 6 HCMV and 4 HSV-1") implied, and the database
contains an experimentally validated negative set that the protocol assumed did not exist.

---

## Gate 1 — label budget

| | Count |
| --- | --- |
| Rows in file | 2,035 |
| Rows with `status == "A"` (active) | 1,696 |
| Active HSV rows | 951 |
| **Unique HSV mutations** | **700** |

Row counts overstate labels: one mutation appears once per drug, reference and assay.
The analysis unit is the unique `(virus, gene, aa_change)` triple.

### Unique mutations by evidence tier

Tier is assigned from `strain_generation`, taking the strongest tier available
for each mutation.

| Tier | Source values | Total | HSV-1 | HSV-2 | Independent refs |
| --- | --- | --- | --- | --- | --- |
| **A — engineered** | Recombinant BAC, Marker Transfer | **133** | **103** | **30** | 61 |
| **B — lab-selected** | Lab induced / derived mutant, drug induced | 89 | 75 | 14 | 31 |
| **C — clinical isolate** | Isolated strain | 478 | 378 | 100 | 42 |

**Gate 1 threshold was ≥50 Tier A HSV-1 mutations. Observed: 103. PASS.**
All 103 are substitutions, so the gate also passes under the substitution-only
modelling scope.

### Tier A by gene

| Virus | Gene | Tier A |
| --- | --- | --- |
| HSV-1 | UL30 (DNA polymerase) | 49 |
| HSV-1 | UL23 (thymidine kinase) | 40 |
| HSV-1 | UL5 (helicase) | 7 |
| HSV-1 | UL52 (primase) | 6 |
| HSV-1 | UL29 (ICP8) | 1 |
| HSV-2 | UL30 | 18 |
| HSV-2 | UL23 | 12 |

Five HSV-1 genes are covered, not four — UL29 has been added since publication.

---

## Finding 1 — the phenotype is quantitative, and negatives exist

The thirteen drug columns hold **numeric fold-change in EC50 relative to a control
strain**, not resistant/susceptible calls. 122 of 133 Tier A mutations carry at least
one numeric value.

| Fold-change | Tier A mutations |
| --- | --- |
| ≤ 2× (no meaningful resistance) | **33** |
| > 2× | 89 |
| > 5× | 56 |

Range 0.90× to 5000×, median 4.7×.

**This is the most consequential finding of Phase 0.** Those 33 mutations were
deliberately engineered into a clean background, phenotyped, and found *not* to confer
resistance. They are experimentally validated negatives — the class the protocol
assumed was unavailable and planned to approximate with conservation-matched
population variants.

Consequences for the plan:

- The supervised arm gains a true negative set, so Aim 4 no longer rests entirely on
  proxy negatives. Proxy negatives remain useful as a second, larger comparison set.
- The primary task can be framed as **regression on log fold-change** rather than
  binary classification, which uses the full information and avoids an arbitrary
  threshold. Binary classification at a declared 2× cutoff becomes the secondary
  analysis.
- The circularity audit becomes stronger: engineered positives and engineered
  negatives come from the same experimental process, so comparing them controls for
  publication and ascertainment effects that proxy negatives cannot.

### Drug coverage, Tier A

Aciclovir 132 · Foscarnet 108 · Cidofovir 61 · Penciclovir 59 · Brivudine 28 ·
Amenamevir 15 · Ganciclovir 13 · Pritelivir 6

Aciclovir and foscarnet are the only drugs with enough Tier A labels to model
independently. The helicase-primase inhibitors have 21 labels between them — enough
for the descriptive baseline of Aim 2, not for a per-drug model.

---

## Finding 2 — HSV-2 is not underpowered at the genome level

NCBI nucleotide counts, retrieved 2026-11-02:

| Query | Count |
| --- | --- |
| HSV-1, all nucleotide records | 10,864 |
| **HSV-1, 140–200 kb (complete genomes)** | **597** |
| HSV-2, all nucleotide records | 8,955 |
| **HSV-2, 140–200 kb (complete genomes)** | **546** |
| HSV-1 UL23 targeted amplicons (0.5–5 kb) | 1,335 |
| HSV-1 UL30 targeted amplicons | 973 |
| HSV-2 UL23 targeted amplicons | 349 |
| HSV-2 UL30 targeted amplicons | 277 |

HSV-2 complete genomes are at near parity with HSV-1. The comparative arm is
therefore **viable rather than exploratory**, and HSV-2 frequency estimates can carry
real confidence intervals.

The targeted amplicons are a trap and must be excluded from the primary frequency
denominator. Single-gene UL23 and UL30 submissions come predominantly from
resistance-testing studies, so they are enriched for treatment failure by
construction — exactly the ascertainment bias the protocol warns about. They are kept
as a clearly separated secondary stratum.

---

## Finding 3 — variant types, and the indel question resolved

Across all 700 unique HSV mutations:

| Type | n | % |
| --- | --- | --- |
| Substitution | 644 | 92.0 |
| Nonsense (premature stop) | 30 | 4.3 |
| Indel / frameshift | 13 | 1.9 |
| Other | 13 | 1.9 |

The concern that indel blindness would cripple the clinical interpretation was
**overstated**. Substitutions are 92% of the catalogue, and every Tier A HSV-1 entry
is a substitution.

Truncating events are nonetheless real and concentrated exactly where expected:

- HSV-1 UL23 — 17 nonsense, 4 deletions, 9 frameshifts (30 events)
- HSV-2 UL23 — 10 nonsense, 4 frameshifts
- HSV-2 UL30 — 6 deletions, mostly in-frame (`del676-677`, `del681-682`, …)

The in-frame UL30 deletions are notable: a polymerase cannot tolerate frameshift, so
deletions there are short and in-frame, while thymidine kinase — being dispensable —
tolerates full inactivation. That contrast is direct evidence for the Aim 3 constraint
asymmetry and should be reported as such rather than discarded as out of scope.

---

## Finding 4 — the co-mutation column cannot be used as a confounder filter

| Group | Rows flagged with `co_gene`/`co_aa` |
| --- | --- |
| Engineered | 0 / 160 (0.0%) |
| Lab-selected | 3 / 149 (2.0%) |
| Clinical isolate | 5 / 642 (0.8%) |

A 0.8% co-mutation rate among clinical isolates is not credible as a true rate —
clinical isolates carrying resistance mutations routinely carry other substitutions.
The column is **sparsely populated rather than reliably negative**, so it cannot be
used to identify clean single-mutation clinical entries.

Tier C therefore stays unresolvable for attribution. This strengthens the case for
treating Tier A as the analysis set and Tier C as context only.

---

## Finding 5 — coordinates and structures are not a risk

Two risks flagged in the protocol are now closed:

**Coordinates.** The database README states that all mutations are mapped to the NCBI
reference strain, so the catalogue is internally harmonised. The Phase 1 gate (>95%
reference-residue agreement) still runs, but as confirmation rather than as a
repair step.

**Structures.** All four primary targets have experimental structures in the PDB; no
AlphaFold models are required.

| Gene | PDB entries | Examples |
| --- | --- | --- |
| UL23 thymidine kinase | 35 | 1KIM, 1VTK, 1KI6 |
| UL30 DNA polymerase | 16 | 8EXX, 9JA3, 9ENQ |
| UL5 / UL52 helicase-primase | 12–14 | 9VLQ, 9YCP, 9NNP |
| UL42 processivity factor | 12 | 1DML, 8OJ6 |

---

## Three changes to the protocol

1. **Primary task becomes regression on log fold-change**, with binary classification
   at a declared 2× threshold as the secondary analysis. Reason: the phenotype is
   quantitative and 122 Tier A mutations carry numeric values.

2. **The negative set becomes experimental rather than proxy.** 33 engineered
   non-resistant mutations form the primary negative class; conservation-matched
   population variants become the secondary comparison. The matched-negative gap is
   still reported.

3. **HSV-2 is promoted from exploratory to a full comparative arm**, on the strength
   of 546 complete genomes and 30 Tier A mutations.

Unchanged: the substitution-only modelling scope, the evidence tiering, the
pre-declared metric and threshold, and the commitment to report a negative Aim 4
result as the finding.

---

## Finding 6 — CHARMD holds no HSV content, and its site is currently offline

Checked 2 November 2026. CHARMD is the Comprehensive Herpesviruses Antiviral drug
Resistance Mutation Database, published in *Antiviral Research* in 2024 by Tilloy and
colleagues at the French CNR Herpesvirus in Limoges.

**Its web interface is down.** Every endpoint under `/cnr-herpesvirus/outils/`
returns HTTP 404, including the entry point the CNR's own "Base de données
résistances antiviraux" page links to. The whole tools tree is gone, not just one
page. Whether this is a migration or a lapse is not determinable from outside.

**It has no HSV data.** A CNR presentation dated October 2025, still hosted on the
site, lists under *Future developments*:

- "Increase the CHARMD working group to integrate HSV and other herpesviruses resistance"
- "Databases for HSV1/2 (+ E Frobert, D Boutolleau, G Andrei...)"
- "ASK for additional funding european? NIH?" and, as the closing line, "SEARCH FOR SPONSORING"

So thirteen months ago HSV was aspirational and explicitly unfunded. Nothing
indicates that changed, and the database is for cytomegalovirus in practice.

### Consequences

**The cross-check source does not exist.** The protocol planned to validate Tier A
labels against CHARMD. That is not available. The Jena compilation becomes the only
external cross-check, and it overlaps HerpesDRG substantially, so label validation
rests mainly on HerpesDRG's own provenance tiering. This should be stated as a
limitation rather than left implicit.

**The timing risk is lower than assumed.** HSV expansion is unfunded and at the
seeking-sponsors stage. The named collaborators are established HSV resistance
researchers, so the group could move quickly if funded — but their stated object is
curating resistance mutations, which is not the population-frequency contribution
this project makes. The gap is not being closed by anyone at present.

**Expect poor concordance if CHARMD HSV ever appears.** A 2025 *Frontiers in
Microbiology* study comparing the two databases for cytomegalovirus reported that
none of its confirmed UL54 mutations were found in both CHARMD and HerpesDRG. Two
expert curations of the same literature disagreeing that completely is itself a
finding about the state of resistance annotation, and worth a sentence in the
discussion.

**All Phase 0 items are now closed.**

---

## Files produced

| Path | Contents |
| --- | --- |
| `data/raw/herpesdrg-db.tsv` | Pinned source file |
| `data/raw/PROVENANCE.md` | Commit SHA, MD5, retrieval date, citation |
| `data/processed/hsv_mutations_tiered.tsv` | 700 unique HSV mutations with tier, parsed position, drugs, max fold-change |
| `reports/phase0/phase0_summary.json` | All counts as machine-readable JSON |
| `scripts/phase0_audit.py` | Reproduces every number above; gates as assertions |
