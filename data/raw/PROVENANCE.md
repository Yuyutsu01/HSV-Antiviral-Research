# Raw data provenance

## herpesdrg-db.tsv
- Source: https://github.com/ojcharles/herpesdrg-db
- Commit SHA (pinned): 22fc18676a1ea49cd6d7cd1c3280db210a7580e6
- Permalink: https://raw.githubusercontent.com/ojcharles/herpesdrg-db/22fc18676a1ea49cd6d7cd1c3280db210a7580e6/herpesdrg-db.tsv
- Retrieved (UTC): 2026-11-01T19:30:05Z
- MD5: d78ad69a5879d854a2d1b22a2698ae46
- Citation: Charles OJ, Venturini C, Goldstein RA, Breuer J. HerpesDRG: a comprehensive
  resource for human herpesvirus antiviral drug resistance genotyping.
  BMC Bioinformatics 2024. doi:10.1186/s12859-024-05885-5
- License: MIT (see repo LICENSE)

---

## Second catalogue: Dähne et al. 2025 (Clinical Microbiology and Infection)

**Citation:** Dähne T, Jaki L, Gosert R, Fuchs J, Krumbholz A, Nägele K,
Pletz MW, Khanna N, Leuzinger K, Panning M. *Herpes simplex virus and drug
resistance — comprehensive update on resistance mutations and implications for
clinical management: a narrative review.* Clin Microbiol Infect.
2025;31(11):1484–1490. doi:10.1016/j.cmi.2025.04.046

**Licence:** CC BY 4.0 (open access). Reuse permitted with attribution.

**Files supplied by the user, not re-downloaded:**
- `mmc3.xlsx` — Table S3, the full dataset. Six sheets
  (HSV-1/HSV-2 × TK, Pol, HPC), 1,594 usable drug × mutation rows.
- Table S2 (pasted as text) — the list of corrections applied to the
  Sauerbrei et al. dataset, ~130 entries.

**Why it matters that the coordinates match.** The paper states that HSV-1
mutations were aligned to strain 17 **NC_001806** and HSV-2 mutations to HG52
**NC_001798**. Those are the same two references this project uses
(`REFSEQ` in `p1_02_gate2.py`), so positions are directly comparable with no
remapping and no risk of an off-by-one from a different start codon.

**How their labels must be read.**
- *Resistance status* is "recorded based on author reports" — transcribed from
  the original publications, not re-measured. It is a second curation, not an
  independent gold standard, and it inherits each original author's EC50
  cut-off.
- *Resistance confirmation in cell culture*: 1 if phenotypic testing was
  explicitly performed, else 0.
- *Clinical correlation*: 1 if treatment failure aligned with phenotypic
  testing, **2 for contradictory results**, 0 if no clinical data. The
  "contradiction" values in the spreadsheet are this code 2.
- Labels are per drug. Collapsing a substitution across drugs loses real
  information, since cross-resistance is partial (TK mutations do not confer
  foscarnet or cidofovir resistance).

**Their own discordance result**, which is complementary to this project's:
of 513 mutation–drug combinations having both phenotypic and clinical data,
461 (89.9%) agreed and 52 (10.1%) disagreed — and **every disagreement was a
phenotypically susceptible isolate in a patient assessed as clinically
resistant.** That is the opposite direction to the failure mode studied here
(catalogued resistance alleles that are common natural polymorphisms), so the
two findings bracket the problem rather than duplicate it.

**Known discrepancy.** The paper's text gives per-gene combination counts of
451 (HSV-1 TK), 722 (HSV-1 POL), 12 (HSV-1 HP), 287 (HSV-2 POL) and 15
(HSV-2 HP), which match the spreadsheet to within one row. For HSV-2 TK the
text reads "212" while the sheet carries 117 data rows. The sheet is treated as
authoritative here; the discrepancy is noted rather than resolved.

**Scope note for writing up.** This dataset already labels many variants as
polymorphisms or as Sensitive. For those, this project's contribution is a
population frequency they do not provide — not the discovery that the variant
is benign. Claims must be worded accordingly.
