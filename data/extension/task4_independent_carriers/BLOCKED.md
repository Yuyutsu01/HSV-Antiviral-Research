# Sub-Query Audit and Blocked Access Notes: Task 4 (Independent Carriers)

## 1. R176Q: Single-Study Dependency in Public Sequence Repositories
- **Query / Endpoint**: Exhaustive scan across all 1,261 HSV-1 UL23 sequence records in NCBI Nucleotide and ENA.
- **Finding**: Zero independent sequence accessions carrying R176Q were identified outside the single dominant study `PRJNA562737` (Shipley et al. 2019, accessions `MN401201.1`, `MN401205.1`, `MN401206.1`, `MN401207.1`).
- **Literature Audit**: Historical virological papers reporting R176Q (Chatis et al. 1991, Saijo et al. 2001, Frobert et al. 2005, van der Beek et al. 2004) characterized isolated strains phenotypically but did not deposit accessioned nucleotide sequences into INSDC databases.
- **Outcome**: The success criterion (>= 1 independent sequence carrier from another study) was **NOT MET**. R176Q remains completely study-dominated (`max_study_share = 1.0`), and the negative finding stands.

## 2. R41H: Universal Absence of High-Level Resistance Phenotype
- **Query / Endpoint**: PubMed and virological database search for `"R41H"` and `("fold change" OR "IC50" OR "EC50" OR "resistance ratio")`.
- **Finding**: Both gold-standard reverse-genetics studies using recombinant BAC-derived clones (Bohn-Wippert et al. 2016, Suzuki et al. 2018) measure fold changes of 1.0x and 0.8x respectively (identical to wild-type sensitivity). Zero studies in the literature report a fold change above 3x for an isolated R41H substitution.
- **Outcome**: The success criterion (>= 1 independent phenotype above 3x) was **NOT MET**. The evidence conclusively proves that R41H is a susceptible circulating polymorphism, not a circulating resistance allele; the negative stands.

## 3. R222C: Broad Population Spread Confirmed
- **Query / Endpoint**: Scan across 1,261 HSV-1 UL23 records.
- **Finding**: Identified 27 additional sequence carriers across 3 independent studies (PMIDs 25926637, 25855744, 37319284), raising the carrier total from 11 to 38 and reducing maximum single-study share from 72.7% to 21.1% (< 0.50).
- **Outcome**: The success criterion for R222C is fully **MET**.
