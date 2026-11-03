# Data Dictionary: Task 2 Accession Tables and HerpesDRG Catalogue

## 1. Accession Tables (`ul5_accession_table.tsv`, `ul52_accession_table.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `accession` | string | Sequence accession identifier (GenBank / ENA / DDBJ). |
| `virus` | string | Viral species (`HSV1` or `HSV2`). |
| `gene` | string | Target gene (`UL5` or `UL52`). |
| `source_db` | string | Discovery database (`NCBI_Nucleotide`, `ENA`, `DDBJ_Amenamevir`). |
| `country` | string | Geographic origin annotated in `/country` or `/geo_loc_name`. Empty if unannotated. |
| `collection_year` | string | Year of isolate collection. Empty if unannotated. |
| `host` | string | Natural host organism (e.g. `Homo sapiens`). Empty if unannotated. |
| `sequence_length` | integer | Length of parent sequence or CDS span in base pairs. |
| `CDS_start` | integer | 1-based start coordinate of the CDS. |
| `CDS_end` | integer | 1-based end coordinate of the CDS. |
| `has_full_CDS` | boolean | `True` if non-truncated, in-frame multiple of 3, no internal stop codons, and no assembly gaps. |
| `lab_strain` | boolean | `True` for laboratory-adapted passaged strains and in vitro drug-selected mutants. Must be excluded from natural population carrier frequency baselines. |
| `already_in_599` | boolean | `True` if accession was analyzed in the manuscript's 599-genome denominator. |
| `retrieved_utc` | string | ISO 8601 UTC timestamp of retrieval. |
| `url` | string | Web link to the sequence record. |

## 2. Active HerpesDRG Catalogue (`herpesdrg_ul5_ul52.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `gene` | string | Target gene (`UL5` or `UL52`). |
| `virus` | string | Target species (`HSV1` or `HSV2`). |
| `position` | integer | Amino acid position of the resistance substitution. |
| `ref_residue` | string | Reference amino acid residue (NCBI RefSeq). |
| `alt_residue` | string | Variant amino acid residue. |
| `phenotype_fold_change` | string/float | In vitro susceptibility fold-change (IC50 or EC50 ratio vs wild-type). |
| `drug` | string | Helicase-primase inhibitor drug tested (`Amenamevir` or `Pritelivir`). |
| `study_doi` | string | DOI of primary publication reporting the phenotype. |
| `notes` | string | Study notes or phenotypic annotations from HerpesDRG. |
