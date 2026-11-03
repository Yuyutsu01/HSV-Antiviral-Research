# Data Dictionary: Task 1 Accession Table

Table file: `accession_table.tsv`

| Column | Type | Description |
| :--- | :--- | :--- |
| `accession` | string | Primary sequence accession (GenBank/ENA/INSDC). |
| `source_db` | string | Database(s) through which the record was retrieved (`NCBI_Nucleotide`, `ENA`, `NCBI_PopSet`, `DDBJ`, `GenBank_Explicit_OR`). Semicolon-delimited if multi-source. |
| `country` | string | Geographic origin as explicitly annotated in the `/country` feature qualifier. Empty if unannotated. |
| `collection_year` | string | Year or date of isolate collection as annotated in `/collection_date` qualifier. Empty if unannotated. |
| `host` | string | Natural host organism as annotated in `/host` qualifier (e.g. `Homo sapiens`). Empty if unannotated. |
| `sequence_length` | integer | Total nucleotide length of the parent sequence record in base pairs. |
| `UL23_CDS_start` | integer | 1-based start coordinate of the UL23 coding sequence in the parent sequence. |
| `UL23_CDS_end` | integer | 1-based end coordinate of the UL23 coding sequence in the parent sequence. |
| `has_full_CDS` | boolean | `True` if feature is non-truncated, starts with ATG, ends with canonical stop, has no internal stop codons, and length >= 900 nt. |
| `retrieved_utc` | string | ISO 8601 UTC timestamp of retrieval. |
| `url` | string | Direct web URL to the sequence record. |
