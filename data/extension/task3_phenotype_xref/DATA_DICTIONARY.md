# Data Dictionary: Task 3 Phenotype Cross-Reference

## 1. Cross-Reference Results (`xref_results.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `variant_id` | string | Unique identifier formatted as `{virus}_{gene}_{aa_change}`. |
| `gene` | string | Target viral gene (`UL23` or `UL30`). |
| `ref` | string | Reference amino acid (NCBI RefSeq). |
| `pos` | integer | Amino acid position coordinate. |
| `alt` | string | Variant amino acid. |
| `found_in` | string | Semicolon-delimited list of databases/sources with records for this variant (`HerpesDRG`, `PubMed`, `None`). |
| `fold_change` | string/float | Quantitative in vitro susceptibility fold-change (IC50 or EC50 ratio vs wild-type). Empty if no numeric measurement exists. |
| `source_url_or_doi` | string | Semicolon-delimited primary DOIs or database URLs. |
| `confidence` | string | Confidence level (`high`, `medium`, `low`) regarding the presence or absence of a validated phenotype. |
| `notes` | string | Detailed explanation of database findings, qualitative annotations, or absence of quantitative assay data. |

## 2. Auditable Absence of Phenotype (`remaining_unphenotyped.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `variant_id` | string | Unique identifier formatted as `{virus}_{gene}_{aa_change}`. |
| `gene` | string | Target viral gene (`UL23` or `UL30`). |
| `ref` | string | Reference amino acid. |
| `pos` | integer | Amino acid position coordinate. |
| `alt` | string | Variant amino acid. |
| `herpesdrg_query` | string | Query expression and outcome for HerpesDRG. |
| `sauerbrei_query` | string | Query and access outcome for Sauerbrei 2016 JAC supplementary database. |
| `herphegen_query` | string | Query and access outcome for herPHEgen 2019 evaluation. |
| `fda_label_query` | string | Inspection outcome of the FDA Acyclovir prescribing information. |
| `pubmed_query` | string | Query expression and outcome for PubMed literature search. |
| `status` | string | Categorical audit outcome (`confirmed_unphenotyped`). |
| `notes` | string | Synthesis of findings explaining why the variant remains unphenotyped. |
