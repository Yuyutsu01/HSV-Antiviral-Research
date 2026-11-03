# Data Dictionary: Task 5 Audit for Missed Phenotyped Circulating Variants

Table file: `newly_phenotyped_observed.tsv`

| Column | Type | Description |
| :--- | :--- | :--- |
| `variant_id` | string | Unique variant identifier formatted as `{virus}_{gene}_{aa_change}`. |
| `gene` | string | Target viral gene (`UL23`, `UL30`, `UL5`, `UL52`). |
| `fold_change` | string/float | Quantitative in vitro susceptibility fold-change (IC50 or EC50 ratio vs wild-type). |
| `source` | string | Database, literature DOI, or accession origin of the newly discovered phenotype. |
| `why_missed` | string | Detailed audit explanation of why the phenotype was omitted from the initial Tier A release. |
