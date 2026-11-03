# Data Dictionary: Task 4 Independent Carriers and Phenotypes

## 1. Carrier Tables (`r222c_carriers.tsv`, `r176q_carriers.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `accession` | string | Nucleotide sequence accession identifier (GenBank / INSDC). |
| `study_doi` | string | BioProject, PubMed ID, or primary study DOI identifying the study origin. |
| `country` | string | Geographic origin annotated in the record. Empty if unannotated. |
| `year` | string | Year of isolate collection. Empty if unannotated. |
| `already_in_599` | boolean | `True` if accession is present in the manuscript's 599-genome set; `False` if independent/new. |
| `source_url` | string | Web URL linking directly to the NCBI sequence record. |

## 2. Phenotype Table (`r41h_phenotypes.tsv`)

| Column | Type | Description |
| :--- | :--- | :--- |
| `study_doi` | string | DOI of the primary publication reporting the phenotypic assay. |
| `fold_change` | string | Quantitative in vitro susceptibility fold-change (IC50 ratio vs wild-type). |
| `assay` | string | Phenotypic susceptibility assay method (e.g. `Plaque reduction assay (PRA)`, `Dye-uptake assay`). |
| `cell_line` | string | Host cell line utilized in the phenotypic assay (e.g. `Vero`). |
| `drug` | string | Antiviral compound tested (`Aciclovir`). |
| `source_url` | string | Primary literature link to the publication. |
