# Sub-Query Audit and Blocked Access Notes: Task 1 (HSV-2 UL23)

## 1. NCBI PopSet E-Utilities Endpoint
- **Query / Endpoint**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=popset`
- **Result / Failure**: NCBI E-utilities returned HTTP error / JSON error: `Invalid db name specified: popset`.
- **Root Cause**: NCBI retired `popset` as an independent queryable Entrez database in the public E-utilities schema.
- **Resolution**: PopSet web interface (`https://www.ncbi.nlm.nih.gov/popset/?term=Herpes+simplex+virus+2+UL23`) was queried directly. Found PopSet cluster UID `1132402808` containing HSV-2 UL23 sequence `KX683521.1`, which was successfully harvested.

## 2. ENA Project PRJEB64252
- **Query / Endpoint**: `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB64252&result=read_run`
- **Finding**: PRJEB64252 consists of 66 raw sequencing read runs (`ERR11681237`–`ERR11681302`). ENA contains 0 assembled analysis/sequence records under this project accession because the authors (Giorgi et al., 2023, Antiviral Res) deposited the assembled consensus coding sequences in GenBank under accession series `OR232861`–`OR232890`.
- **Resolution**: The assembled HSV-2 UL23 consensus sequences (`OR232889` and `OR232890`) were retrieved and analyzed via the GenBank query path (1e). Raw read re-assembly is outside the scope of consensus harvesting.

## 3. DDBJ API Endpoint Transition
- **Query / Endpoint**: `http://api.ddbj.nig.ac.jp/search`
- **Result / Failure**: Connection dropped / TLS unrecognized name on legacy endpoints.
- **Root Cause**: DDBJ transitioned their search architecture to `https://ddbj.nig.ac.jp/search/api/entries/`.
- **Resolution**: Queried the modern REST API endpoint. DDBJ index returned matching BioProject/SRA entries (e.g. `ERP163787`), whose corresponding INSDC sequence accessions were cross-checked.
