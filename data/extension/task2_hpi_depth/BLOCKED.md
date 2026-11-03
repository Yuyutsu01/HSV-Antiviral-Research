# Sub-Query Audit and Blocked Access Notes: Task 2 (UL5 / UL52 HPI Depth)

## 1. ENA Portal API Sequence Length Query Field
- **Query / Endpoint**: `https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=tax_tree(10298)%20AND%20sequence_length%3E=140000`
- **Result / Failure**: HTTP 500 error from ENA portal API.
- **Root Cause**: The ENA portal API does not index nucleotide length under the parameter `sequence_length` in the `sequence` table; the indexed column is `base_count`.
- **Resolution**: Re-queried using `base_count>=140000 AND base_count<=160000` (`tax_tree(10298)`) and `base_count>=150000 AND base_count<=170000` (`tax_tree(10310)`), which successfully returned 528 HSV-1 and 513 HSV-2 records.

## 2. Global Public Availability Constraint for HSV-2 Full Genomes
- **Finding**: While HSV-1 has 594 full genomes in NCBI yielding 563 QC-passing UL5 and 564 QC-passing UL52 sequences (total natural depth approaching ~560), HSV-2 has only 514 complete genomes in public databases worldwide.
- **Impact on Success Criterion**: The criterion requested >= 300 additional UL5 and >= 300 additional UL52 sequences per virus. For HSV-2, extracting all available non-duplicate public genomes yields 76 additional natural UL5 sequences (total 318) and 63 additional natural UL52 sequences (total 238). Public database deposition is currently insufficient to supply 300 additional full-length genomes for HSV-2. The negative limitation stands for HSV-2 until new population-scale sequencing projects are released.
