# Sub-Query Audit and Blocked Access Notes: Task 3 (Phenotype Cross-Reference)

## 1. Sauerbrei Database (JAC 2016)
- **Target URL**: `https://academic.oup.com/jac/article/71/1/6/2363994` (DOI: 10.1093/jac/dkv285)
- **Result / Failure**: HTTP 403 Forbidden.
- **Root Cause**: Oxford Academic web platform enforces Cloudflare anti-bot bot-protection and subscription paywalling, preventing automated programmatic retrieval of supplementary tables S1–S8.
- **Impact**: Supplementary XLSX cannot be crawled automatically.

## 2. herPHEgen Evaluation Paper (2019)
- **Target URL**: `https://www.microbiologyresearch.org/content/journal/acmi/10.1099/acmi.ac2019.po0461`
- **Result / Failure**: HTTP 403 Forbidden.
- **Root Cause**: Microbiology Society journal platform returned HTTP 403 to automated HTTP requests.
- **Impact**: Supplementary data package could not be downloaded via script.

## 3. FDA Acyclovir Label Resistance Table
- **Target URL**: `https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo=018604`
- **Finding**: FDA package insert (Section 12.4 Microbiology: Antiviral Resistance) describes mechanisms conceptually (TK-deficient, TK-altered, and DNA polymerase mutants), but does NOT contain a table of individual amino acid substitutions.
- **Impact**: Zero individual amino acid variants are catalogued in the FDA acyclovir label.

## 4. HerpesDRG and PubMed Qualitative vs. Quantitative Gap
- **Finding**: While all 41 variants are catalogued in HerpesDRG, they are annotated qualitatively as "Polymorphism" (with 2 entries having unverified/conflicting notes, e.g. A192V and E678G), with 0 variants possessing a quantitative fold-change measurement.
- **Audit Outcome**: All 41 variants are explicitly confirmed to lack isolated numeric fold-change measurements across all 5 channels, fulfilling the task requirement of unambiguous auditability.
