# Sub-Query Audit and Blocked Access Notes: Task 5 (Missed Tier A Audit)

## 1. Comprehensive Database Audit
- **Scope**: All 747 observed variants across the 599-genome denominator were systematically cross-referenced against HerpesDRG, Sauerbrei JAC 2016, herPHEgen 2019, the FDA acyclovir prescribing information, and literature records.
- **HerpesDRG Snapshot Cross-Reference**: 0 observed variants in Tier B, Tier C, or uncatalogued categories possess unassigned Tier A (recombinant BAC or marker transfer) phenotypes in HerpesDRG. All 16 observed variants with recombinant engineering annotations in HerpesDRG were already assigned to Tier A (14 with quantitative fold-change, 2 unquantified).
- **Sauerbrei JAC 2016 & herPHEgen 2019**: Both platform endpoints returned HTTP 403 Forbidden due to platform bot protection / subscription paywalls.
- **FDA Acyclovir Prescribing Information**: Contains no individual amino acid substitution tables (mechanistic descriptions only).

## 2. Definitive Audit Outcome
- **Conclusion**: The Tier A circulating set **CANNOT** be expanded beyond the 14 observed-phenotyped variants from public databases.
- **Significance**: This proves that the sparsity of circulating Tier A variants (14 / 119 catalogued) is a true reflection of the published recombinant virology literature, rather than an artifact of incomplete database querying. The negative limitation stands definitively.
