# 05 — Validation

This stage evaluates the fixed five-protein RPPA keratinization classifier against orthogonal RNA measurements and independent SCC cohorts.

**Final panel:** ZAP-70, EVI1, EPPK1, ANNEXIN1 (ANXA1), CD171.

- `01_RPPA_vs_RNA_classifier.R`: matched pathology-defined SCC-K/SCC-NK comparison of RPPA, RNA expression of the five corresponding genes, and KRT10/KRT13/KRT16/IVL.
- `02_DeLong_comparison.R`: paired DeLong tests on identical patients.
- `03_TCGA_validation.R`: fixed-model application to supplied TCGA-HNSC, TCGA-LUSC, TCGA-ESCA and TCGA-CESC RPPA matrices. No refitting.
- `04_external_SCC_validation.R`: fixed-model application to a supplied independent SCC RPPA cohort.

Historical manuscript values are not hard-coded: the direct comparison reported AUC 0.877 for the five-protein RPPA classifier, 0.704 for RNA expression of the corresponding five genes, and 0.572 for the conventional four-gene panel; paired DeLong RPPA versus conventional RNA gave ΔAUC 0.306, P=0.0015. The scripts recompute analyses from supplied inputs.

External cohorts require appropriately matched feature names; pathology labels are optional and are only used for descriptive validation, never for model fitting.
