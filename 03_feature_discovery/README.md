# 03 — Feature discovery

This stage characterises proteomic consequences associated with three molecular alteration classes in the harmonised cervical cancer cohort:

1. recurrent somatic mutation,
2. gene-associated DNA methylation state, and
3. gene-associated RNA-expression state.

The seven recurrent genes used in the study are **CREBBP, DMD, KMT2D, LRP2, RYR2, SYNE1 and USH2A**. The default analysis population is `SCC_NOS`, matching the feature-discovery analyses described in the study materials. Change `STATE_FILTER` to `NULL` where a full-cohort analysis is required.

## Order of execution

```text
01_RPPA_consequence_analysis.R
02_methylation_consequence_analysis.R
03_RNA_consequence_analysis.R
04_cross_omics_feature_integration.R
05_pathway_TF_analysis.R
```

## Statistical conventions

- RPPA consequence testing uses limma empirical-Bayes linear models.
- Mutation comparisons are mutated versus wild-type.
- Methylation and RNA comparisons use the upper and lower quartiles, with the middle 50% excluded from the binary comparison.
- `gene_specific_FDR` is calculated within each alteration-specific protein test set.
- `global_FDR` is calculated across the complete set of tests produced by that analysis stage.
- The study's principal supported-feature criterion is **FDR < 0.05 and |RPPA difference| >= 0.20**.
- Both gene-specific and global-FDR/effect-size tables are retained so that exploratory within-driver findings are not confused with the more stringent cross-test evidence.

## Required inputs

The scripts expect:

- `data/processed/CESC_multiomics_harmonised.rds`
- `data/processed/CESC_sample_annotations.csv` with `Sample` and `State`
- `data/processed/CESC_gene_methylation_features.csv` with `DriverGene` and `MethylationFeature` for methylation analysis
- `data/processed/RPPA_feature_to_gene_symbol.csv` with `RPPA_feature` and `gene_symbol` for pathway/TF analysis

The last two mappings are intentionally explicit. RPPA and methylation identifiers can differ between data releases, so the repository does not silently guess mappings.

## Interpretation

The feature-discovery stage is association-based. Significant RPPA differences identify proteomic consequences associated with molecular states; they do not establish that a mutation, methylation change or RNA state causally produces the protein change.

The cross-omics stage is designed to expose recurrent and state-specific protein signatures. Pathway and transcription-factor enrichment should be interpreted as enrichment of represented proteins, not direct evidence of pathway activation or transcription-factor activity.
