# ============================================================
# 05_multiomics_cohort_harmonisation.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Harmonise RPPA, RNA, DNA methylation and MAF data using shared
#   TCGA sample identifiers.
#
# Expected processed inputs:
#   data/processed/CESC_RPPA_processed.rds
#   data/processed/CESC_RNA_processed.rds
#   data/processed/CESC_methylation_processed.rds
#   data/processed/CESC_MAF_processed.rds
#
# Output:
#   data/processed/CESC_multiomics_harmonised.rds
#   data/processed/CESC_multiomics_sample_manifest.csv
#
# This script does not impose the final diagnosis filter used in
# downstream SCC-state analyses. It establishes the common cohort;
# diagnosis/state filtering belongs to the analysis scripts.
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c("dplyr", "tibble")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(dplyr)
library(tibble)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  RPPA = "data/processed/CESC_RPPA_processed.rds",
  RNA = "data/processed/CESC_RNA_processed.rds",
  Methylation = "data/processed/CESC_methylation_processed.rds",
  MAF = "data/processed/CESC_MAF_processed.rds"
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(
    "Run preprocessing scripts 01-04 first. Missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

rppa <- readRDS(input_files["RPPA"])
rna <- readRDS(input_files["RNA"])
meth <- readRDS(input_files["Methylation"])
maf <- readRDS(input_files["MAF"])

extract_samples <- function(x, object_name) {
  candidates <- list(
    x$metadata$Sample,
    colnames(x$matrix),
    colnames(x$logCPM),
    colnames(x$beta),
    colnames(x$mutation_matrix)
  )
  for (z in candidates) {
    if (!is.null(z) && length(z) > 0) {
      return(as.character(z))
    }
  }
  stop("Could not identify sample identifiers in ", object_name)
}

sample_sets <- list(
  RPPA = extract_samples(rppa, "RPPA"),
  RNA = extract_samples(rna, "RNA"),
  Methylation = extract_samples(meth, "methylation"),
  MAF = extract_samples(maf, "MAF")
)

common_samples <- Reduce(intersect, sample_sets)

if (length(common_samples) < 2) {
  stop(
    "Fewer than two samples are shared across all four omics layers. ",
    "Check TCGA barcode formatting before harmonisation."
  )
}

# Preserve the original order of each layer while selecting the shared cohort.
rppa_idx <- match(common_samples, colnames(rppa$matrix))
rna_idx <- match(common_samples, colnames(rna$counts))
meth_idx <- match(common_samples, colnames(meth$beta))
maf_idx <- match(common_samples, colnames(maf$mutation_matrix))

rppa_matrix <- rppa$matrix[, rppa_idx, drop = FALSE]
rna_counts <- rna$counts[, rna_idx, drop = FALSE]
rna_logCPM <- rna$logCPM[, rna_idx, drop = FALSE]
rna_voom <- rna$voom$E[, rna_idx, drop = FALSE]
meth_beta <- meth$beta[, meth_idx, drop = FALSE]
maf_matrix <- maf$mutation_matrix[, maf_idx, drop = FALSE]

colnames(rppa_matrix) <- common_samples
colnames(rna_counts) <- common_samples
colnames(rna_logCPM) <- common_samples
colnames(rna_voom) <- common_samples
colnames(meth_beta) <- common_samples
colnames(maf_matrix) <- common_samples

sample_manifest <- tibble(
  Sample = common_samples,
  RPPA = common_samples,
  RNA = common_samples,
  Methylation = common_samples,
  MAF = common_samples
)

saveRDS(
  list(
    samples = common_samples,
    RPPA = rppa_matrix,
    RNA_counts = rna_counts,
    RNA_logCPM = rna_logCPM,
    RNA_voom = rna_voom,
    Methylation = meth_beta,
    MAF = maf_matrix,
    sample_manifest = sample_manifest
  ),
  "data/processed/CESC_multiomics_harmonised.rds"
)

write.csv(
  sample_manifest,
  "data/processed/CESC_multiomics_sample_manifest.csv",
  row.names = FALSE
)

message("Multi-omics cohort harmonisation complete.")
message("Common samples: ", length(common_samples))
message("RPPA features: ", nrow(rppa_matrix))
message("RNA features: ", nrow(rna_logCPM))
message("Methylation probes: ", nrow(meth_beta))
message("MAF genes: ", nrow(maf_matrix))
