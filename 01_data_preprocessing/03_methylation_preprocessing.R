# ============================================================
# 03_methylation_preprocessing.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Prepare TCGA-CESC DNA methylation data for harmonised
#   multi-omics analysis.
#
# Input:
#   Beta-value matrix with CpGs/Illumina probe IDs as rows and
#   TCGA samples as columns.
#
# Output:
#   data/processed/CESC_methylation_processed.rds
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(2026)

required_packages <- c("data.table")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(data.table)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

input_file <- Sys.getenv(
  "CESC_METHYLATION_FILE",
  unset = "data/raw/TCGA-CESC_methylation_beta.tsv"
)

if (!file.exists(input_file)) {
  stop(
    "Methylation input file not found: ", input_file,
    "\nSet CESC_METHYLATION_FILE or place the source file at ",
    "data/raw/TCGA-CESC_methylation_beta.tsv."
  )
}

meth_raw <- fread(input_file, data.table = FALSE, check.names = FALSE)

if (nrow(meth_raw) == 0 || ncol(meth_raw) < 3) {
  stop("The methylation input does not contain a usable CpG-by-sample matrix.")
}

# First column is treated as the CpG/probe identifier unless a standard
# identifier column is available.
id_candidates <- c("IlmnID", "Probe_ID", "probe_id", "CpG", "cg_id")
id_col <- id_candidates[id_candidates %in% names(meth_raw)][1]

if (is.na(id_col)) {
  id_col <- names(meth_raw)[1]
}

sample_cols <- setdiff(names(meth_raw), id_col)

is_numeric_like <- vapply(
  meth_raw[sample_cols],
  function(x) any(!is.na(suppressWarnings(as.numeric(x)))),
  logical(1)
)
sample_cols <- sample_cols[is_numeric_like]

if (length(sample_cols) < 2) {
  stop("Fewer than two numeric methylation sample columns were detected.")
}

probe_id <- make.unique(as.character(meth_raw[[id_col]]))

beta <- as.matrix(meth_raw[, sample_cols, drop = FALSE])
storage.mode(beta) <- "numeric"
rownames(beta) <- probe_id

# Keep probes with sufficient observed beta values.
observed_fraction <- rowMeans(!is.na(beta))
beta <- beta[observed_fraction >= 0.80, , drop = FALSE]

# Guard against accidental non-beta values.
range_ok <- apply(
  beta,
  1,
  function(x) all(is.na(x) | (x >= 0 & x <= 1))
)

if (any(!range_ok)) {
  warning(sum(!range_ok), " rows contain values outside [0,1] and were removed.")
  beta <- beta[range_ok, , drop = FALSE]
}

# Remove duplicate probe identifiers if the input contains them.
beta <- beta[!duplicated(rownames(beta)), , drop = FALSE]

sample_metadata <- data.frame(
  Sample = colnames(beta),
  TCGA_barcode = colnames(beta),
  row.names = colnames(beta),
  check.names = FALSE
)

saveRDS(
  list(
    beta = beta,
    metadata = sample_metadata
  ),
  "data/processed/CESC_methylation_processed.rds"
)

message("Methylation preprocessing complete.")
message("Retained CpGs: ", nrow(beta))
message("Samples: ", ncol(beta))
