# ============================================================
# 01_RPPA_preprocessing.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Preprocess TCGA-CESC RPPA data for downstream multi-omics
#   harmonisation and keratinization-state analyses.
#
# Expected input:
#   RPPA matrix containing protein features and TCGA sample IDs.
#   The exact source/download is intentionally not hard-coded because
#   TCGA/GDC/MD Anderson source files may change over time.
#
# Output:
#   data/processed/CESC_RPPA_processed.rds
#   data/processed/CESC_RPPA_sample_metadata.csv
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(2026)

required_packages <- c("data.table", "dplyr", "tibble")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(data.table)
library(dplyr)
library(tibble)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

input_file <- Sys.getenv(
  "CESC_RPPA_FILE",
  unset = "data/raw/TCGA-CESC_RPPA.tsv"
)

if (!file.exists(input_file)) {
  stop(
    "RPPA input file not found: ", input_file,
    "\nSet CESC_RPPA_FILE or place the source file at data/raw/TCGA-CESC_RPPA.tsv."
  )
}

rppa_raw <- fread(input_file, data.table = FALSE, check.names = FALSE)

if (nrow(rppa_raw) == 0 || ncol(rppa_raw) < 3) {
  stop("The RPPA input does not contain a usable feature-by-sample matrix.")
}

# Detect a likely feature identifier column.
id_candidates <- c("AGID", "gene_symbol", "protein_name", "peptide_target",
                   "Composite.Element.REF", "Composite.Element.Ref")
id_col <- id_candidates[id_candidates %in% names(rppa_raw)][1]

if (is.na(id_col)) {
  id_col <- names(rppa_raw)[1]
  message("No standard RPPA identifier found; using first column: ", id_col)
}

# Preserve metadata columns where present. Numeric sample columns are retained.
metadata_cols <- intersect(
  c("AGID", "lab_id", "catalog_number", "set_id", "peptide_target",
    "gene_symbol", "protein_name"),
  names(rppa_raw)
)

sample_cols <- setdiff(names(rppa_raw), metadata_cols)

# Keep columns containing at least one numeric value.
is_numeric_like <- vapply(
  rppa_raw[sample_cols],
  function(x) any(!is.na(suppressWarnings(as.numeric(x)))),
  logical(1)
)
sample_cols <- sample_cols[is_numeric_like]

if (length(sample_cols) < 2) {
  stop("Fewer than two numeric RPPA sample columns were detected.")
}

rppa <- rppa_raw[, c(metadata_cols, sample_cols), drop = FALSE]

# Convert sample columns explicitly to numeric.
rppa[sample_cols] <- lapply(
  rppa[sample_cols],
  function(x) suppressWarnings(as.numeric(x))
)

# Remove completely missing protein rows.
keep <- rowSums(!is.na(rppa[, sample_cols, drop = FALSE])) > 0
rppa <- rppa[keep, , drop = FALSE]

# Ensure a stable protein identifier.
if ("peptide_target" %in% names(rppa)) {
  feature_id <- rppa$peptide_target
} else if ("gene_symbol" %in% names(rppa)) {
  feature_id <- rppa$gene_symbol
} else {
  feature_id <- rppa[[id_col]]
}

feature_id <- make.unique(as.character(feature_id))
rownames(rppa) <- feature_id

# Numeric RPPA matrix used by downstream analyses.
rppa_matrix <- as.matrix(rppa[, sample_cols, drop = FALSE])
storage.mode(rppa_matrix) <- "numeric"

# Sample metadata: TCGA barcode is retained exactly as supplied.
sample_metadata <- tibble(
  Sample = colnames(rppa_matrix),
  TCGA_barcode = colnames(rppa_matrix)
)

saveRDS(
  list(
    matrix = rppa_matrix,
    metadata = sample_metadata,
    feature_metadata = rppa[, setdiff(names(rppa), sample_cols), drop = FALSE]
  ),
  "data/processed/CESC_RPPA_processed.rds"
)

write.csv(
  sample_metadata,
  "data/processed/CESC_RPPA_sample_metadata.csv",
  row.names = FALSE
)

message("RPPA preprocessing complete.")
message("Features: ", nrow(rppa_matrix))
message("Samples:  ", ncol(rppa_matrix))
