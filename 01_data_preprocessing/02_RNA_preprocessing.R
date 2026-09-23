# ============================================================
# 02_RNA_preprocessing.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Preprocess TCGA-CESC RNA-seq count data for downstream
#   harmonisation and expression-state analyses.
#
# Input:
#   Raw/processed TCGA-CESC RNA STAR counts matrix.
#
# Output:
#   data/processed/CESC_RNA_processed.rds
#   data/processed/CESC_RNA_sample_metadata.csv
#
# Notes:
#   This script performs feature filtering and TMM normalization.
#   Voom transformation is generated when limma/edgeR are available.
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(2026)

required_packages <- c("data.table", "edgeR", "limma")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(data.table)
library(edgeR)
library(limma)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

input_file <- Sys.getenv(
  "CESC_RNA_FILE",
  unset = "data/raw/TCGA-CESC_RNA_counts.tsv"
)

if (!file.exists(input_file)) {
  stop(
    "RNA input file not found: ", input_file,
    "\nSet CESC_RNA_FILE or place the source file at data/raw/TCGA-CESC_RNA_counts.tsv."
  )
}

rna_raw <- fread(input_file, data.table = FALSE, check.names = FALSE)

if (nrow(rna_raw) == 0 || ncol(rna_raw) < 3) {
  stop("The RNA input does not contain a usable gene-by-sample matrix.")
}

# GDC STAR-count files can contain annotation columns before sample counts.
annotation_candidates <- c(
  "gene_id", "gene_name", "gene_type", "Symbol", "Ensembl_ID",
  "Ensembl", "Gene"
)

annotation_cols <- intersect(annotation_candidates, names(rna_raw))
sample_cols <- setdiff(names(rna_raw), annotation_cols)

is_numeric_like <- vapply(
  rna_raw[sample_cols],
  function(x) any(!is.na(suppressWarnings(as.numeric(x)))),
  logical(1)
)
sample_cols <- sample_cols[is_numeric_like]

if (length(sample_cols) < 2) {
  stop("Fewer than two numeric RNA sample columns were detected.")
}

gene_id_col <- if ("gene_name" %in% names(rna_raw)) {
  "gene_name"
} else if ("Symbol" %in% names(rna_raw)) {
  "Symbol"
} else {
  annotation_cols[1]
}

if (is.na(gene_id_col) || length(gene_id_col) == 0) {
  gene_id <- names(rna_raw)[seq_len(nrow(rna_raw))]
} else {
  gene_id <- as.character(rna_raw[[gene_id_col]])
}

counts <- as.matrix(rna_raw[, sample_cols, drop = FALSE])
storage.mode(counts) <- "numeric"

gene_id <- make.unique(gene_id)
rownames(counts) <- gene_id

# Remove genes with no meaningful counts.
dge <- DGEList(counts = counts)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge <- calcNormFactors(dge, method = "TMM")

# log2 CPM for downstream exploratory/association analyses.
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# limma-voom transformation.
voom_obj <- voom(dge, plot = FALSE)

sample_metadata <- data.frame(
  Sample = colnames(counts),
  TCGA_barcode = colnames(counts),
  LibrarySize = dge$samples$lib.size,
  NormFactor = dge$samples$norm.factors,
  row.names = colnames(counts),
  check.names = FALSE
)

saveRDS(
  list(
    counts = dge$counts,
    dge = dge,
    logCPM = logCPM,
    voom = voom_obj,
    metadata = sample_metadata
  ),
  "data/processed/CESC_RNA_processed.rds"
)

write.csv(
  sample_metadata,
  "data/processed/CESC_RNA_sample_metadata.csv",
  row.names = FALSE
)

message("RNA preprocessing complete.")
message("Retained genes: ", nrow(dge))
message("Samples: ", ncol(dge))
