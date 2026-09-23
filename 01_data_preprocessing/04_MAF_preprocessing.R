# ============================================================
# 04_MAF_preprocessing.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Convert a TCGA-CESC MAF into a gene-by-sample binary mutation
#   matrix suitable for multi-omics integration.
#
# Input:
#   TCGA-CESC MAF file.
#
# Output:
#   data/processed/CESC_MAF_processed.rds
#   data/processed/CESC_MAF_gene_frequencies.csv
#
# Important:
#   The binary matrix records whether a gene has at least one
#   retained mutation in a sample. It is not a mutation burden
#   or functional-impact score.
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(2026)

required_packages <- c("data.table", "dplyr")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(data.table)
library(dplyr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

input_file <- Sys.getenv(
  "CESC_MAF_FILE",
  unset = "data/raw/TCGA-CESC.maf"
)

if (!file.exists(input_file)) {
  stop(
    "MAF input file not found: ", input_file,
    "\nSet CESC_MAF_FILE or place the source file at data/raw/TCGA-CESC.maf."
  )
}

maf <- fread(
  input_file,
  data.table = FALSE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
missing_cols <- setdiff(required_cols, names(maf))

if (length(missing_cols) > 0) {
  stop("MAF is missing required columns: ", paste(missing_cols, collapse = ", "))
}

maf <- maf %>%
  filter(
    !is.na(Hugo_Symbol),
    Hugo_Symbol != "",
    !is.na(Tumor_Sample_Barcode),
    Tumor_Sample_Barcode != ""
  )

# Remove obvious non-somatic/empty entries when a Variant_Classification
# column is available.
if ("Variant_Classification" %in% names(maf)) {
  maf <- maf %>%
    filter(
      !is.na(Variant_Classification),
      Variant_Classification != ""
    )
}

maf <- maf %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode)

genes <- sort(unique(maf$Hugo_Symbol))
samples <- sort(unique(maf$Tumor_Sample_Barcode))

mutation_matrix <- matrix(
  0L,
  nrow = length(genes),
  ncol = length(samples),
  dimnames = list(genes, samples)
)

mutation_matrix[
  cbind(
    match(maf$Hugo_Symbol, genes),
    match(maf$Tumor_Sample_Barcode, samples)
  )
] <- 1L

gene_frequencies <- data.frame(
  Gene = rownames(mutation_matrix),
  Mutated_Samples = rowSums(mutation_matrix),
  Frequency = rowMeans(mutation_matrix),
  row.names = NULL,
  check.names = FALSE
) %>%
  arrange(desc(Mutated_Samples), Gene)

saveRDS(
  list(
    mutation_matrix = mutation_matrix,
    gene_frequencies = gene_frequencies,
    maf_pairs = maf
  ),
  "data/processed/CESC_MAF_processed.rds"
)

write.csv(
  gene_frequencies,
  "data/processed/CESC_MAF_gene_frequencies.csv",
  row.names = FALSE
)

message("MAF preprocessing complete.")
message("Genes: ", nrow(mutation_matrix))
message("Samples: ", ncol(mutation_matrix))
