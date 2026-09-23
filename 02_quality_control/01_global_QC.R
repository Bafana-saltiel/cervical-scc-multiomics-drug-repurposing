# ============================================================
# 01_global_QC.R
# Cervical SCC Multi-Omics and Drug Repurposing
# ============================================================
# Purpose:
#   Distributional and completeness QC for the harmonised
#   TCGA-CESC multi-omics cohort.
#
# Input:
#   data/processed/CESC_multiomics_harmonised.rds
#
# Output:
#   results/QC/global_QC_summary.csv
#   results/QC/global_QC_distributions.png
#   results/QC/global_QC_missingness.png
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "data.table")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

library(ggplot2)
library(data.table)

input_file <- "data/processed/CESC_multiomics_harmonised.rds"

if (!file.exists(input_file)) {
  stop("Missing harmonised cohort: ", input_file,
       "\nRun 01_data_preprocessing/05_multiomics_cohort_harmonisation.R first.")
}

obj <- readRDS(input_file)

dir.create("results/QC", recursive = TRUE, showWarnings = FALSE)

qc_summary <- data.frame(
  Omics = c("RPPA", "RNA", "Methylation", "MAF"),
  Features = c(
    nrow(obj$RPPA),
    nrow(obj$RNA_logCPM),
    nrow(obj$Methylation),
    nrow(obj$MAF)
  ),
  Samples = c(
    ncol(obj$RPPA),
    ncol(obj$RNA_logCPM),
    ncol(obj$Methylation),
    ncol(obj$MAF)
  ),
  Missing_fraction = c(
    mean(is.na(obj$RPPA)),
    mean(is.na(obj$RNA_logCPM)),
    mean(is.na(obj$Methylation)),
    mean(is.na(obj$MAF))
  ),
  stringsAsFactors = FALSE
)

write.csv(qc_summary, "results/QC/global_QC_summary.csv", row.names = FALSE)

# Sample-level feature coverage.
coverage <- data.frame(
  Sample = obj$samples,
  RPPA = colSums(!is.na(obj$RPPA)),
  RNA = colSums(!is.na(obj$RNA_logCPM)),
  Methylation = colSums(!is.na(obj$Methylation)),
  MAF_mutated_genes = colSums(obj$MAF > 0, na.rm = TRUE)
)

coverage_long <- rbind(
  data.frame(Sample = coverage$Sample, Omics = "RPPA",
             Value = coverage$RPPA),
  data.frame(Sample = coverage$Sample, Omics = "RNA",
             Value = coverage$RNA),
  data.frame(Sample = coverage$Sample, Omics = "Methylation",
             Value = coverage$Methylation),
  data.frame(Sample = coverage$Sample, Omics = "MAF mutated genes",
             Value = coverage$MAF_mutated_genes)
)

p1 <- ggplot(coverage_long, aes(x = Value)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ Omics, scales = "free") +
  theme_bw() +
  labs(
    title = "Global multi-omics QC",
    x = "Per-sample feature count",
    y = "Samples"
  )

ggsave(
  "results/QC/global_QC_distributions.png",
  p1, width = 12, height = 8, dpi = 600
)

missing_by_sample <- data.frame(
  Sample = obj$samples,
  RPPA = colMeans(is.na(obj$RPPA)),
  RNA = colMeans(is.na(obj$RNA_logCPM)),
  Methylation = colMeans(is.na(obj$Methylation)),
  MAF = colMeans(is.na(obj$MAF))
)

missing_long <- reshape(
  missing_by_sample,
  varying = c("RPPA", "RNA", "Methylation", "MAF"),
  v.names = "Missing_fraction",
  timevar = "Omics",
  times = c("RPPA", "RNA", "Methylation", "MAF"),
  direction = "long"
)

p2 <- ggplot(missing_long, aes(x = Omics, y = Missing_fraction)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Global sample-level missingness",
    x = NULL,
    y = "Missing fraction"
  )

ggsave(
  "results/QC/global_QC_missingness.png",
  p2, width = 8, height = 6, dpi = 600
)

message("Global QC complete.")
