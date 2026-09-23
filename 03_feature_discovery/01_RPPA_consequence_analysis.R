# ============================================================
# 01_RPPA_consequence_analysis.R
# Proteomic consequences of recurrent somatic alterations
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ------------------------- CONFIG ----------------------------
INPUT_RDS <- "data/processed/CESC_multiomics_harmonised.rds"
ANNOTATION_FILE <- "data/processed/CESC_sample_annotations.csv"
OUTDIR <- "results/feature_discovery/RPPA"
DRIVER_GENES <- c("CREBBP", "DMD", "KMT2D", "LRP2", "RYR2", "SYNE1", "USH2A")
STATE_FILTER <- "SCC_NOS"   # set NULL to analyse the full harmonised cohort
MIN_MUTATED <- 10
MIN_WT <- 10
FDR_CUTOFF <- 0.05
EFFECT_CUTOFF <- 0.20
# -------------------------------------------------------------

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(INPUT_RDS)
ann <- fread(ANNOTATION_FILE, data.table = FALSE)
stopifnot(all(c("Sample", "State") %in% names(ann)))

get_component <- function(x, candidates) {
  if (is.list(x)) {
    nms <- names(x)
    hit <- candidates[candidates %in% nms]
    if (length(hit)) return(x[[hit[1]]])
  }
  NULL
}

rppa <- get_component(obj, c("RPPA", "rppa", "RPPA_matrix", "rppa_matrix"))
maf  <- get_component(obj, c("MAF", "maf", "MAF_matrix", "maf_matrix"))
if (is.null(rppa) || is.null(maf)) stop("Harmonised RDS must contain RPPA and MAF components.")

rppa <- as.matrix(rppa)
if (is.null(colnames(rppa))) stop("RPPA matrix must have sample IDs as column names.")
if (!is.numeric(rppa)) mode(rppa) <- "numeric"
maf <- as.matrix(maf)
if (is.null(rownames(maf)) || is.null(colnames(maf))) stop("MAF matrix must have gene and sample names.")

samples <- intersect(colnames(rppa), colnames(maf))
if (!is.null(STATE_FILTER)) samples <- intersect(samples, ann$Sample[ann$State == STATE_FILTER])
if (!length(samples)) stop("No samples remain after sample matching/state filtering.")
rppa <- rppa[, samples, drop = FALSE]
maf <- maf[, samples, drop = FALSE]

# Collapse duplicated RPPA feature labels if necessary.
rownames(rppa) <- make.unique(rownames(rppa))

all_results <- list()
for (gene in DRIVER_GENES) {
  if (!gene %in% rownames(maf)) next
  y <- as.numeric(maf[gene, ] > 0)
  names(y) <- colnames(maf)
  n_mut <- sum(y == 1, na.rm = TRUE)
  n_wt <- sum(y == 0, na.rm = TRUE)
  if (n_mut < MIN_MUTATED || n_wt < MIN_WT) next

  design <- model.matrix(~ factor(y, levels = c(0, 1)))
  colnames(design) <- c("Intercept", "Mutated")
  fit <- eBayes(lmFit(rppa, design))
  tt <- topTable(fit, coef = "Mutated", number = Inf, sort.by = "none")
  tt$protein <- rownames(tt)
  tt$driver_gene <- gene
  tt$n_mutated <- n_mut
  tt$n_wildtype <- n_wt
  tt$RPPA_difference <- tt$logFC
  tt$gene_specific_FDR <- tt$adj.P.Val
  all_results[[gene]] <- tt[, c("driver_gene","protein","n_mutated","n_wildtype","RPPA_difference","P.Value","gene_specific_FDR")]
}

res <- bind_rows(all_results)
if (!nrow(res)) stop("No driver gene met the minimum group-size requirements.")
res$global_FDR <- p.adjust(res$P.Value, method = "BH")
res$FDR_effect_supported <- res$gene_specific_FDR < FDR_CUTOFF & abs(res$RPPA_difference) >= EFFECT_CUTOFF
res$global_FDR_effect_supported <- res$global_FDR < FDR_CUTOFF & abs(res$RPPA_difference) >= EFFECT_CUTOFF

fwrite(res, file.path(OUTDIR, "RPPA_mutation_consequences_all.csv"))
fwrite(filter(res, FDR_effect_supported), file.path(OUTDIR, "RPPA_mutation_consequences_gene_FDR_effect.csv"))
fwrite(filter(res, global_FDR_effect_supported), file.path(OUTDIR, "RPPA_mutation_consequences_global_FDR_effect.csv"))

summary_tbl <- res %>% group_by(driver_gene) %>% summarise(
  n_mutated = first(n_mutated), n_wildtype = first(n_wildtype),
  nominal = sum(P.Value < 0.05),
  gene_FDR = sum(gene_specific_FDR < FDR_CUTOFF),
  gene_FDR_effect = sum(FDR_effect_supported),
  global_FDR = sum(global_FDR < FDR_CUTOFF),
  global_FDR_effect = sum(global_FDR_effect_supported), .groups = "drop")
fwrite(summary_tbl, file.path(OUTDIR, "RPPA_mutation_consequence_summary.csv"))

pdat <- res %>% group_by(driver_gene) %>% slice_min(global_FDR, n = 10, with_ties = FALSE) %>% ungroup()
p <- ggplot(pdat, aes(x = -log10(global_FDR), y = reorder(protein, -global_FDR))) +
  geom_point(aes(size = abs(RPPA_difference), shape = global_FDR_effect_supported)) +
  facet_wrap(~ driver_gene, scales = "free_y") +
  labs(x = "-log10(global FDR)", y = "RPPA feature", size = "|RPPA difference|", shape = "Global FDR + effect") +
  theme_bw()
ggsave(file.path(OUTDIR, "RPPA_mutation_consequences_top10.png"), p, width = 12, height = 9, dpi = 600)

message("Completed RPPA mutation-consequence analysis: ", nrow(res), " tests.")
