# ============================================================
# QC template for a defined cervical SCC subgroup
# ============================================================
# This script is sourced by the subgroup-specific QC scripts.
# It intentionally does not infer SCC-K/SCC-NK/SCC-NOS labels.
# Those labels must be present in a cohort annotation file.
# ============================================================

run_subset_qc <- function(
    state_name,
    annotation_file = "data/processed/CESC_sample_annotations.csv",
    input_file = "data/processed/CESC_multiomics_harmonised.rds"
) {

  required_packages <- c("ggplot2")
  missing <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) > 0) {
    stop("Install required packages: ", paste(missing, collapse = ", "))
  }

  if (!file.exists(input_file)) {
    stop("Missing harmonised cohort: ", input_file)
  }

  if (!file.exists(annotation_file)) {
    stop(
      "Missing sample annotation file: ", annotation_file,
      "\nThe QC scripts do not infer pathological keratinization labels. ",
      "Provide the study's final sample annotation table."
    )
  }

  obj <- readRDS(input_file)
  annotation <- read.csv(annotation_file, check.names = FALSE)

  if (!all(c("Sample", "State") %in% names(annotation))) {
    stop("Annotation file must contain columns: Sample and State")
  }

  annotation <- annotation[
    !is.na(annotation$Sample) &
      annotation$Sample %in% obj$samples &
      annotation$State == state_name,
    ,
    drop = FALSE
  ]

  if (nrow(annotation) == 0) {
    stop("No samples found for state: ", state_name)
  }

  samples <- annotation$Sample

  get_matrix <- function(x) {
    x[, samples, drop = FALSE]
  }

  matrices <- list(
    RPPA = get_matrix(obj$RPPA),
    RNA = get_matrix(obj$RNA_logCPM),
    Methylation = get_matrix(obj$Methylation),
    MAF = get_matrix(obj$MAF)
  )

  outdir <- file.path("results/QC", state_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  summary <- do.call(
    rbind,
    lapply(names(matrices), function(layer) {
      x <- matrices[[layer]]
      data.frame(
        State = state_name,
        Omics = layer,
        Features = nrow(x),
        Samples = ncol(x),
        Missing_fraction = mean(is.na(x)),
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(
    summary,
    file.path(outdir, paste0(state_name, "_QC_summary.csv")),
    row.names = FALSE
  )

  coverage <- data.frame(
    Sample = samples,
    RPPA = colSums(!is.na(matrices$RPPA)),
    RNA = colSums(!is.na(matrices$RNA)),
    Methylation = colSums(!is.na(matrices$Methylation)),
    MAF_mutated_genes = colSums(matrices$MAF > 0, na.rm = TRUE)
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

  p <- ggplot2::ggplot(coverage_long, ggplot2::aes(x = Value)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(~ Omics, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0(state_name, " multi-omics QC"),
      x = "Per-sample feature count",
      y = "Samples"
    )

  ggplot2::ggsave(
    file.path(outdir, paste0(state_name, "_QC_distributions.png")),
    p, width = 12, height = 8, dpi = 600
  )

  invisible(summary)
}
