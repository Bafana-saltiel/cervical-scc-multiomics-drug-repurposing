# ============================================================
# 02_feature_selection.R
# Prespecify and document the final five-protein RPPA panel
# ============================================================

suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble)})

INPUT <- "results/classifier/SCC_K_SCC_NK_training_cohort.csv"
OUTDIR <- "results/classifier"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

FINAL_MARKERS <- c("ZAP-70","EVI1","EPPK1","ANNEXIN1","CD171")

if (!file.exists(INPUT)) stop("Run 01_SCC_K_SCC_NK_training.R first: ", INPUT)
dat <- readr::read_csv(INPUT, show_col_types=FALSE)
miss <- setdiff(c("Patient","classifier_group",FINAL_MARKERS),names(dat))
if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))

# This file documents the final study panel rather than silently replacing it
# with a new data-driven selection. Direction and approximate fitted weights
# are reported by the final model in script 03.
manifest <- tibble(
  feature = FINAL_MARKERS,
  direction = c("Nonkeratinizing","Nonkeratinizing","Keratinizing","Keratinizing","Keratinizing"),
  role = "Final five-marker RPPA classifier"
)

write_csv(manifest,file.path(OUTDIR,"five_marker_feature_manifest.csv"))

summary <- dat %>%
  summarise(across(all_of(FINAL_MARKERS), list(mean=~mean(.x), sd=~sd(.x)), .names="{.col}_{.fn}"))
write_csv(summary,file.path(OUTDIR,"five_marker_training_summary.csv"))

message("Final prespecified panel: ",paste(FINAL_MARKERS,collapse=", "))
