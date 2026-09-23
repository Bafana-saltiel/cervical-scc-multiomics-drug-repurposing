# ============================================================
# 06 Drug repurposing: rank L1000 candidates
# ============================================================
# Input: results/drug_repurposing/L1000_query_results.csv
# Expected columns: compound/drug identifier and a reversal score.
# More-negative scores represent stronger predicted reversal in
# the historical analysis.
# ============================================================

suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr) })

IN <- "results/drug_repurposing/L1000_query_results.csv"
OUT <- "results/drug_repurposing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN)) stop("Missing L1000 query result: ", IN)

dat <- read_csv(IN, show_col_types = FALSE)

score_candidates <- c("score", "tau", "connectivity_score", "reversal_score", "mean_score")
score_col <- score_candidates[score_candidates %in% names(dat)][1]
if (is.na(score_col)) stop("No supported L1000 score column found. Expected one of: ", paste(score_candidates, collapse = ", "))

name_candidates <- c("compound", "drug", "perturbagen", "perturbagen_name", "compound_name", "name")
name_col <- name_candidates[name_candidates %in% names(dat)][1]
if (is.na(name_col)) stop("No compound-name column found.")

ranked <- dat %>%
  mutate(reversal_score = as.numeric(.data[[score_col]]),
         compound_name = as.character(.data[[name_col]])) %>%
  filter(is.finite(reversal_score), !is.na(compound_name), compound_name != "") %>%
  arrange(reversal_score) %>%
  mutate(rank = row_number())

write_csv(ranked, file.path(OUT, "L1000_ranked_candidates.csv"))
write_csv(slice_head(ranked, n = 15), file.path(OUT, "L1000_top15_candidates.csv"))

message("Ranked ", nrow(ranked), " compounds; strongest predicted reversals are the most negative scores.")
