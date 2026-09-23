# ============================================================
# 06 Drug repurposing: candidate summary
# ============================================================
# Produces the manuscript-facing candidate table from ranked
# L1000 results. No therapeutic efficacy is inferred here.
# ============================================================

suppressPackageStartupMessages({ library(readr); library(dplyr) })

IN <- "results/drug_repurposing/L1000_ranked_candidates.csv"
OUT <- "results/drug_repurposing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN)) stop("Run 03_L1000_ranked_candidates.R first.")

dat <- read_csv(IN, show_col_types = FALSE)

summary <- dat %>%
  select(rank, compound_name, reversal_score, everything()) %>%
  slice_head(n = 15)

write_csv(summary, file.path(OUT, "L1000_candidate_summary_top15.csv"))

# A compact figure-ready table.
write_csv(summary %>% select(rank, compound_name, reversal_score),
          file.path(OUT, "L1000_candidate_summary_figure.csv"))

message("Candidate summary written. Negative reversal scores indicate predicted transcriptional opposition to the keratinizing RNA programme; they do not establish drug efficacy.")
