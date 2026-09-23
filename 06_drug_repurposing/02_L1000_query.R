# ============================================================
# 06 Drug repurposing: L1000 query
# ============================================================
# Purpose:
# Prepare the keratinizing-state signature for L1000/CMap-style
# perturbational-signature matching.
#
# This repository does not silently fabricate API credentials,
# URLs, or query results. The script supports a local result file
# supplied after querying the L1000 resource through the permitted
# access route.
# ============================================================

suppressPackageStartupMessages({ library(readr); library(dplyr) })

IN <- "results/drug_repurposing"
OUT <- "results/drug_repurposing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

up_file <- file.path(IN, "L1000_Keratinizing_UP_symbols.txt")
down_file <- file.path(IN, "L1000_Keratinizing_DOWN_symbols.txt")
if (!file.exists(up_file) || !file.exists(down_file)) {
  stop("Run 01_keratinizing_signature.R first.")
}

up <- read_lines(up_file)
down <- read_lines(down_file)

write_lines(up, file.path(OUT, "@L1000_Keratinizing_UP_symbols.txt"))
write_lines(down, file.path(OUT, "@L1000_Keratinizing_DOWN_symbols.txt"))

manifest <- tibble(
  signature = c("keratinizing_UP", "keratinizing_DOWN"),
  n_genes = c(length(up), length(down)),
  file = c(basename(up_file), basename(down_file)),
  query_direction = c("up", "down")
)
write_csv(manifest, file.path(OUT, "L1000_query_manifest.csv"))

message("Signature prepared. Submit UP/DOWN lists to the approved L1000/CMap query service and save the returned compound-level result as:")
message(file.path(OUT, "L1000_query_results.csv"))
