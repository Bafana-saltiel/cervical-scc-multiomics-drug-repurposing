# ============================================================
# 06 Drug repurposing: construct keratinizing RNA signature
# ============================================================
# Purpose:
# Construct the transcriptional programme associated with the
# RPPA-defined keratinizing state for downstream L1000 querying.
#
# Required input:
# results/classifier/SCC_NOS/SCC_NOS_classification.csv
# plus an RNA expression matrix and sample metadata.
#
# The historical analysis used 150 upregulated and 150
# downregulated genes. This script calculates those lists from
# the supplied expression data rather than hard-coding genes.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(limma)
})

OUT <- "results/drug_repurposing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

state_file <- "results/classifier/SCC_NOS/SCC_NOS_classification.csv"
rna_file <- Sys.getenv("KERATINIZING_RNA_FILE", "data/processed/CESC_RNA_for_L1000.csv")

if (!file.exists(state_file)) stop("Missing classifier state file: ", state_file)
if (!file.exists(rna_file)) stop("Missing RNA matrix: ", rna_file)

states <- read_csv(state_file, show_col_types = FALSE)
required_state <- c("Patient", "classifier_group")
if (!all(required_state %in% names(states))) {
  stop("Classifier file must contain: ", paste(required_state, collapse = ", "))
}

# Use terminal molecular states only; intermediate SCC-NOS states are
# excluded from construction of the keratinizing-versus-opposing signature.
states <- states %>%
  filter(classifier_group %in% c("Keratinizing", "Nonkeratinizing")) %>%
  select(Patient, classifier_group) %>%
  distinct()

rna <- read_csv(rna_file, show_col_types = FALSE)
if (!"gene" %in% names(rna)) names(rna)[1] <- "gene"

expr <- as.data.frame(rna)
gene <- expr$gene
expr$gene <- NULL
expr <- as.matrix(expr)
rownames(expr) <- gene
storage.mode(expr) <- "numeric"

common <- intersect(states$Patient, colnames(expr))
if (length(common) < 10) stop("Too few matched RNA samples: ", length(common))

states <- states[match(common, states$Patient), ]
expr <- expr[, common, drop = FALSE]

group <- factor(states$classifier_group, levels = c("Nonkeratinizing", "Keratinizing"))
design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr, design))

coef_name <- "groupKeratinizing"
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  rownames_to_column("gene")

up <- tt %>% arrange(desc(logFC)) %>% slice_head(n = 150) %>% mutate(direction = "UP")
down <- tt %>% arrange(logFC) %>% slice_head(n = 150) %>% mutate(direction = "DOWN")
sig <- bind_rows(up, down)

write_csv(sig, file.path(OUT, "keratinizing_state_signature_full.csv"))
write_lines(up$gene, file.path(OUT, "L1000_Keratinizing_UP_symbols.txt"))
write_lines(down$gene, file.path(OUT, "L1000_Keratinizing_DOWN_symbols.txt"))
write_csv(up, file.path(OUT, "L1000_Keratinizing_UP_evidence.csv"))
write_csv(down, file.path(OUT, "L1000_Keratinizing_DOWN_evidence.csv"))

message("Constructed L1000 signature: 150 UP + 150 DOWN genes.")
