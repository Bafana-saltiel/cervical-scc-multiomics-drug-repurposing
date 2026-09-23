# ============================================================
# 01_SCC_K_SCC_NK_training.R
# Define the pathology-labelled SCC-K / SCC-NK training cohort
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
})

SEED <- 20260916
set.seed(SEED)

INPUT <- "data/processed/labelled_SCC_RPPA_classifier_data.csv"
ANNOT <- "data/processed/CESC_sample_annotations.csv"
OUTDIR <- "results/classifier"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# The final classifier is trained on pathology-defined keratinizing (SCC-K)
# and nonkeratinizing (SCC-NK) reference tumours. Expected reference size:
# 16 SCC-K + 27 SCC-NK = 43 tumours.

MARKERS <- c("ZAP-70", "EVI1", "EPPK1", "ANNEXIN1", "CD171")

if (!file.exists(INPUT)) {
  stop("Training RPPA table not found: ", INPUT,
       "\nProvide a patient-level table containing Patient, classifier_group and the five markers.")
}

dat <- readr::read_csv(INPUT, show_col_types=FALSE) %>% tibble::as_tibble()
req <- c("Patient", "classifier_group", MARKERS)
miss <- setdiff(req, names(dat))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

training <- dat %>%
  mutate(
    Patient = substr(as.character(Patient),1,12),
    classifier_group = case_when(
      grepl("nonkeratinizing", classifier_group, ignore.case=TRUE) ~ "SCC-NK",
      grepl("keratinizing", classifier_group, ignore.case=TRUE) ~ "SCC-K",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(classifier_group)) %>%
  distinct(Patient, .keep_all=TRUE)

for (m in MARKERS) training[[m]] <- as.numeric(training[[m]])
if (anyNA(training[,MARKERS])) stop("Missing values occur in training markers.")
if (anyDuplicated(training$Patient)) stop("Duplicate patient IDs remain.")

counts <- training %>% count(classifier_group)
print(counts)

if (nrow(training) != 43L) warning("Expected 43 labelled tumours; found ", nrow(training), ".")
if (!all(c("SCC-K","SCC-NK") %in% training$classifier_group)) stop("Both SCC-K and SCC-NK classes are required.")

training$classifier_group <- factor(training$classifier_group, levels=c("SCC-NK","SCC-K"))
training$y <- as.integer(training$classifier_group=="SCC-K")

write_csv(training, file.path(OUTDIR,"SCC_K_SCC_NK_training_cohort.csv"))
write_csv(counts, file.path(OUTDIR,"SCC_K_SCC_NK_training_counts.csv"))

message("Training cohort written to ", OUTDIR)
