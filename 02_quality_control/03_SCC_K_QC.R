# ============================================================
# 03_SCC_K_QC.R
# ============================================================
# QC for the SCC-K (keratinizing SCC) training cohort.
# Requires the final study sample annotation:
#   data/processed/CESC_sample_annotations.csv
# with columns: Sample, State
# and State == "SCC_K".
# ============================================================

source("02_quality_control/_subset_QC_helper.R")

run_subset_qc(
  state_name = "SCC_K"
)
