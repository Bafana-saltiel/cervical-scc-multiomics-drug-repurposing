# ============================================================
# 02_SCC_NOS_QC.R
# ============================================================
# QC for the SCC-NOS cohort used for subsequent state resolution.
# Requires the final study sample annotation:
#   data/processed/CESC_sample_annotations.csv
# with columns: Sample, State
# and State == "SCC_NOS".
# ============================================================

source("02_quality_control/_subset_QC_helper.R")

run_subset_qc(
  state_name = "SCC_NOS"
)
