# ============================================================
# 04_SCC_NK_QC.R
# ============================================================
# QC for the SCC-NK (nonkeratinizing SCC) training cohort.
# Requires the final study sample annotation:
#   data/processed/CESC_sample_annotations.csv
# with columns: Sample, State
# and State == "SCC_NK".
# ============================================================

source("02_quality_control/_subset_QC_helper.R")

run_subset_qc(
  state_name = "SCC_NK"
)
