# Drug repurposing

This stage links the RPPA-defined cervical SCC keratinization state to perturbational-signature matching.

## Workflow

1. `01_keratinizing_signature.R` constructs the keratinizing RNA programme from terminal keratinizing versus nonkeratinizing molecular states and exports the 150-UP/150-DOWN signature used in the study.
2. `02_L1000_query.R` prepares the signature for L1000/CMap querying and records the query manifest.
3. `03_L1000_ranked_candidates.R` standardises and ranks returned compound-level results by reversal score.
4. `04_candidate_summary.R` creates the manuscript-facing top-candidate tables.

## Interpretation

The L1000 step is a computational signature-reversal analysis. More-negative reversal scores indicate stronger predicted opposition to the keratinizing transcriptional programme. They do **not** establish drug sensitivity, mechanism, clinical benefit, or therapeutic efficacy.

The repository deliberately does not fabricate external L1000 query results. Returned results should be placed at:

`results/drug_repurposing/L1000_query_results.csv`

The exact column names can vary between L1000/CMap interfaces; script 03 accepts common score and compound-name fields and otherwise stops with an explicit error.
