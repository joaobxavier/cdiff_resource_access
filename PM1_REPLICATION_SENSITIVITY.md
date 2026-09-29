# PM1 one-plate sensitivity (2026-09-20)

run_all now runs assess_pm1_equal_replication after rebuild_pm1. It consumes
30 reconstructed plate call tables, reconstructed mouse scores and supplied
chemical-class annotations. It verifies the 700-call union, enumerates all
81 selections and writes results under results/tables/pm1_equal_replication/.
There are nine distinct rank patterns, not 81 independent experiments.
Raw rho ranges 0.67482–0.71830 and adjusted rho 0.43653–0.48038.

The tables retain both the former residual-corr p and the corrected
partial-correlation approximation (AdjustedRankP_ControlDF, 18 df).
The latter is authoritative for inference. The current Figure 4 generator
and expected association table use corrected pooled p = 0.0364940843264458.
Other model tests are unchanged. No primary calls are changed.

This is an incremental code update. Earlier full-package test reports and
release archives describe their dated builds, not a newly completed full run.

Update, 20 September: the complete standard-route rerun now passes all 21
stages and 35 regression checks. These sensitivity tables match the approved
analysis byte for byte. See TEST_REPORT.md for the separate S3B p-value
consistency finding and the distinction between this rerun and the older ZIP.
