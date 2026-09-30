# PM1 equal replication sensitivity

`assess_pm1_equal_replication` runs after `rebuild_pm1`. It reads the 30
reconstructed plate-call tables, mouse disease/protection scores and chemical
classes. The primary matrix combines replicate calls by OR and contains 700
positive calls.

ST1-12, ST1-68, ST1-75 and VPI10463 each have three plates; other strains have
one. Selecting one plate from each replicated strain gives 81 combinations.
The script enumerates every combination without altering the primary matrix
and writes results to `results/tables/pm1_equal_replication/`.

The combinations produce nine distinct rank patterns, not 81 independent
experiments. Raw amino-acid-access/protection correlations range from 0.67482
to 0.71830; disease-adjusted rank correlations range from 0.43653 to 0.48038.
Significance is not retained in every selection.

For inference, `AdjustedRankP_ControlDF` uses the partial-rank approximation
with 18 residual degrees of freedom, accounting for the disease covariate.
The residual-correlation diagnostic is also recorded separately and is not
the inferential p-value. Figure 4E and S2B report the same pooled adjusted
test, p = 0.0364940843264458. Formula and cross-figure equality checks are
included in `test_release`.
