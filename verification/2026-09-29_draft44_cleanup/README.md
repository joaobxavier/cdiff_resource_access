# Draft 44 affected-stage verification

An isolated package copy regenerated Figure 1, Figure 4/S2 and Figure 5/S3
using verified upstream outputs from the 20 September corrected test. The
complete regression suite passed again; this is not a new uninterrupted
21-stage analysis run or genome reconstruction.

Included records:

- `regression_tests.csv`: 35 passing checks, 11,026 numerical values.
- `static_analysis.csv`: diagnostics from checking all 48 MATLAB files;
  no parser errors, 19 nonfatal diagnostics.
- `tree_sensitivity_comparison.csv`: retained side-by-side tree-sensitivity
  estimates, with reported conclusions checked by the suite.
- `figure1_ks65_mouse_identity.csv`: source-verified KS65 mouse mapping.
- `figure1_numerical_regression.csv`: three passing additional Figure 1
  checks, covering 266 numerical values.

All 158 original input hashes passed. The suite also confirmed exactly eight
active figures with editable exports, corrected flow-axis labels and the
two-panel S3. Repeating the retirement helper with the current figure set
preserved the active exports and the regression pass. Figure 4, S2, Figure 5
and S3 were visually inspected. Numerical expectations are unchanged.

See the release's TEST_REPORT.md for historical full-run evidence and scope.
