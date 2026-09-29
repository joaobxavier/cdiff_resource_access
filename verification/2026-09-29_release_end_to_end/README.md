# Complete standard-route release test

Date: 29 September 2026. Tested archive SHA-256:
`a95171d5d043bc3d8fc335109701bac1db016f5f04910d15c4a84975baac676c`.

The ZIP was extracted to
`an isolated temporary package copy`.
The initial results directory did not exist. No previous experimental outputs
were supplied. MATLAB used only this copy plus installed dependencies and
the preserved MCP infrastructure. `run_all('standard')` completed all 21
stages, totaling 393.397 seconds. This route uses the distributed derived-tree
cache; no new genome reconstruction was performed.

- `run_status.csv`: all 21 stages PASS.
- `regression_tests.csv`: 35 passing checks, 11,026 numerical values.
- `static_analysis.csv`: 19 nonfatal diagnostics across 48 MATLAB files;
  no parser errors.
- `bootstrap_summary.csv`: 1,999 null-refitted replicates, p = 0.013.
- `tree_sensitivity_comparison.csv`: reported sensitivity conclusions checked.
- `figure1_numerical_regression.csv`: three extra checks, 266 values.
- `figure1_ks65_mouse_identity.csv`: source-verified identity mapping.

All 158 input hashes passed. Exactly eight active figures and their editable
exports were generated; their panel composition and labels were visually
inspected. The standard call exceeded the MCP request's 300-second timeout,
but MATLAB continued and wrote completion evidence. A later MCP request
confirmed the 21 PASS records and reran regression/identity checks successfully.

The final release differs only by license/documentation/verification files;
its MATLAB code, scientific inputs, expectations and tree cache were compared
byte for byte with this tested copy. No repository has been published.
