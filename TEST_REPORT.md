# Validation report — development snapshot 1.1.0-rc.1

## Standard route

A fresh package copy, containing no experimental result outputs, completed
all 22 standard-route stages on 1 October 2026 with MATLAB R2024b on macOS ARM64.
The supplied phylogenetic cache was used; this was not a new genome reconstruction.
The complete numerical and identity suites passed afterward. All executable
code in the tested copy matches this snapshot. Final comment-only cleanups
were checked for executable-line identity; static checks were refreshed.

- All 161 manifest-listed input/control hashes passed.
- The original 35 regression checks passed for all 11,026 numerical values.
- All nine new rank associations matched MATLAB's native routines to 1e-12.
- All eight auxiliary class-adjusted tests matched MATLAB `partialcorr`,
  including one-covariate degrees of freedom and the eight-class BH correction.
- New frozen checks cover coefficients, p/q values, intervals, degrees of
  freedom, bootstrap counts, and exact strain/input matching.
- All 189 strain-deletion results and all 81 plate selections were retained.
- All 10,521 model-grid outcomes matched the former parameterization exactly.
- KS65 mouse identities and Figure 1 numerical checks passed.
- All eight figure exports were present with editable PDF, SVG and FIG files.
- Static analysis covered 54 MATLAB files with no parser errors and 20
  retained nonfatal diagnostics (style, performance or unused-code messages).

`verification/` contains stage status, regression results, static diagnostics,
identity checks and ecological-analysis checks. These records are not calculation
inputs. Sixteen auxiliary class-adjusted p/q entries in the original expectations
were corrected after comparison with native MATLAB tests (18 rather than 19 df).
No coefficients, raw class tests, displayed statistics or conclusions changed.
All other original expectations are unchanged; overlap expectations are in
`expected/resource_overlap/`. Figure 4 and S2 were visually inspected, and their
rendered pixels match the manuscript artwork exactly. Saved-figure checks verify
the current panel composition and retained-event cytometry labels.

## Interpretation and limits

The standard route regenerates experimental analyses and downstream tree tests
from supplied inputs. It does not repeat raw-read assembly, GC-MS peak
integration, FCS gating or full genome-tree reconstruction. The prior full-route
tree-reference checks remain available in verification/genome_reconstruction_reference.csv.

Bootstrap intervals resample 21 strain-level point estimates, not individual
mice. New overlap tests do not account for phylogenetic relatedness; the
separate existing phylogeny-adjusted relative-amino-acid analysis is retained.
One constant bootstrap draw was invalid for amino-acid overlap (4,999 valid).
No invalid one-plate selections occurred. Neither the overlap tests nor the
existing associations establish a validated out-of-sample prediction rule.

This expanded package is a **development snapshot**, not a new tagged release.
Identify it by its Git commit. Public version 1.0.0 remains unchanged; no release
tag has been created or moved.
