# Release validation

Package version 1.0.0 was verified on 29 September 2026 in a fresh standalone
copy with no generated experimental results. MATLAB R2024b Update 6
(24.2.0.2923080) on macOS ARM64 was used with Statistics and Machine Learning
Toolbox and Bioinformatics Toolbox.

## Standard analysis

`run_all` completed all 21 stages, totaling 332.241 seconds of recorded stage
time. It recalculated the experimental analyses, downstream tree/trait tests,
models and figures using the explicitly supplied derived-tree cache.
`run_all('tests')` and `test_figure1_identity` also passed when rerun after
completion.

- All 158 manifest-listed input/control files passed SHA-256 verification.
- All 35 numerical regression checks passed, covering 11,026 finite values.
- Figure 1 verification passed three additional comparisons covering 266
  values; the largest absolute difference was 9.9476 × 10⁻¹⁴.
- The disease-adjusted partial-rank formula and Figure 4E/S2B agreement checks
  passed. The synthetic tied-alignment parser fixture also passed.
- All eight main/supplementary PNG figures and their PDF, SVG and MATLAB FIG
  exports were generated and visually inspected. Content checks confirmed
  the retained-event cytometry labels and two-panel microbiota supplement.
- The 1,999-replicate phylogeny/disease-adjusted bootstrap reproduced p = 0.013.
- Static analysis covered 47 MATLAB files with no parser errors. Nineteen
  nonfatal style/performance diagnostics are retained in the supplied report.

Records are in `verification/`: stage status, numerical regressions, Figure 1
identity checks, static analysis, bootstrap summary and tree-sensitivity
comparisons. Expected numerical tables are verification references, never
inputs for fitting the study models.

## Genome reconstruction

The standard run above is not a fresh genome-to-tree reconstruction.
`run_all('full')` performs that separate route with the external dependencies
specified in `GENOMICS.md`. Genome reconstruction was independently tested on
19 September 2026 using the supplied assemblies and pinned macOS ARM64
environment. The five complete tree-input alignments matched their reference
SHA-256 checksums; downstream tests retained the reported tree-sensitivity
directions and significance conclusions. The numerical tree comparisons are
provided in `verification/genome_reconstruction_reference.csv`.

A fresh installation of the genomics environment and execution on other
operating systems have not been validated. Parallel tree searches can produce
small secondary-estimate differences; code checks the reported conclusions
and records both estimates. `expected/alignment_checksums.tsv` specifies the
full-rebuild alignment checks.
