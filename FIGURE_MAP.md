# Figures and producing analyses

| Figure | Inputs | Producing code and outputs |
|---|---|---|
| 1: ST1-75 protection and enrichment | Mouse protection screen; KS65 weights and calibrated qPCR quantities | Mouse-score/survival analysis, `rebuild_qpcr`, `generate_figure1`; tables prefixed `figure1_` |
| 2: Disease and protection screen | Mono-colonization and co-colonization animal tables | Mouse-score and scheduled-CFU analyses, `generate_figure2`; tables prefixed `figure2_` |
| 3: Metabolic programs | `intra.xlsx`, `secretome.xlsx`; reconstructed PM1 calls for cross-assay comparisons | Matched-batch GC-MS analyses and `generate_figure3`; tables prefixed `figure3_` and detailed GC-MS tables |
| 4: Phylogeny and relative resource access | Genomes/derived trees, reconstructed PM1 calls, chemical classes and mouse scores | PM1 traits/models/specificity, tree tests and `generate_figure4`; tables prefixed `figure4_` |
| 5: Host and microbiota context | RAG1 workbook, 90 cytometry CSVs and genus abundances | Cytometry fractions, residual-microbiota analysis and `generate_host_figures`; tables prefixed `figure5_` |
| S1: Fecal CFU | Both animal tables | `analyze_fecal_cfu_screen`; tables prefixed `figureS1_` |
| S2: PM1 robustness | Reconstructed PM1 calls, chemical classes and mouse scores | Specificity, well/strain influence, resampling and prediction tests; `generate_figure4`; tables prefixed `figureS2_` |
| S3: Residual microbiota | Genus-abundance table | `generate_host_figures`; A: Bray–Curtis community differences, B: day-1 genus contrasts |

The package produces eight PNG figures with editable PDF, SVG and MATLAB FIG
exports. Figure 4 contains the phylogeny; Figure 5E contains residual Shannon
diversity. Figure 5C/D fractions use all retained events in the relevant staining
panel as their denominator.

Figure 1C/D mouse identities 3–5 are linked by explicit identifiers and raw/relative
weights in `input/qpcr/Analysis/KS65.xlsx`, sheet Biometrics. The generator checks
these against the plotted observations, uses matching colors, and writes
`results/tables/main/figure1_ks65_mouse_identity.csv`.

The disease-and-phylogeny-adjusted amino-acid analysis has no additional figure.
It writes model comparisons, score and shared-branch covariance matrices,
leave-one-strain-out results, and a 1,999-replicate null-refitted parametric
bootstrap under `results/analyses/phylogeny_adjusted_amino_access/`.

`run_all.m` lists all producing steps in order. Detailed diagnostic plots are
written separately from the main and supplementary figures.
