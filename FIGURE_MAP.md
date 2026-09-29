# Figure-to-analysis map

| Figure | Main inputs | Producing analysis |
|---|---|---|
| 1: ST1-75 protection and enrichment | Mouse protection screen; KS65 weights and original qPCR quantity estimates | Mouse scores/survival, KS65 adapter, Figure 1 generator |
| 2: Disease/protection screen | Both animal tables | Mouse scores, scheduled CFU extraction, Figure 2 generator |
| 3: Metabolic programs | `intra.xlsx`, `secretome.xlsx`; regenerated PM1 calls for cross-assay comparisons | Matched-batch GC-MS workflows and Figure 3 generator |
| 4: Phylogeny and relative resource access | Genomes/derived trees, regenerated PM1 calls, chemical classes, mouse scores | PM1 traits/models/specificity, tree tests, Figure 4 generator |
| 5: Host and microbiota context | RAG1 workbook, 90 cytometry CSVs, genus-abundance table | Flow fractions, residual microbiota, host figure generator |
| S1: Fecal CFU | Both animal tables | Scheduled CFU analysis |
| S2: PM1 robustness | Rebuilt PM1 calls, classes, mouse scores | Figure 4 workflow; specificity, influence, resampling and held-out prediction |
| S3: Residual microbiota | Genus-abundance table | Host workflow; A: Bray–Curtis group effects, B: day-1 genus-level contrasts, excluding Clostridioides and renormalizing remaining genera |

These are the eight active Draft 44 figures. Phylogenetic context appears in
Figure 4 rather than a duplicate supplementary figure. Shannon diversity
appears only in Figure 5E. Figure 5C,D label fractions of retained events.
Some numerical-table filenames retain historical figure numbers, as documented
in README.md; they do not indicate additional active figures.

Figure 1C/D mouse identities 3–5 are linked by the explicit identifiers and
weights in `input/qpcr/Analysis/KS65.xlsx`, sheet Biometrics. The generator
checks these against the plotted weights and uses the same mouse colors in
both panels. It writes `results/tables/main/figure1_ks65_mouse_identity.csv`.

The additional disease-and-ancestry-adjusted amino-acid analysis writes model
comparisons, full score covariance, shared-branch covariance, leave-one-strain-out
results, and the 1,999-replicate null-refitted parametric bootstrap. It has no
additional figure. `run_all.m` is the authoritative execution order.

Historical version numbers retained in a few MATLAB filenames identify their
verified implementation origin, not additional manuscript versions or inputs.
No old manuscript text, hand-edited Illustrator file, or saved MATLAB figure
is read by the analysis.
