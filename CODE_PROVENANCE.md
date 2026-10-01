# Analysis provenance

`run_all.m` specifies the execution order for the MATLAB code in `code/analysis/`.
The workflow reads the distributed inputs and produces its own measurements,
models, tables and figures. It does not read manuscript text, hand-edited artwork
or previously generated figure objects.

| Analysis | Implementation and source |
|---|---|
| Mouse disease, protection, survival and fecal CFU | `paper/analyze_mouse_weight_survival_screen.m` and `paper/analyze_fecal_cfu_screen.m`, using Kevin Sia's supplied animal records |
| Fecal strain-specific qPCR | `qpcr/rebuild_qpcr.m`, implementing the quantity, pellet-mass and relative-fraction calculations in Kevin's KS65 workbook |
| Figure 1 mouse identities | `paper/label_ks65_mice.m`, matching explicit IDs and raw/relative weights in KS65 Biometrics to the plotted observations |
| PM1 plate calls | `biolog/analyze_biolog_plate.m`, from Vishwas Mishra's upstream repository at commit `4cb8408438b62d88605c640e3248899f58dd264b`; `biolog/rebuild_pm1.m` combines replicate calls by OR |
| GC-MS | `gcms/` matched-batch intracellular and extracellular analyses, independently fitted descriptive PLS-DA, and cross-assay comparisons |
| PM1 protection associations | `paper/prepare_figure4.m`, `paper/compare_pm1_protection_models.m` and `paper/analyze_amino_acid_breadth_specificity.m` |
| Equal-replication sensitivity | `biolog/assess_pm1_equal_replication.m`, enumerating one plate per strain |
| Protein-family and nutrient overlap | `ecology/calculate_pgfam_overlap.m` and `ecology/analyze_resource_overlap.m`; independently implemented MATLAB set coverage following Spragge et al. (2023), with strain-level rank tests, bootstrap intervals and plate/strain sensitivities |
| Resource-competition model | `ecology/resource_competition_outcomes.m`; overlap and total access parameterization, with exact legacy-equation equivalence checked for every strain and limitation-grid point |
| Phylogeny and ancestry adjustment | `phylogeny/` genome-coordinate alignment, recombination filtering, tree construction and disease/phylogeny-adjusted model tests |
| Cytometry | `host/flow_adaptive.m` and `host/flow_innate.m`, implementing the supplied range-filtered retained-event fraction calculations |
| Residual microbiota | `host/compute_residual_microbiome_stats.m`, excluding *Clostridioides* and renormalizing remaining genera |
| Main and supplementary figures | `paper/generate_figure1.m`, `paper/generate_figure2.m`, `gcms/generate_figure3.m`, `paper/generate_figure4.m`, `host/generate_host_figures.m`, and the fecal-CFU workflow |

Paths in this table are relative to `code/analysis/`. The original PM1 code is
available at https://github.com/vim4007/C.difficile_Protection; its license notice
is retained in `NOTICE.md`.

Strain order is fixed where resampling requires deterministic assignment of
random draws. Figure 1C/D identities are established by explicit KS65 identifiers
and weight records; the labeling helper asserts unchanged plotted observations.

Genome alignments exclude equally best NUCmer alignments with identical
endpoints but different gap placements, preserving the unique-alignment rule
regardless of parallel output order. Inferior alignments at another locus do
not invalidate a unique best alignment. IQ-TREE resume checks require a
completion log, not merely a provisional treefile. Source selection requires
the numbered ST1 assembly filenames. These rules and the alignment checksums
are documented in `GENOMICS.md` and code.

`expected/numerical/` contains reference results for regression tests, not
calculation inputs. `expected/source/` contains supplied PM1 and qPCR comparison
controls, also used only after the corresponding quantities are recalculated.
The synthetic alignment fixture tests the parser and is never included in a
study alignment. `cache/` is a labeled derived-tree shortcut; the full route
reconstructs it from genome inputs.

`test_code_consistency` independently checks the eight chemical-class tests
adjusted for total breadth against MATLAB `partialcorr`, using 18 degrees
of freedom. It also checks the equal-plate inferential p-values and the saved
Figure 4C/F and S2D contents. The auxiliary total-breadth-adjusted class p/q
columns were corrected from a plain residual-correlation test to the
one-covariate test; only those 16 frozen entries were updated after the native
cross-check. Raw class coefficients/p/q values, disease-adjusted amino-acid
statistics, all displayed values and scientific conclusions are unchanged.
Unused plotting functions for the replaced panels are absent from this
implementation; their prior versions remain in Git history.

The code/documentation and data licenses have separate scopes. See `LICENSE`,
`NOTICE.md` and `DATA_LICENSE.md`. Genomics executables are external dependencies
and are not redistributed.
