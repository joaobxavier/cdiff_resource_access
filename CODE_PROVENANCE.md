# Code provenance

Draft 44 changes only figure presentation and export checks: PM1 robustness is
S2; the residual-microbiota supplement is S3A,B without the repeated Shannon
panel; the duplicate phylogeny supplement is retired. Figure 5C,D axes name
the retained-event denominator already implemented by the fraction calculation.
All scientific inputs, calculations, and frozen numerical expectations are
unchanged. Historical figure-number stems in numerical table filenames remain
stable. The retired supplementary tree generator is preserved in the authors'
historical source repository and omitted from this active package.

Draft 43 adds mouse-identity labels and matching qPCR colors to Figure 1C.
`label_ks65_mice.m` verifies the explicit mouse identifiers in the original
`input/qpcr/Analysis/KS65.xlsx` Biometrics sheet against both its raw and
relative weights, then matches those values to the plotted observations.
It asserts that no plotted X/Y values change. No statistical calculation or
numerical expectation is changed. The existing KS65 input is sufficient;
no derived author table or historical figure is needed by the package.

The maintained release code is the `code/` directory, executed by `run_all.m`.
It preserves the verified MATLAB statistical calculations and replaces private
workspace paths, copied earlier-draft result tables and saved figure objects
with explicit producing steps. Version-numbered function names retain their
implementation provenance; only the active analysis chain is distributed.

| Release component | Verified implementation origin |
|---|---|
| Mouse scores, survival and CFU | Study MATLAB `analysis/manuscript_draft/analyze_mouse_weight_survival_screen.m` and `analyze_fecal_cfu_screen.m` |
| Figure 1–3 displays | Author-approved Draft 39 MATLAB generators, with the current single-comparison Figure 1 survival annotation |
| PM1 plate calls | Vishwas GitHub snapshot `4cb8408438b62d88605c640e3248899f58dd264b`, `analysis/fig5/analyze_biolog_plate.m`; package combines replicate calls directly |
| GC-MS | Verified MATLAB matched-batch workflows in `analysis/gcms_exploration/` |
| PM1 associations and model | Draft 42 Figure 4 preparation/generator plus the verified simple-model and amino-acid-specificity workflows |
| Phylogeny and ancestry adjustment | Verified 2026-09-19 MATLAB genome-coordinate, recombination and tree-test workflows |
| Cytometry | Fraction-calculation portions of `analysis/fig2/fig2_adaptive_analysis.m` and `fig2_innate_analysis.m`; no unused tSNE/PCA |
| Microbiota | Verified residual-community statistical functions, originally implemented in the Draft 11 figure workflow and retained in subsequent drafts |
| Host figures | Five-panel main display and two-panel microbiota S3, composed from regenerated fractions/results, with Draft 44 presentation corrections |
| qPCR | Direct implementation of Kevin's supplied KS65 worksheet formulas and sample assignments |

Packaging changes do not change scoring, models, thresholds, multiplicity
families or experimental observations. They remove the retired accessory-gene
PCA and earlier figure/result dependencies. The fixed original strain order is
recorded in code to reproduce historical resampling draws. One supplementary
plot label is moved clear of its interval; its numbers and data are unchanged.
The CFU display uses the current S1 colors and footer from
`generate_theory_led_thirty_first_draft_figureS5.m`, rather than the earlier
exploration's colors. Its statistical calculations are identical.

Clean-room testing identified order-dependent treatment of two equally scored
NUCmer gap placements. The package now explicitly excludes tied alignments,
implementing the original unique-alignment requirement deterministically. The
resulting 23-genome alignment was compared base by base with the original and
was identical. An inferior alternative at another locus is not excluded when
the best alignment is unique. IQ-TREE resume checks also require its completion
log, not merely a provisional treefile created during bootstrap refinement.
When resuming, source-assembly selection explicitly requires numbered ST1
filenames so that case-insensitive macOS matching cannot include a lower-case
Gubbins output as an additional genome.

`expected/numerical/` contains frozen manuscript-result checks, not independent
inputs to scientific calculations. `expected/source/` contains the supplied
combined PM1 matrix and plotted qPCR worksheet, also used only after rebuilding
their quantities. `cache/` is a separately labeled derived phylogeny shortcut.

The study software is licensed under MIT, retaining the verified upstream
Vishwas Mishra notice. See LICENSE and NOTICE.md. DATA_LICENSE.md separately
defines the CC BY 4.0 scope for original experimental and study-derived data,
without relicensing supplied/public assemblies or third-party software.
The pinned environment references external software; no executables are
redistributed. Public deposition has not yet occurred.
