# Resource-access breadth marks protective strains against toxigenic *Clostridioides difficile*

Reproducibility package for Draft 44, retaining the Draft 42 analyses and
source-verified mouse labels/colors added to Figure 1C in Draft 43. Draft 44
removes the duplicate phylogeny supplement and Shannon-diversity panel and
labels flow-cytometry fractions by their retained-event denominator. Source data are
unchanged copies supplied by Vishwas Mishra/Kevin Sia or accession-verified
public assemblies. `input_manifest.tsv` gives provenance and SHA-256 hashes.

## Run

In a fresh MATLAB R2024b session with Statistics and Machine Learning Toolbox and Bioinformatics
Toolbox, change to this folder and run `run_all` from the Command Window.
The default run recomputes
the experimental analyses and figures, using the explicitly labeled derived
tree cache. It does not require the authors' workspace or shared drive.

`run_all('full')` additionally rebuilds the genome-coordinate alignment,
recombination filtering, and maximum-likelihood trees from the supplied genomes.
Set `CDIFF_PHYLO_BIN` to the directory containing the genomics executables.
See `GENOMICS.md` for dependencies and details. Full reconstruction can take
hours and produces large temporary alignments; these are not release inputs.

`run_all('tests')` checks an existing run. Results are written under `results/`.
The run stops on errors rather than silently substituting stored results.

`test_figure1_identity` additionally checks the Figure 1 numerical outputs
against frozen expectations after its generating stages have run. The plotting
helper itself checks mouse IDs against the original KS65 Biometrics worksheet
and asserts unchanged observations. See `verification/2026-09-28_figure1_identity/`
for this affected-stage check; it does not replace the earlier full-run record.

## Contents and scope

- `input/`: experimental data and 24 assembled genomes.
- `code/`: MATLAB analysis and plotting functions.
- `expected/`: independent supplied PM1/qPCR controls and numerical regression
  expectations, plus a clearly synthetic alignment-parser test fixture. These
  are checked after calculation, never used to fit study models.
- `cache/`: derived phylogenetic trees for the default run, explicitly separate
  from primary inputs and reproducible with the full route.
- `verification/`: retained test summaries, not inputs to any analysis.
- `results/`: generated tables, figures, logs, and test results; safe to omit
  when redistributing the source package.

Starting data include cleaned animal records, integrated GC-MS peaks, exported
gated cytometry events, genus abundances, and assembled genomes. This package
reproduces the analyses from those products, not mass-spectral integration,
FCS gating, raw-read processing, or genome assembly. KS65 uses Kevin's fecal
experiment, not the dissertation's in-vitro qPCR calibration.

All random seeds, resample counts, scoring and missing-value rules are stated
in code. Terminal-event zeros are analytical scores at recorded terminal rows,
never measured body weights or invented later observations. PM1 uses the
audited 700-call rule, not the legacy 881-call matrix. Supplementary figures
end at S3: S1 fecal CFU, S2 PM1 robustness, and S3 residual microbiota
(A: community composition; B: genus abundances). Shannon diversity appears
once in Figure 5E; phylogenetic context appears in Figure 4. The package
generates exactly five main and three supplementary PNG/PDF/SVG/FIG sets.
Retired accessory-gene and sequential-challenge analyses are omitted.

Some numerical table filenames retain older figure numbers (for example,
`figureS3_amino_nonamino_specificity.csv` supplies current S2). Their names
preserve calculation provenance; the current figure exports and checks use
Draft 44 numbering. Frozen numerical expectations have not changed.
On rerunning an older results directory, superseded figure exports are moved
to `results/historical_figure_exports/` before the active figures are rebuilt.

## Licensing and release status

Study analysis code and its documentation use the MIT license in `LICENSE`,
with Vishwas's upstream notice retained in `NOTICE.md`. Original experimental
and study-derived data use CC BY 4.0 within the scope in `DATA_LICENSE.md`;
public/supplied genome assemblies retain their original provenance and terms.
External dependencies retain their own licenses and are not bundled.

Public repository: https://github.com/joaobxavier/cdiff_resource_access.
This source package corresponds to Draft 44; use a specific commit when citing
or reproducing a frozen version.
Consult TEST_REPORT.md for the actual verification status.
Generated outputs are excluded by `.gitignore`; source data, code, explicit
tree cache and test expectations are retained. The supplied genomics environment
is for macOS ARM64; other operating systems have not been validated.
