# Resource-access breadth marks protective strains against toxigenic *Clostridioides difficile*

This repository reproduces the study's mouse disease and protection scores,
strain-specific qPCR, GC-MS, BIOLOG PM1, phylogenetic, cytometry and microbiota
analyses. It generates five main figures and three supplementary figures from
the supplied experimental data and genome assemblies.

## Requirements

The standard analysis requires MATLAB R2024b with Statistics and Machine
Learning Toolbox and Bioinformatics Toolbox. It has been tested on macOS ARM64.
Other operating systems have not been validated. Genome reconstruction also
requires the external programs listed in [GENOMICS.md](GENOMICS.md).

## Run the analysis

Download a release or clone this repository into a new folder. In a fresh
MATLAB session, change to that folder and run:

```matlab
run_all
```

The standard route recalculates the experimental analyses, downstream
phylogenetic tests and all figures. It uses the derived trees supplied in
`cache/phylogeny/`, so it does not require genomics executables, the authors'
workspace or the shared drive. No generated experimental results are supplied
as calculation inputs.

To reconstruct the alignment, recombination filtering and trees from the
supplied genome assemblies, first install the genomics dependencies and set
`CDIFF_PHYLO_BIN` as described in [GENOMICS.md](GENOMICS.md), then run:

```matlab
run_all('full')
```

Full reconstruction takes hours and creates large intermediate alignments.
Use a fresh package copy for a from-scratch reconstruction.

Both routes stop on errors and write outputs to `results/`. To check an
existing run without repeating its calculations:

```matlab
run_all('tests')
test_figure1_identity
```

Checks cover input integrity, numerical results, statistical definitions,
Figure 1 mouse identities, and figure completeness. [TEST_REPORT.md](TEST_REPORT.md)
states the validation scope and the distinction between the standard and full
routes. Running the tests does not substitute stored results for calculations.

## Files and outputs

| Location | Contents |
|---|---|
| `input/` | Experimental records and 24 genome assemblies |
| `code/analysis/` | MATLAB analyses, figure generators and tests |
| `cache/phylogeny/` | Derived trees and associated records used by the standard route |
| `expected/` | Numerical reference values, supplied comparison controls and a synthetic alignment-parser fixture; used only for verification |
| `verification/` | Validation records for this release |
| `results/tables/` | Recalculated measurements, scores, tests and figure data |
| `results/analyses/` | PM1 model comparisons and disease/phylogeny-adjusted analyses |
| `results/figures/` | Five main and three supplementary PNG figures |
| `results/vectors/` | Editable PDF, SVG and MATLAB FIG exports |
| `results/diagnostic_figures/` | Additional analysis diagnostics, not additional manuscript figures |
| `results/logs/` | Execution logs |

[FIGURE_MAP.md](FIGURE_MAP.md) links each figure to its inputs and producing
code. [DATA_DICTIONARY.md](DATA_DICTIONARY.md) explains scoring, missing values,
assay conventions and sampling. [PM1_REPLICATION_SENSITIVITY.md](PM1_REPLICATION_SENSITIVITY.md)
describes the equal-replication check.

Starting data include cleaned animal records, integrated GC-MS peaks,
exported cytometry events, genus abundances and assembled genomes. The package
does not reproduce mass-spectral peak integration, FCS gating, raw-read
processing or genome assembly. Strain-specific fecal qPCR uses Kevin Sia's
KS65 records and calibration, not the separate in-vitro competition assay.
Its mixed inoculum targeted approximately 1:5 ST1-75:VPI10463; the plotted
inoculum guide is not a measured starting strain frequency.

Random seeds, resample counts, scoring rules and missing-value handling are
specified in code. Terminal-event zeros are analytical scores at recorded
terminal rows, not measured weights or added later observations. The PM1
matrix contains 700 positive calls across the study panel.

## Provenance and licensing

`input_manifest.tsv` lists every supplied input, its provenance, distributed
SHA-256 checksum and original source checksum. Public genomic accessions are
listed in [GENOME_ACCESSIONS.md](GENOME_ACCESSIONS.md).
Seven XLSX copies omit only Excel's remembered private folder location; their
other XML members, including cells, formulas and styles, are unchanged.
[CODE_PROVENANCE.md](CODE_PROVENANCE.md) documents the analysis implementations.

Code and documentation use the MIT license in [LICENSE](LICENSE), with the
upstream notice retained in [NOTICE.md](NOTICE.md). Original experimental and
study-derived data use CC BY 4.0 within the scope of [DATA_LICENSE.md](DATA_LICENSE.md).
Public genome assemblies and external dependencies retain their own terms.

Repository: https://github.com/joaobxavier/cdiff_resource_access.
The package version is recorded in `RELEASE_VERSION`. Cite a release tag or
specific commit to identify a reproducible snapshot. Generated results are
excluded from the source distribution.
