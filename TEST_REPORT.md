# Verification report

## Public-package metadata preparation: 29 September 2026

Public copies omit private machine and shared-drive locations from five IQ-TREE
reports, a MATLAB usage comment, provenance columns and historical test-folder
descriptions. Seven XLSX copies also omit Excel's remembered private folder
location from `xl/workbook.xml`; all other workbook members are unchanged,
including every worksheet, formula, style and recorded value. Original source
hashes remain in `input_manifest.tsv` and original files in the authors'
historical archive. These edits do not change expected numerical values,
fitted trees, model parameters or executable MATLAB statements.
The byte-for-byte description below refers to the earlier,
pre-redaction archive; the publication copy includes these metadata-only edits
and refreshed release checksums.

After cleanup, MATLAB verified all 158 distributed input hashes, reconstructed
the eight measured qPCR fractions, and passed all 35 numerical checks (11,026
values) plus the KS65 identity checks. Static analysis still found no parser
errors across 48 files. This was affected-stage verification using outputs
from the complete run below, not another uninterrupted 21-stage run.

## Latest complete standard-route test: 29 September 2026

The exact Draft 44 release ZIP was extracted to a new temporary folder with
no `results/` directory and no seeded experimental outputs. Its SHA-256 was
`a95171d5d043bc3d8fc335109701bac1db016f5f04910d15c4a84975baac676c`.
All 21 standard-route stages passed, totaling 393.397 seconds of stage time.
All 35 numerical checks passed (11,026 values); all 158 supplied input hashes
matched. The additional KS65 checks passed three comparisons (266 values).
All eight PNG/PDF/SVG/FIG figure sets were regenerated and visually inspected.
MATLAB static analysis found no parser errors across 48 files; 19 nonfatal
diagnostics remain. The 1,999-replicate ancestry-adjusted bootstrap reproduced
p = 0.013. Evidence is in `verification/2026-09-29_release_end_to_end/`.

This was a complete default-route run using the explicitly distributed tree
cache, not a fresh genome-to-tree reconstruction. Author workspace paths were
removed while preserving MATLAB MCP infrastructure. The MCP request timed
out after 300 seconds, but MATLAB continued: all 21 PASS records and the final
regression report were saved. A subsequent MCP call confirmed completion and
reran both regression and KS65 checks successfully.

Content review found no credential-pattern hits, private correspondence,
generated results, manuscript drafts, manually edited AI files, bundled
executables or installed dependency trees. All archive entries were enumerated
and checksummed. The largest file was 6,010,705 bytes. The XLSX structure scan
found no macros, embedded objects, external-link parts or comment parts; a
targeted input scan found no email addresses or common clinical-identity labels.
These checks are not a legal clearance or a guarantee of absence of personal
information. One stale figure-map document was corrected to S1–S3.

MIT was selected for study code/documentation, retaining the verified original
Vishwas Mishra notice. DATA_LICENSE.md defines CC BY 4.0 for author-supplied
experimental and study-derived data, excluding public/supplied genome assembly
rights and external software. LICENSE and NOTICE.md record the code terms.
No institutional ownership determination or GitHub upload was made.

The final ZIP adds license files, the corrected figure map and these verification
records. All executable MATLAB code, scientific inputs, expected values and
tree-cache files are byte-identical to the exact tested ZIP. Final archive
integrity is checked after rebuilding. The tested pre-license ZIP is retained
as `cdiff_resource_access_reproducibility_2026-09-29_prelicense_tested.zip`.

## Latest release verification: 29 September 2026 — Draft 44

The maintained package and release archive now follow Draft 44: five main
figures and three supplements. S2 is the PM1 robustness figure; S3 contains
only the Bray–Curtis and genus-level microbiota panels. Figure 5C,D use
“Fraction of retained events.” The duplicate phylogeny supplement and repeated
Shannon panel are no longer generated. The scientific analyses, supplied
inputs and numerical expectations are unchanged.

Verification used an isolated package copy at
`an isolated temporary package copy`, with author-analysis
paths excluded and the MATLAB MCP infrastructure preserved. Upstream results
were seeded from the verified 20 September corrected run. Figure 1 was rerun
to check the source-verified KS65 mouse labels; the Figure 4/S2 and Figure 5/S3
generators were then rerun. This is affected-stage verification, **not** a new
uninterrupted 21-stage run or genome reconstruction.

- All 158 supplied input files passed SHA-256 validation.
- All 35 numerical regression checks passed (11,026 values).
- Figure 1 identity verification passed three additional comparisons
  (266 values; maximum absolute difference 9.9476e-14).
- MATLAB static checks covered 48 distributed files, with no parser errors
  and 19 retained style/performance diagnostics.
- Presentation checks confirmed exactly eight active PNGs with editable
  PDF/SVG exports, the two retained-event axis labels and the two-panel S3.
- The regenerated Figure 4, S2, Figure 5 and S3 were visually inspected.

Reports are in `verification/2026-09-29_draft44_cleanup/`. Generated outputs
are preserved separately in the workspace under
`reproducibility/development_test_outputs/2026-09-29_draft44_synchronized/`.
The release checksum manifest and ZIP have been refreshed; the previous ZIP
is preserved as a dated historical archive. Historical evidence below describes
its own test scope and figure numbering, not the current release.

## Latest correction verification: 20 September 2026

The S3 inconsistency identified below is resolved. The specificity and
Figure 4/S2/S3 stages were rerun in a separate copy using the verified upstream
products of the complete standard-route run. S3B now reports amino-acid
p = 0.0364940843264457 and non-amino-acid p = 0.730049043748984, using 18 df.
The two S3 expectation tables were updated only in their partial-rank p columns,
including the leave-one-well tests; all other numerical columns were verified
unchanged within the existing tolerance. Historical expectations are retained
under verification/2026-09-20_s3_correction/.

All 35 regression checks passed again (11,026 numerical values). New independent
formula and Figure 4E/S3B equality assertions passed. A negative-control test
temporarily substituted the stale S3 values in the isolated copy and confirmed
that the new check rejected them; the corrected output was restored and the
full suite passed. Static checks found no parser errors. The two modified local
analysis/generator scripts also passed MATLAB Code Analyzer without issues.
The regenerated S3 was visually checked and promoted to the manuscript together
with its editable exports. No other official figure was replaced.

Full corrected outputs: reproducibility/development_test_outputs/2026-09-20_s3_corrected/.
This was targeted post-fix verification, not another full genome reconstruction
or a new uninterrupted 21-stage run. The source-folder checksum manifest was
refreshed; the older release ZIP was not rebuilt.

## Earlier internal rerun: 20 September 2026 (pre-correction)

The current package completed the standard route in a new isolated copy at
`an isolated temporary package copy`, with MATLAB's search path reset to
defaults before adding that copy. All 21 stages passed, including the new PM1
equal-replication stage. All 35 regression checks passed (11,026 finite numerical
values), all nine PNG/PDF/SVG figure sets were generated, and all 46 distributed
MATLAB files were statically checked without parser errors (22 retained
diagnostics). The original input-integrity checks passed. The ancestry-adjusted
bootstrap reproduced p = 0.013 with 1,999 null replicates. The four numerical
PM1 sensitivity/baseline/predictor tables matched the approved 20 September
analysis byte for byte.

This was a standard-route rerun using the labeled tree cache, not a new
genome-to-tree reconstruction. The full-route evidence below remains dated
19 September. The MATLAB MCP request timed out after 300 seconds, but MATLAB
continued and wrote all 21 PASS records and the final regression report;
the analysis stages totaled approximately 312 seconds. A subsequent optional
diagnostic request encountered an MCP `mcpEval` signature error and did not
execute. Completion is established by the saved run products, not that request.

**Consistency finding still open:** both the official and regenerated S3B
show the former adjusted amino-acid p = 0.0315, whereas Figure 4E and the current
manuscript use the corrected p = 0.0364940843264458 (18 df). The S3 numerical
expectation also retains the former value, so a regression pass does not resolve
this cross-figure inconsistency. Its non-amino-acid partial-rank p uses the same
older residual-correlation convention and should be reconciled in the same
targeted correction. Correlations and the biological conclusion are unaffected
by the known amino-acid correction. No manuscript, figure or analysis-code
correction was made during this check.

Fresh test records are under `verification/2026-09-20_standard/`. Full generated
outputs are retained in the workspace at
`reproducibility/development_test_outputs/2026-09-20_standard_rerun/`.
The detailed manuscript-consistency report is
an internal review retained by the authors.
The existing release ZIP and release checksums have not been refreshed; this is
verification of the current folder, not certification of that older archive.
Resolve the S3 discrepancy before generating a new release archive.

## Historical full-route and standard-route verification

Date: 19 September 2026. Tested on macOS ARM64 using MATLAB R2024b Update 6
(24.2.0.2923080), Statistics and Machine Learning Toolbox, and Bioinformatics
Toolbox. Genome reconstruction used the package builds recorded in the explicit
environment file. A fresh installation of that environment and other operating
systems have not been tested.

## Completed checks

| Check | Outcome |
|---|---|
| Supplied input integrity | All 158 supplied files match their SHA-256 manifest; 180,295,708 input bytes, including the two independent supplied validation files |
| Analysis stages | All 20 stages completed successfully in the isolated full-route test copy |
| Numerical regression | 33 reference tables and two cytometry arrays passed: 11,026 finite numerical values |
| Missingness and identifiers | Missing-value positions and the tested strain, metabolite, model and comparison identifiers preserved |
| PM1 reconstruction | All 700 positive calls recovered directly from the 30 supplied plate workbooks and matched to Vishwas's supplied combined matrix |
| KS65 qPCR | Eight measured mixed-infection fractions recovered from Kevin's original quantity estimates; mouse 5's missing day-1 observation preserved |
| Genome reconstruction | All five complete tree-input FASTA files match the original analysis byte for byte, verified by SHA-256 |
| Phylogenetic sensitivity | All 20 secondary tree/trait comparisons retain the reported association directions and significance conclusions; actual and reference estimates are recorded side by side |
| Figures | Five main and four supplementary PNGs regenerated and visually inspected; PDF and SVG exports exist for all nine |
| MATLAB static analysis | All 45 MATLAB files checked; no parser errors. The 22 retained diagnostics concern style, unused arguments/variables and performance |
| Alignment parser | Synthetic test distinguishes equally best alternative gap placements from identical duplicates and uniquely better alignments |

Numerical tolerances are 1e-7 absolute plus 1e-7 times the reference magnitude.
The reference tables and arrays are used only for checks, never to fit models
or populate plotted results. Main-tree manuscript results retain these strict
checks; secondary tree sensitivity checks assess the qualitative claims stated
in the manuscript and retain both sets of estimates for inspection.

The full reconstruction recovered 3,133,615 retained primary-tree alignment
positions and 45 non-recombinant ST1-75/ST1-68 SNPs across jointly callable
positions. The ancestry-adjusted parametric-bootstrap result reproduced its
reference p-value of 0.013 with 1,999 null replicates.

## Test execution and fixes

Tests ran in temporary package copies outside the manuscript workspace, with
the manuscript analysis directories removed from MATLAB's search path. The
full-route test rebuilt the genomes, alignments, recombination filtering and
five trees, without substituting the supplied tree cache.

Testing exposed and corrected two reproducibility problems: order-dependent
NUCmer handling of equally scored alignments, and case-insensitive filename
matching that could select a Gubbins product when resuming on macOS. The fixes
recover the original alignment exactly and do not alter study observations or
reported statistical results. See CODE_PROVENANCE.md and GENOMICS.md.
The regression reader also uses an explicit CSV delimiter and rejects checks
that unexpectedly contain no numeric columns, preventing punctuation in prose
fields from weakening numeric verification.

Stalls in the long-lived MATLAB/MCP graphics session required restarting
MATLAB. The completed genome reconstruction was preserved and the
remaining workflow resumed outside interactive figure capture. Following the
restart-selection correction, the final complete numerical and alignment
checks passed. The full-route evidence is therefore a completed, resumed
clean-room reconstruction, not a claim that the first uninterrupted attempt
succeeded.

The default workflow subsequently completed all 20 stages without interruption
in a fresh MATLAB session and a new package copy, passing all 35 checks and
11,026 numerical comparisons. Final visual comparison with the official artwork
identified the older CFU color scheme and footer in the starting analysis file;
these were updated to the current S1 presentation. The CFU stage and complete
regression suite were then rerun for both test copies, with unchanged numbers.
The retained `verification/` reports document the final checks and stage runs.

Gubbins assigned 41 event rows rather than the original 40, while retaining
identical tree-input sequences and the same 319,251-base interval-union length.
Parallel IQ-TREE support calculations also differed slightly. These unreported
intermediate differences are disclosed, not forced to match cached values.
The main-tree statistics, ancestry-adjusted analyses, and all reported
phylogenetic sensitivity conclusions passed their stated checks.

## Boundaries

This package starts from supplied processed animal records, integrated GC-MS
signals, gated cytometry exports, genus-abundance tables and assembled genomes.
It does not reproduce upstream spectral integration, FCS gating or raw-read
assembly. The standard run intentionally uses the labeled derived tree cache;
the full route rebuilds it. Neither route needs the authors' shared drive.

All nine regenerated PNGs were checked for panel composition, labels, clipping
and preservation of the author-approved presentation. The existing manuscript,
official figure assets, original analysis files and manually edited Illustrator
files were not changed. No repository was published and no licensing decision
was made.
