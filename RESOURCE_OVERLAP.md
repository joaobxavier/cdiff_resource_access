# Protein-family and nutrient overlap

## Inputs and identity

The two feature CSVs in `input/pgfam/` are unchanged author-supplied BV-BRC
exports, with PATRIC annotations and cross-genus protein-family (PGFam)
assignments. Genome IDs are read as strings. The 21 ST1 identities come from
the supplied strain mapping; its obsolete VPI row is not used. The VPI export
is genome **1496.1555**, assembly **GCA_001995155.1**, contigs MUJV01000001
and MUJV01000002. Exported coordinates are checked against that supplied FASTA.
This is not the primary phylogenetic outgroup GCF_015238635.1. Annotation
snapshot dates/versions were not supplied in the exports.

Protein-family overlap is the number of distinct VPI-assigned CDS PGFams also
present in the candidate divided by the 2,983 distinct VPI PGFams. Repeated
families count once; missing assignments are excluded; assigned hypothetical
proteins are included. No metabolic-function filtering is applied.
This implements Spragge et al.'s pathogen-centered coverage with one candidate
as the resident set (Science 2023, doi:10.1126/science.adj3502). It does not
repeat the community experiments or Bakkeren et al.'s candidate-centered
unique-family analysis (Nature Microbiology 2025, doi:10.1038/s41564-025-02162-w).
The set implementation was independently checked against all 156 communities
and 468 values in the Spragge code/data repository at commit
f9181dc2088094f71a91d24b6323b54d92b9d67d (maximum difference about 5e-16).
These upstream validation data are not required study inputs.

## Analyses

`analyze_resource_overlap` reads freshly reconstructed PM1 calls and disease/
protection estimates, not frozen expectations. It verifies 700 positive calls,
95 substrates, 16 amino-acid-class substrates and all 21 strain identities.
Nutrient overlap divides candidate/VPI shared positive calls by VPI-positive
calls: 23 overall, 3 amino-acid-class and 20 remaining calls. The amino-acid
metric has four possible values: 0, 1/3, 2/3 and 1.

Nine rank tests are retained. BH correction is separate for the three PGFam
tests (raw overlap, disease-adjusted overlap, relative amino-acid access
adjusted for disease and PGFam overlap) and six nutrient tests (three classes,
raw and disease-adjusted). Existing correction families are unchanged.
The extra PGFam-adjusted amino-acid test is exploratory, not a validation claim.

Average ranks, intercept-plus-covariate residualization, and n-minus-covariates-
minus-two degrees of freedom match MATLAB corr/partialcorr. Alphabetically
ordered strains are resampled together 5,000 times with seed 20261001; ranks
are recomputed per draw. Percentile intervals omit only explicitly invalid
draws. One draw is constant for amino-acid overlap; its raw and adjusted
estimates are recorded as missing, leaving 4,999 valid draws for those tests.
These intervals resample strain-level point estimates and do not propagate
mouse-score estimation covariance or adjust for phylogeny.

Outputs retain all 189 leave-one-strain-out results, the nine full bootstrap
series, native MATLAB checks, and all 81 one-plate selections. Each selection
recomputes its VPI denominator from the selected VPI plate; zero denominators
would be reported as invalid rather than imputed. All 81 observed selections
are valid for all three measures. The pooled matrix remains the primary analysis.

## Model and verification

The model inputs are the shared-resource count and each strain's total access.
Access outside overlap is derived internally. All 10,521 strain-by-grid
outcomes agree exactly with the previous equations, including feasibility
before invasion signs. The original regression checks are retained, with
the auxiliary class-adjusted p/q correction documented in CODE_PROVENANCE.md;
new overlap expectations supplement rather than replace these checks.

BV-BRC annotations are third-party source data and retain their source terms;
the package's CC BY notice does not relicense them. Original and distributed
hashes are recorded in input_manifest.tsv.
