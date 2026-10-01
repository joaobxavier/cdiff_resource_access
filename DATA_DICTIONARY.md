# Data and analysis conventions

The exact source filenames, roles, hashes and provenance are enumerated in
`input_manifest.tsv`. All inputs are read-only to the workflow.

| Input family | Starting measurement and retained structure |
|---|---|
| Mouse screen CSVs | Animal, experiment, treatment/strain, day, relative body weight, terminal event and fecal CFU/g. Experiment/condition/mouse jointly identify animals. KS65 is separate from the panel-wide screen. |
| KS65 workbooks | Original weight observations and calibrated strain-specific qPCR quantity estimates, pellet masses, instrument and setup records. Recalculate fractions/copies per gram from Kevin's quantities; do not substitute in-vitro calibration or refit standard curves. |
| GC-MS | Intracellular and extracellular peak signals for five strains and media, three matched batches. Suffix `_i` matches Media_i. Intracellular blank correction subtracts that medium; extracellular continuous features use log2(sample/medium). Missing or nonpositive signals are not silently imputed. |
| PM1 raw plates | Time-series absorbance, well/substrate labels and chemical-class annotations. Keep Time < 24; subtract each well's initial absorbance and clip negatives to zero; fit the original mixed model. Positive substrate-related coefficients at p <= 0.05/95 define plate calls, combined by OR across replicate plates. |
| Cytometry | Supplied all-event and gated-event CSV exports for mice 1, 6, 7, 11, 12, 13. Match events with mouse identity, use the original range-filtered denominator, and calculate per-mouse fractions before group summaries. No re-gating from FCS is performed. |
| RAG1 | Observed post-challenge weight trajectories for previously uncolonized wild-type, ST1-75-colonized wild-type and ST1-75-colonized RAG1-deficient mice. Both effects share the uncolonized wild-type reference. No uncolonized RAG1-deficient arm is inferred. |
| Microbiota | Supplied genus-abundance table and animal/day/treatment metadata. Exclude Clostridioides and renormalize the remaining community for the reported comparisons. |
| Genomes | Twenty-one supplied ST1 assemblies, accession-verified R20291 and primary VPI references, and one alternative VPI assembly. See GENOMICS.md. |
| PGFam exports | Unmodified BV-BRC PATRIC feature tables for 21 mapped ST1 genomes and VPI10463 genome 1496.1555 (GCA_001995155.1). Deduplicate assigned CDS PGFams within genomes; retain assigned hypothetical proteins. See RESOURCE_OVERLAP.md. |

Terminal events include death or humane removal. Only a source row with missing
weight and a recorded terminal event receives a zero in the primary mouse
score; no later rows are invented. Plotted weights remain observed weights.
The RAG1 model uses observed weights and does not acquire zero-score rows.

Fecal non-detects remain zeros in the source and appear in the ND band, not at
an invented detection limit. Positive burdens are modeled on the log scale.
Strain-alone counts resolve the inoculated strain; co-colonization counts are
total C. difficile and do not resolve composition. Sampling schedules are
explicit in the CFU code. Late measurements are conditional on observation.

Sample order is fixed for resampling. In particular, the original 21-strain
order is retained for non-phylogenetic PM1 bootstrap/permutation tests; the
phylogenetic tests use their documented alphabetical order. This metadata is
not an independently authored biological result table. The new overlap
bootstraps use alphabetical strain order and seed 20261001; they do not
replace the original relative-amino-acid resamples.

In the extracellular analysis, zero signals use half the feature-specific
minimum positive signal only where a logarithm is required. Original signals
remain available for medium-range comparisons. Intracellular transformed
sensitivity analyses and cross-validated diagnostic classifiers are separate
from the above-background, full-data PLS-DA shown in Figure 3A.

The PLS-DA projections use all individual samples, are fitted independently
by compartment and are descriptive, not validated classifiers. The plotted
ellipses summarize observed sample dispersion, not confidence regions.

Stable numerical-output identifiers retain `Shared_With_VPI` for the number
of overlapping calls, `ST1_Private` for additional candidate access,
`VPI_Private` for pathogen access outside overlap, and `Private_Advantage`
for total relative access. The count `Shared_With_VPI` is not itself a
fraction: nutrient overlap divides it by the corresponding VPI-positive count.
Scalar-model diagnostic tables retain their identifiers for comparison across
releases; they do not designate the primary overlap-plus-access interpretation.
`DiseaseAdjustedRankP` in the equal-plate tables is the ordinary correlation
p-value for two residual vectors, retained as a diagnostic only. Use
`AdjustedRankP_ControlDF` (18 degrees of freedom) for inference, as in Figure 4E.
