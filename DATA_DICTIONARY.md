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
not an independently authored biological result table.

The PLS-DA projections use all individual samples, are fitted independently
by compartment and are descriptive, not validated classifiers. The plotted
ellipses summarize observed sample dispersion, not confidence regions.
