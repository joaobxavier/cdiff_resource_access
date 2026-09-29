# Genome reconstruction

Public strain-level read accessions, reference assemblies, and their analytical
roles are listed in [GENOME_ACCESSIONS.md](GENOME_ACCESSIONS.md), corresponding
to manuscript Supplementary Table S1. Exact input-file sources and SHA-256
checksums remain in `input_manifest.tsv`.

The default run uses the derived trees in `cache/phylogeny/tree_tests`. Those
files were built from the genomes in `input/genomes`; they are not independent
experimental inputs. `run_all('full')` reconstructs them and the full downstream
analysis. Use a fresh package copy for a from-scratch reconstruction.

External dependencies: MUMmer 4.0.1, SKA 2.5.1, Gubbins 3.4.3, and IQ-TREE 2.4.0.
`environment-osx-arm64-explicit.txt` records the original macOS ARM environment.
It is platform-specific, not a portable Linux lockfile. On other systems install
these versions with their dependencies. Set `CDIFF_PHYLO_BIN` to their `bin`
directory. MATLAB additionally needs Bioinformatics Toolbox.

On macOS ARM64, create an environment with a Conda-compatible installer using
`micromamba create -p /your/environment/path --file environment-osx-arm64-explicit.txt`.
Then in MATLAB use `setenv('CDIFF_PHYLO_BIN','/your/environment/path/bin')`.
The explicit file pins package builds; environment installation still requires
network access. The analysis itself uses the supplied local input files.

The workflow maps 21 ST1 assemblies and two alternative VPI assemblies to the
accession-verified R20291 reference. It retains unique one-to-one alignments
of at least 500 bp and 95% identity, masks ambiguous/overlapping positions and
indel flanks, reference repeats, and disagreements with independent SKA calls.
Equally best alignments with identical endpoints but different gap placements
are explicitly excluded. MUMmer's one-to-one filter otherwise handles these
ties differently when parallel output changes their order. This explicit rule
recovers the original study alignment without adding or changing any base call.
No missing base is filled with the reference. Gubbins runs on the ST1 alignment
only. Its inferred recombinant interval union is then removed from every taxon.
Complete-column alignments retain invariant sites. IQ-TREE uses ModelFinder,
1,000 ultrafast bootstraps with refinement and 1,000 SH-aLRT replicates, four
threads, and seed 20260919.

Gubbins reached the 20-iteration cap under its weighted-tree stopping rule;
the retained recombination interval union was stable. The code records this
diagnostic rather than treating the cap as proof of convergence. Five tree
variants assess filtering, outgroup, and common-column sensitivity.
Branch assignments can yield a different count of event rows on poorly resolved
branches while producing identical excluded positions and retained sequences.
The release tests the scientific outputs and all five complete tree-input
alignments, not equality of that unreported event-row count.
Parallel tree searches and support calculations can differ slightly despite
the fixed seed. Main-tree manuscript statistics have strict numerical checks;
the 20 secondary tree/trait comparisons additionally verify the reported
association directions and significance conclusions, with both sets of
estimates written to `results/tree_sensitivity_comparison.csv`.

Gubbins helper programs need a whitespace-free scratch path. The default is a
new temporary directory; optionally set `CDIFF_PHYLO_WORK` to an empty scratch
directory. Existing outputs are reused for interrupted runs. The workflow never
deletes user directories. Large intermediate alignments remain in `results/`
and scratch, not in the distributed source package.
