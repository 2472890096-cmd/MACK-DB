# PRJEB81106: complete archived 33-Run matrix inputs

This small local release-candidate bundle supplies a complete BioProject matrix starting point: 53 species × 33 archived Run samples and 38 historical genus groups × the same 33 samples. It is not the current full public BioProject cohort. No FASTQ, SRA, Sketch, large reference database, model checkpoint or saved neural-network prediction is included.

## Inputs and order

| File | Meaning |
|---|---|
| `inputs/species/abundance.matrix` | Original tab-separated 53 × 33 abundance values; first header is `OTU_id`, sample columns are IDs `0`–`32`. |
| `inputs/species/taxa.tsv` | Headerless 0-based taxon ID and exact original species label. |
| `inputs/genus/abundance.matrix` | Original tab-separated 38 × 33 historical genus-group values, with the same sample order. |
| `inputs/genus/taxa.tsv` | Headerless 0-based taxon ID and exact original genus-group label. |
| `metadata/runs.txt` | Ordered 33 Run accessions; column ID `i` is line `i+1`. |
| `metadata/sample_mapping.tsv` | Explicit 0-based `sample_index` → `run` correspondence. |
| `reference/species_taxonomy.tsv` | Original GTDB R226-derived taxonomy header and 53 exact matching species records. |
| `reference/genus_taxonomy.tsv` | Original GTDB R226-derived taxonomy header and 38 exact matching genus records. |

The matrices and mappings were copied byte for byte. No normalization, transposition, relabeling or rounding was performed. `INPUT_FORMAT_CHECKS.json` records dimensions, name/ID checks, finite nonnegative values and the 33 ordered accessions. It is a format check, not a new scientific run.

Genus values retain the author-adopted historical species-label token aggregation. This groups the token extracted by `label.split(' ')[0].split('_')[1]`, including folding some GTDB suffixes. It is not strict grouping by complete GTDB genus strings. This bundle does not regenerate or alter the archived genus matrix.

## Minimal taxonomy references

The accepted Zi-Pi program reads semicolon-separated taxonomy with the first column as row names. All nine original columns and their values are retained unchanged. For species, matching removes `s__` and replaces spaces with underscores, exactly as the program does; the first full-table occurrence of each requested species is retained. For genus, the first original row whose `Genus` equals `g__` plus a requested historical genus label is retained. Requested names all have one matching reference. Source row order is preserved.

Separate references matter: the species subset contains GTDB names such as `Bacillus_A`, whereas the historical matrix can contain `Bacillus`. Supplying the genus-specific subset follows the original full-table lookup; it does not change taxonomic names to force a match. `SOURCE_PROVENANCE.tsv` identifies the full source table SHA256 and the filtering rule. The two small tables are GTDB-derived reference excerpts: retain GTDB R226 attribution and its upstream terms. No new license or ownership claim is assigned here.

## Existing expected results

`expected/species/` and `expected/genus/` contain small outputs of the already accepted Stage2 formal-parameter run: FastSpar correlation/covariance/P-value matrices and selected edges, Zi-Pi node/edge tables, and EPI success summaries. Empty genus Zi-Pi tables and zero EPI candidates are legitimate successful scientific outcomes, not fabricated completion files. FastSpar edge counts are 25 species and 3 genus; the Zi-Pi selected graph has 6 species nodes/4 edges and 0 genus nodes/edges; EPI has 3 eligible species and 2 eligible genus groups, with zero candidates in each metric.

`expected/ARCHIVED_EXECUTION_SUMMARY.json` gives path-free accepted checks and historical runtimes. These are reference values, not evidence that the current checkout was executed on a different server. FastSpar, modularity and numerical libraries can vary across builds: comparisons must distinguish file identity, numerical tolerance and network/module equivalence rather than silently rewrite expected outputs.

No saved `qtrn`/`qtst`, checkpoint or DKI-scoring replay input is included. The complete abundance matrix supports fresh perturbation generation and a separate cNODE2 training/prediction workflow, but this bundle does not claim fresh-training acceptance or guarantee that a new stochastic training run produces identical historical candidates. A replay from saved predictions is a separate scope.

## Provenance, sharing and integrity

`SOURCE_PROVENANCE.tsv` records source collection-relative paths, scientific origin, transformations, file sizes and source/output SHA256. Source aliases are relative to the existing project collections, not machine-specific execution paths. `SHA256SUMS.tsv` covers every other file in this bundle. No patient identifiers or clinical attributes are bundled; the sample metadata contains public Run accessions only.

The author explicitly authorized preparation of this complete-project matrix example in the independent local candidate. Local example preparation and use can proceed. No upload or public release has occurred; the author must still determine the data license before public distribution. GTDB attribution and applicable upstream terms remain separate. The original files remain unchanged.
