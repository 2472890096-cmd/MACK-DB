# BioProject abundance inputs

Species inputs are ordered `PROJECT_uniq.txt` Run lists and two-column Profile files. `shape_sparcc_format.py -l RUNS -i PROFILE_DIR -f 0.1 -p ABUNDANCE -s MAPPING` writes a taxon×sample TSV with leading ID column, `OTU_id` header and0-based sample columns, plus headerless0-based ID/name mapping. The adaptive shell wrapper increments f by0.1 while>2000taxa. The actual retained cutoff is Python `round(sample_count*f)`, not a new ceil rule; values are not renormalized here.

PRJEB81106 species53×33 values/order were independently reconstructed and read back in Stage2 (**ARCHIVED_REPLAY** provenance). The genus38×33 input uses the author-adopted original name-token aggregation: take the genus-prefix token after the species-label numeric prefix, sum matching rows, and format4decimals. This folds suffix genera such as Pseudomonas_E/Pseudomonas and Dorea_A/Dorea. All38labels/order/1254values match the specified converter; grouping fullGTDBnames would give40groups. See the [exact original-source contract](../../docs/GENUS_AGGREGATION_STATUS.md). This is the adopted historical grouping, not a strict fullGTDBgenus aggregation.

`count_sparcc.R -i ABUNDANCE -n 1000000 -o COUNT` supplies FastSpar's accepted scaled input; it is not an estimator of original read counts. Other branches consume abundance directly. Preserve exact names, including GTDB suffixes; do not transpose or trim names automatically.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.

The byte-identical [original converter](SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py) is now included. The default [complete-project runner](../../examples/real_bioproject/PRJEB81106/README.md) stages its legacy filenames and actually checks all1254genus cells; users need not locate a private converter path.
