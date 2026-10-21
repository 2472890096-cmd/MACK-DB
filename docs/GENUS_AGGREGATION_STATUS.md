# Adopted species-name-based genus aggregation

The author confirmed on 21 October 2026 that documentation should follow the original `SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py`. Historical genus matrices remain the adopted inputs; no alternative grouping or downstream recalculation is requested.

The second tab-separated mapping field is parsed exactly as:

```python
genus = label.split(' ')[0].split('_')[1]
```

For example, `11_Escherichia_coli` maps to `Escherichia`, `374_Bacillus_A_cereus` to `Bacillus`, and labels with `Pseudomonas_E` or `Dorea_A` to `Pseudomonas` or `Dorea`. Rows sharing the extracted token are summed separately for each sample; output uses `.4f`. This source performs no taxonomy lookup, extra normalization or independent genus prevalence filter. Its genus IDs follow first encounter in mapping order and sample columns remain unchanged.

This **species-label-derived genus-prefix grouping** folds some GTDB genus suffixes. It is not strict grouping by the complete GTDB `g__` field. The description follows actual source behavior, not a renamed scientific algorithm.

## Original CLI and filenames

`-i/--input` is one BioProject per line; `-d/--dir` must contain `{project}_species.txt` and `{project}_profile.txt`. The species mapping is headerless ID/name; its row position is used as the species ID, so it must be ordered0..N-1 to match the matrix. Here `_profile.txt` is the **BioProject abundance matrix with a header**, not a single Run's two-column Profile. `-o/--outdir` must already contain `genus/` and `abundance/`; the converter writes `{project}_genus.txt` and `{project}_abundance.txt`. `-t/--tool` is parsed but unused. Dependencies imported by the original source include pandas/numpy.

Use filename-only staging for existing `.matrix` inputs and use the included byte-identical original source:

```bash
GENUS_CONVERTER=pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py
mkdir -p work/genus_input work/genus_output/genus work/genus_output/abundance
cp -n data/species/species/PRJEXAMPLE_species.matrix work/genus_input/PRJEXAMPLE_species.txt
cp -n data/species/abundance/PRJEXAMPLE_abundance.matrix work/genus_input/PRJEXAMPLE_profile.txt
printf '%s\n' PRJEXAMPLE > work/genus_projects.txt
python "$GENUS_CONVERTER" -i work/genus_projects.txt -d work/genus_input -o work/genus_output
```

The portable complete-project runner implements this filename-only staging and actually checks every genus abundance cell against the bundled archive. Always use fresh output locations: the original script directly opens output files for writing and can overwrite matching names.

The specified source already reproduced all 38×33archived PRJEB81106 values, labels and order. Strict complete-GTDB-name grouping would produce40groups; that difference is explained by the adopted suffix folding, not by parameter tuning. Only the old `run_stage.py` archival entry retains `--allow-archived-genus-replay`. The new default `run_workflow.py` validates the bundle and runs the original converter directly; it does not introduce a second genus definition.

The source SHA is in [GENUS_CONVERTER_SOURCE.tsv](GENUS_CONVERTER_SOURCE.tsv); its byte-identical file is included in the current 82-script candidate. Matching output is distinct from proving the exact original execution for every historical BioProject. No original source, matrix, paper or result was modified.

See [core workflow](CORE_WORKFLOW.md), [matrix guide](../pipeline/abundance/README.md) and [PRJEB81106 example](../examples/real_bioproject/PRJEB81106/README.md).
