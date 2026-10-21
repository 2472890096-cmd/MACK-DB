# DKI: matrix training and unified-R scoring

For the complete project example, use [PRJEB81106 run_workflow.py](../../examples/real_bioproject/PRJEB81106/README.md). `--steps dki-score` generates perturbations, performs a new CPU cNODE2 training/prediction run, and calls the accepted unified-R scorer, taxon-name mapping and candidate selector. It does not use saved predictions. On2026-10-21 both ranks actually passed:53/38taxa×33samples, new training15.05/11.45seconds and3candidates each. [Stage execution records](../../tests/expected/PORTABLE_MATRIX_EXECUTION_20261021.tsv) distinguish real execution from environment installation recipes.

The scientific source and scoring release `MACKDB_DKI_R_20261012_v1` are unchanged. The remaining guide below describes historical interfaces and saved-prediction replay; its Stage2/3 statements are historical, not a limitation of the new runner. Keep the separate orchestration Python and model Python selected in the [environment guide](../../environments/portable-matrix/README.md).

## Historical interfaces: DKI — training versus unified-R scoring

Abundance/map inputs are independent of FastSpar and EPI. `generate_dki_input_from_abundance.py` creates normalized Ptrain (taxa×samples), focal-removal Ptest (taxa×perturbations) and1-based sample/taxon indices. `DKI_formula_fixed.py` is the retained cNODE2 training/prediction implementation, with existing upstream MIT notice. Formal overrides are1000epochs/patience80,lr0.005,batch32,relu,automaticbottleneck,normalizedZ,endpoint,rk4/step0.05,seed20260629; raw model CLI defaults3000/200 are not these production overrides. New training/prediction is **NOT_VERIFIED** in Stage2/3.

The default scoring entry is `calculate_structural_keystoneness_rk4_fast.R` with matched saved Ptrain/Ptest/IDs and qtrn/qtst. It uses focal-removed, renormalized q_null, Bray–Curtis scaled by(1-observed abundance), per-taxon median and complete-background sample SD. Map exact names with `map_keystoneness_*_names.py`; export candidates with `prepare_dki_keys_for_mysql.py --sd-mode sample`, applying TopK=min(N,max(3,ceil(0.01N))) and Ks≥mean+2sampleSD. These are local exports, not database writes.

**PASS, scoring replay only**: Stage2 reproduced all archived R medians/ranks, independent background thresholds and three candidates per rank. Genus remains archive-input replay. This consumes saved neural predictions; it neither trains nor predicts and cannot demonstrate new model convergence. Check project status DONE versus MISSING_INPUT/NO_VALID_ROWS/FAILED; process exit0 alone is insufficient.

The specialized streaming revision adapter retains its existing prior-legacy-reproduction PASS gate. Never fabricate that gate or make it the default training/scoring path. The old full-qtrn Lowmem formula is not the accepted scoring definition. No source formula, model parameter or scientific output was changed by documentation work.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.
