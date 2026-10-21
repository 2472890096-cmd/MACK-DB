# EPI v3 — abundance-based D1, D2 and Q

Use `tools/epi_project.py --level species|genus --project ... --abundance ... --mapping ... --run-root NEW_ROOT --rscript ...` as the isolated entry. Input is a finite nonnegative taxon×sample TSV with exact ordered0-based mapping; no FastSpar/fastCCLasso/Zi-Pi input is required.

The unchanged shared core retains inclusive25–75% prevalence, focal deletion and renormalization, accepted D1/D2, Q pQ=0.1 from the unique off-diagonal sample-pair discrete ECDF with ties/fallback, strict score>mean+2sampleSD and final≥2of3. The core's raw default mode=all is not the formal worker override: the accepted entry explicitly uses prevalence. Never infer scientific defaults from CLI alone.

Outputs: three keys, summary, full success/status/log/hash evidence and2of3 keys. **PASS** is Stage2 species single-project execution and prior synthetic numerical tests; genus is **ARCHIVED_REPLAY** only. Zero candidates are valid with success evidence. Production full taxon scores were not independently available for the Stage2 formal-key comparison.

Compatibility copies are intentionally retained: core equals genus/scripts/core byte-for-byte; species/workflow equals genus/species_production_workflow. Do not merge them casually. Genus prepare_inputs/genus_workflow preserve exact1575/pilot/resource gates and are not the recommended small-example entry. Capture wrappers export diagnostics from the same core and do not redefine Q. Resume requires the frozen setup; no complete-cohort scheduling is launched by Quick Start.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.
