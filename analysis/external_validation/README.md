# External literature validation and final-manuscript scope

This module preserves actual statistical computation and its independent checks, not a default plotting pipeline. Current full runs are **BLOCKED** by authorized data/helper availability; Stage3 performed no null simulations, full-project inference or scientific reanalysis. Previous accepted output evidence and static call relationships are separate from current runnable acceptance.

1. `literature_networks/scripts/01_parse_external_networks.py` parses original study files; `02_harmonize_external_taxonomy.py`, stage2b and stage2c apply documented GTDB/NCBI/SGB mappings.
2. stage3a computes baseline overlap/sign; stage3b builds source provenance and excludes source BioProject support before requiring≥2remaining projects. Its exact external `get_link_info_fastspar.py` dependency is held; `pipeline/environment/get_link_info5.py` is not a substitute.
3. stage3c calculates random-pair/degree-preserving nulls and verifies summaries; the shared Chen/Wallen protocol supplies imported/AST-extracted functions.
4. stage4a–d calculate publication recurrence, source-excluded recovery, logistic/exact tests, project support and ranks, with separate recalculators. They support final MethodsS5/S7 and relevant TablesS1/S2/S4–S6.
5. Retained Stage5B code reconstructs sample/column identity, count/reference provenance and the rewire function used by final complete-project comparisons. The older subgroup main is **not** final MethodsS9. Likewise the old stage5a builder supplies provenance, not final residual-denominator720-order statistics.

Final S8 and S9 main programs, Figure2 numerical source and exact upstream taxonomy/source helpers are among the14sources held pending per-file rights confirmation. They were not copied, replaced or advertised as runnable. S7/S8 are Methods, not nonexistent supplementary tables. TableS3 literature evidence remains author curation.

Configure repository-local input paths via MACKDB_PATHS_CONFIG and config/paths.example.json (paths resolve relative to repository root). Each script needs its own exact raw/frozen file layout and SHA evidence, not just one generic candidate TSV. Read the complete script-reference input/dependency columns and pending-source list. Most frozen-path programs have no argparse CLI: `--help` can be ignored and must not be used as a harmless probe. No generic command here claims that missing scientific inputs are supplied. External data need source-specific rights/access; public availability does not automatically grant redistribution.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.
