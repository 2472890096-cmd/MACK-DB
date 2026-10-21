# Signed-network Zi-Pi

Input: significant FastSpar corr/P matrices, the matching abundance/map and GTDB taxonomy. EPI/DKI are not prerequisites. `FastSpar_signed_spinglass_network.R` performs abundance>0.001 filtering, |r|≥0.3/P≤1e-4, isolate removal, signed negative spin-glass optimization (spins25,gamma1,gamma-minus1,seed12345), unweighted Zi/Pi and roles, then plots/exports.

Role boundaries use the accepted Zi2.5/Pi0.62 comparisons and sample SD; valid NA is not fabricated as a hub. Complete node tables contain z/p and taxa_roles. `extract_signed_zipi_keys_by_environment.py` selects Connectors/Module hubs/Network hubs; `prepare_signed_keys_for_mysql.py` writes local TSV only, not SQL.

**PASS**: Stage2 species full small network6nodes/4edges/3modules, numerical/NA/role comparison; genus archive replay yields0nodes legitimately. Exporter boundary tests are synthetic. A plot failure may occur after complete scientific output, but a fallback node table without valid core roles is not sufficient evidence. Preserve statuses/logs; do not decide success from nonempty keys alone.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.
