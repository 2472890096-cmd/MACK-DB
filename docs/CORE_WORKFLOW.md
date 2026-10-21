# Core scripts, functions and execution order

The database uses **264,767 successfully processed samples from 1,575 BioProjects**. Missing or unsuccessful source samples do not redefine this count. The primary keystone comparison uses the separately fixed 1,571-project list. Historic inventory totals are not new acceptance targets.

This table is a reading and execution guide, not a request to rerun the cohort. Commands are in [PIPELINE](../PIPELINE.md); the [complete index](SCRIPT_REFERENCE.tsv) describes every retained helper and CLI.

| Stage | Entry and necessary helpers | Function / input → output |
|---|---|---|
| FASTQ → Profile | [run_single_run.py](../pipeline/metakssd_profile/run_single_run.py), [profile.R](../pipeline/metakssd_profile/profile.R) | Paired FASTQ + configured MetaKSSD2.21/L3K11/R226 MarkerDB → counted Sketch, Query and Species Profile; sharekmer12 |
| Species matrix | [shape_sparcc_format.py](../pipeline/abundance/shape_sparcc_format.py), [batch_matrix.sh](../pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh) | Ordered Run list + Profiles → BioProject taxon×sample abundance and ID/name mapping |
| Genus matrix | [original converter](../pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py) | Extract genus-prefix token from species labels; sum rows per sample; write4decimals. See [exact source/CLI](GENUS_AGGREGATION_STATUS.md) |
| FastSpar input | [count_sparcc.R](../pipeline/abundance/count_sparcc.R) | Abundance → scaled count input, `-n1000000`; does not recover original read counts |
| FastSpar | [run_fastspar.sh](../pipeline/associations/fastspar/run_fastspar.sh), [network_matrix_io.py](../tools/network_matrix_io.py), [sparcc_corr_pvalue.R](../pipeline/associations/fastspar/sparcc_corr_pvalue.R), [corr_inter.R](../pipeline/associations/fastspar/corr_inter.R) | Count → correlation/covariance/P → final edges and project links; abs(r)>=.3,P<=1e-4 |
| fastCCLasso | [species wrapper](../pipeline/associations/fastcclasso/run_fastcclasso_parallel.R), [genus wrapper](../pipeline/associations/fastcclasso/run_fastcclasso_parallel_genus.R), [refilter](../pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py) | Abundance/map + external official core → r/P and final associations; kmax20,nboot100; actual P definition in [dependency guide](../EXTERNAL_DEPENDENCIES.md) |
| Zi-Pi | [signed network](../pipeline/zipi/FastSpar_signed_spinglass_network.R), [extract candidates](../pipeline/zipi/extract_signed_zipi_keys_by_environment.py), [prepare TSV](../pipeline/zipi/prepare_signed_keys_for_mysql.py) | FastSpar network + abundance/map/taxonomy → modules, Zi/Pi, roles and named candidates |
| EPI v3 | [epi_project.py](../tools/epi_project.py) → [shared R core](../pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R) and preserved workers | Abundance/map → D1/D2/Q keys and2-of-3 candidates; prevalence25–75%,pQ.1. Keep workers, not just the R file |
| DKI input/training | [input generator](../pipeline/dki/generate_dki_input_from_abundance.py) → [DKI_formula_fixed.py](../pipeline/dki/DKI_formula_fixed.py) | Abundance/map → perturbations → cNODE2 predictions qtrn/qtst; training is distinct from scoring replay |
| DKI unified-R score | [R scorer](../pipeline/dki/calculate_structural_keystoneness_rk4_fast.R), [species mapping](../pipeline/dki/map_keystoneness_species_names.py), [genus mapping](../pipeline/dki/map_keystoneness_genus_names.py), [selection](../pipeline/dki/prepare_dki_keys_for_mysql.py) | Matched perturbations/predictions → focal-deleted renormalized q_null → median Ks, ranks, thresholds and candidates; old full-qtrn Lowmem is not the default |
| Environment support | [link_to_project.py](../pipeline/environment/link_to_project.py), [get_link_info5.py](../pipeline/environment/get_link_info5.py), [fastCCLasso readers](../pipeline/environment/fastcclasso/) | Project links + environment lists + r/P → environment associations supported by>=2BioProjects, with coefficients and source IDs |
| Keystone concordance | [concordance_release.py](../analysis/keystone_concordance/concordance_release.py) | Six candidate tables + fixed scopes + accepted four-project R replacement → Broad/Strict combinations and unique taxa |
| Literature validation | [external-validation guide](../analysis/external_validation/README.md) | Authorized study networks + taxonomy/source mappings/support → overlap, source exclusion, null and recurrence statistics |

FastSpar and fastCCLasso are parallel branches. **Zi-Pi needs a network; EPI and DKI start directly from abundance.** Concordance and environment integration are later summaries, not inputs required by EPI/DKI.

The original genus converter is source-pinned, not silently replaced by a new implementation. Its legacy `_species.txt`/`_profile.txt` naming differs from production `.matrix` filenames; the [genus guide](GENUS_AGGREGATION_STATUS.md) explains filename-only staging and exact mapping order. The byte-identical source is included; the new matrix runner stages its filenames and verifies all genus cells against the archive.

## Existing execution evidence

- ERR14036375: real paired FASTQ→Sketch→Query→Profile, byte-identical final Profile.
- PRJEB81106: bundled complete33Run matrices; formal-parameter FastSpar/fastCCLasso, Zi-Pi/EPI and new CPU DKI training→unified-R scoring all passed at both ranks. See the [runnable guide](../examples/real_bioproject/PRJEB81106/README.md).
- Nine small scientific groups and twelve safety cases passed as recorded in [TESTING](../TESTING.md).

The optional historical SRA stream is not the tested default. A full-cohort rerun, clean installation on a second machine and complete manuscript re-execution are not claimed. Additional finalS8/S9/Figure2/helper sources remain listed in [PENDING_FINAL_STATISTICS_CODE.tsv](PENDING_FINAL_STATISTICS_CODE.tsv) pending individual release rights; no new methods are invented to fill those gaps.
