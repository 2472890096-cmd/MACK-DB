# Script reference

This index covers all **90 active Python, R and Shell files** after the24November2026 structure revision. Entries include scientific implementations, required wrappers, validators and tests; they are not90 independent analysis steps. [The TSV](SCRIPT_REFERENCE.tsv) records inputs, outputs, CLI contracts, dependencies and source identities. Current file hashes are protected by the [release manifest](RELEASE_SCRIPT_MANIFEST.tsv).

Start from the [real matrix example](../examples/real_bioproject/PRJEB81106/README.md) or [paper reproduction](PAPER_REPRODUCTION.md). The current manuscript uses main Figures1–4, Supplementary FiguresS1/S2 and MethodsS1–S8. A historical Stage identifier does not name a manuscript panel.

Verification columns describe recorded execution scope; the dated current restructuring checks are reported separately. Full-data historical entries still require their explicit external inputs and held helpers.

## Profiling and abundance

| Script | Function | Paper/method scope |
|---|---|---|
| [pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py](../pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py) | Original species-label token genus aggregation; sums within sample and writes4decimals | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/abundance/count_sparcc.R](../pipeline/abundance/count_sparcc.R) | Multiply archived abundances by the accepted fixed scale for the FastSpar count input | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/abundance/shape_sparcc_format.py](../pipeline/abundance/shape_sparcc_format.py) | Merge per-Run Profile records into taxon-by-sample matrices with original prevalence/order behavior | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh](../pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh) | Apply the original adaptive species prevalence loop and delegate matrix generation | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/metakssd_profile/historical_linked_wrappers/query.sh](../pipeline/metakssd_profile/historical_linked_wrappers/query.sh) | Query a saved sample Sketch against a separately supplied compatible MarkerDB | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/metakssd_profile/profile.R](../pipeline/metakssd_profile/profile.R) | Filter query shared-k-mer support and calculate Poisson-corrected abundance percentages rounded to4decimals | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/metakssd_profile/run_single_run.py](../pipeline/metakssd_profile/run_single_run.py) | Bounded portable orchestration of fixed ERR14036375 pairedFASTQ to counted Sketch, pinned R226Query and originalRProfile; independent full-output acceptance | MethodS1 / database profiling, abundance and association workflow |

## Association inference and environment support

| Script | Function | Paper/method scope |
|---|---|---|
| [pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py](../pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py) | Apply the accepted final fastCCLasso absolute-correlation0.3 and P threshold to saved outputs | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/associations/fastcclasso/run_fastcclasso_parallel.R](../pipeline/associations/fastcclasso/run_fastcclasso_parallel.R) | Infer species associations using the separately supplied official fastCCLasso core and bootstrap protocol | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/associations/fastcclasso/run_fastcclasso_parallel_genus.R](../pipeline/associations/fastcclasso/run_fastcclasso_parallel_genus.R) | Infer genus associations with the same official core and rank-specific input/seed layout | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/associations/fastspar/corr_inter.R](../pipeline/associations/fastspar/corr_inter.R) | Convert selected correlation/P pairs to undirected project link files | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/associations/fastspar/run_fastspar.sh](../pipeline/associations/fastspar/run_fastspar.sh) | Run fixed main estimation,10000 bootstrap fits and10000 P-value permutations using external FastSpar | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/associations/fastspar/sparcc_corr_pvalue.R](../pipeline/associations/fastspar/sparcc_corr_pvalue.R) | Map correlation/P matrix IDs to taxon labels and select final significant pairs | MethodS1 / database profiling, abundance and association workflow |
| [pipeline/environment/fastcclasso/get_link_info_fastcclasso.py](../pipeline/environment/fastcclasso/get_link_info_fastcclasso.py) | Join environment support to exact fastCCLasso project coefficients/P/sample counts with fail-fast checks | MethodsS1,S5 / environment association inputs to Figure2 and literature validation |
| [pipeline/environment/fastcclasso/link_to_project_fastcclasso.py](../pipeline/environment/fastcclasso/link_to_project_fastcclasso.py) | Build environment pair membership and filter minimum independent BioProject support | MethodsS1,S5 / environment association inputs to Figure2 and literature validation |
| [pipeline/environment/get_link_info5.py](../pipeline/environment/get_link_info5.py) | Annotate FastSpar environment pairs with rounded project coefficients/P and sample counts; silent-error risks documented | MethodsS1,S5 / environment association inputs to Figure2 and literature validation |
| [pipeline/environment/link_to_project.py](../pipeline/environment/link_to_project.py) | Build project-support belonging files and optional minimum-support filtering | MethodsS1,S5 / environment association inputs to Figure2 and literature validation |

## Keystone methods

| Script | Function | Paper/method scope |
|---|---|---|
| [pipeline/dki/DKI_formula_fixed.py](../pipeline/dki/DKI_formula_fixed.py) | Train and predict cNODE2 compositions with accepted model/ODE implementation; not scoring replay | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/calculate_structural_keystoneness_rk4_fast.R](../pipeline/dki/calculate_structural_keystoneness_rk4_fast.R) | Calculate unified-R structural Ks from matched saved predictions and perturbations | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/generate_dki_input_from_abundance.py](../pipeline/dki/generate_dki_input_from_abundance.py) | Generate normalized Ptrain, focal-deletion Ptest and1-based sample/taxon perturbation indices | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/map_keystoneness_genus_names.py](../pipeline/dki/map_keystoneness_genus_names.py) | Map DKI genus score IDs to exact archived genus names | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/map_keystoneness_species_names.py](../pipeline/dki/map_keystoneness_species_names.py) | Map DKI species score IDs to exact archived species names | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/prepare_dki_keys_for_mysql.py](../pipeline/dki/prepare_dki_keys_for_mysql.py) | Export local named Top-K and mean+2sampleSD candidates; performs no database write | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/dki/r_scoring_revision/score_R_streaming.R](../pipeline/dki/r_scoring_revision/score_R_streaming.R) | Replay accepted R scoring in chunks after a validated legacy-reproduction gate; specialized audit adapter | MethodS4 / Figure3 DKI candidates; unified-R scoring |
| [pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R](../pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R) | Compute unchanged EPIv3 D1/D2/discrete-ECDF Q and per-metric candidate keys | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/genus/scripts/capture_v3.R](../pipeline/epi_v3/genus/scripts/capture_v3.R) | Source the pinned EPI core and save score/candidate objects for diagnostics | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/genus/scripts/genus_workflow.py](../pipeline/epi_v3/genus/scripts/genus_workflow.py) | Preserve genus production pilot, resource-aware scheduling, resume and verified final export orchestration | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/genus/scripts/prepare_inputs.py](../pipeline/epi_v3/genus/scripts/prepare_inputs.py) | Validate genus matrices and pin the unmodified species worker with exact production input/pilot gates | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/genus/scripts/run_one_project.sh](../pipeline/epi_v3/genus/scripts/run_one_project.sh) | Run a single isolated EPI worker with the preserved production argument contract | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/genus/scripts/timed_rscript.sh](../pipeline/epi_v3/genus/scripts/timed_rscript.sh) | Record R execution resource usage for the genus production adapter | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/scripts/workflow.py](../pipeline/epi_v3/scripts/workflow.py) | Shared EPI production worker, successful-output verification, resume and2-of-3 selection for both ranks | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/epi_v3/species/scripts/run_one_project.sh](../pipeline/epi_v3/species/scripts/run_one_project.sh) | Run a single isolated EPI worker with the preserved production argument contract | MethodS3 / Figure3 EPI candidates; shared mathematical implementation |
| [pipeline/zipi/FastSpar_signed_spinglass_network.R](../pipeline/zipi/FastSpar_signed_spinglass_network.R) | Build the filtered signed FastSpar network, optimize modules, calculate Zi/Pi and roles, then render | MethodS2 / Figure3 Zi-Pi candidates |
| [pipeline/zipi/extract_signed_zipi_keys_by_environment.py](../pipeline/zipi/extract_signed_zipi_keys_by_environment.py) | Extract Connector/Module hub/Network hub rows while preserving valid zero-candidate projects | MethodS2 / Figure3 Zi-Pi candidates |
| [pipeline/zipi/prepare_signed_keys_for_mysql.py](../pipeline/zipi/prepare_signed_keys_for_mysql.py) | Export local named Zi-Pi candidate tables; performs no database write | MethodS2 / Figure3 Zi-Pi candidates |

## Paper statistics

| Script | Function | Paper/method scope |
|---|---|---|
| [analysis/associations/original/plot_coabundance_panel_ABCDE_ARIAL.py](../analysis/associations/original/plot_coabundance_panel_ABCDE_ARIAL.py) | Figure2 selection and rendering | Figure2A–E / MethodS1 |
| [analysis/associations/original/plot_species_genus_coabundance_complete.py](../analysis/associations/original/plot_species_genus_coabundance_complete.py) | Figure2 statistics | Figure2A–E / MethodS1 |
| [analysis/associations/render_environment_figure.py](../analysis/associations/render_environment_figure.py) | Figure2可配置绘图入口；只适配字体/输出/summary路径 | Figure2A–E / MethodS1 |
| [analysis/associations/reproduce_environment_statistics.py](../analysis/associations/reproduce_environment_statistics.py) | Figure2可移植统计入口；核验原源码SHA再调用build_summaries/global_summary | Figure2A–E / MethodS1 |
| [analysis/associations/verify_environment_statistics.py](../analysis/associations/verify_environment_statistics.py) | 独立逐值核对90个环境/层级统计及最终面板选择 | Figure2A–E / MethodS1 |
| [analysis/keystone_concordance/concordance_release.py](../analysis/keystone_concordance/concordance_release.py) | Compute Broad/Strict support, exclusive method combinations and unique taxa after accepted four-project R replacement | Figure3A–D / MethodsS2–S4 / TableS3 support |
| [analysis/keystone_concordance/from_final_candidates.py](../analysis/keystone_concordance/from_final_candidates.py) | 从已经统一R的六候选表复算，检查61行修正已存在，禁止重复替换 | Figure3A–D / MethodsS2–S4 / TableS3 support |
| [analysis/paper_network_comparison/fixed_order_recovery.py](../analysis/paper_network_comparison/fixed_order_recovery.py) | 固定Chen→Wallen→Li→Guo→Lacruz→Zhou顺序累积恢复；独立集合复核 | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/original/complete_project/common.py](../analysis/paper_network_comparison/original/complete_project/common.py) | Complete-project helper | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/original/complete_project/compare_networks.py](../analysis/paper_network_comparison/original/complete_project/compare_networks.py) | Complete-project species statistics | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/original/genus_groups/common.py](../analysis/paper_network_comparison/original/genus_groups/common.py) | Genus-group helper | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/original/genus_groups/genus_overlap.py](../analysis/paper_network_comparison/original/genus_groups/genus_overlap.py) | Species-edge genus grouping | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/original/matched_subgroups/build_reference.py](../analysis/paper_network_comparison/original/matched_subgroups/build_reference.py) | Freeze GTDB R226 vocabulary and reference fingerprints for external comparisons | Figure4B / MethodS8; subgroup input, statistics or numerical validation as specified |
| [analysis/paper_network_comparison/original/matched_subgroups/independent_recalculation.py](../analysis/paper_network_comparison/original/matched_subgroups/independent_recalculation.py) | Independently verify historical subgroup overlap, signs, ranks and saved random graphs | Figure4B / MethodS8; subgroup input, statistics or numerical validation as specified |
| [analysis/paper_network_comparison/original/matched_subgroups/run_degree_null.py](../analysis/paper_network_comparison/original/matched_subgroups/run_degree_null.py) | Provide the rewire function used by final complete-project nulls; historical main uses older input scope | Figure4B / MethodS8; subgroup input, statistics or numerical validation as specified |
| [analysis/paper_network_comparison/original/matched_subgroups/stage5b2_analysis.py](../analysis/paper_network_comparison/original/matched_subgroups/stage5b2_analysis.py) | Restore numeric taxonomy, common-universe mappings and historical network comparison reference outputs | Figure4B / MethodS8; subgroup input, statistics or numerical validation as specified |
| [analysis/paper_network_comparison/run_original_statistics.py](../analysis/paper_network_comparison/run_original_statistics.py) | Figure4B/C/D共享便携入口；仅替换路径/状态，保留SHA固定原函数AST | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/test_io_adapter.py](../analysis/paper_network_comparison/test_io_adapter.py) | 适配层小测试：SHA、路径重定位幂等、输出隔离、缺失/重复archive拒绝、AST禁止顶层运行 | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_network_comparison/validate_statistics.py](../analysis/paper_network_comparison/validate_statistics.py) | 对A/B/C/D新统计与40个最终显示值/分母逐项比较 | Figure4A–D / MethodsS5,S8; panel-specific role as specified |
| [analysis/paper_support/replay.py](../analysis/paper_support/replay.py) | 从已存22000次null结果及pair/tested-set输入复算S1/S2和表汇总 | Supplementary FiguresS1/S2 / Table1 / TablesS2,S4 / MethodsS5,S7 |
| [analysis/paper_support/test_replay.py](../analysis/paper_support/test_replay.py) | 五项小测试检查原summary、源SHA拒绝、重复键/目录逃逸/已有输出拒绝 | Supplementary FiguresS1/S2 / Table1 / TablesS2,S4 / MethodsS5,S7 |

## Literature preparation and validation

| Script | Function | Paper/method scope |
|---|---|---|
| [analysis/external_validation/literature_networks/scripts/01_parse_external_networks.py](../analysis/external_validation/literature_networks/scripts/01_parse_external_networks.py) | Parse published Zhou, Guo and Lacruz edges and sample provenance without inventing taxonomy matches | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/scripts/02_harmonize_external_taxonomy.py](../analysis/external_validation/literature_networks/scripts/02_harmonize_external_taxonomy.py) | Resolve external labels through GTDB/NCBI and documented genome mappings; preserve ambiguity | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage2b/02b_taxonomy_rescue.py](../analysis/external_validation/literature_networks/stage2b/02b_taxonomy_rescue.py) | Rescue Zhou labels and Lacruz SGB mappings using recorded evidence | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage2c/li2024/scripts/02c_li2024_taxonomy.py](../analysis/external_validation/literature_networks/stage2c/li2024/scripts/02c_li2024_taxonomy.py) | Map Li2024 taxa, remove mapping self-loops and deduplicate undirected edges | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage3a/scripts/03a_baseline_overlap.py](../analysis/external_validation/literature_networks/stage3a/scripts/03a_baseline_overlap.py) | Compute baseline external-pair recovery and sign agreement | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_provenance.py](../analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_provenance.py) | Build study-specific source-project exclusion configuration from accession evidence | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_source_exclusion.py](../analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_source_exclusion.py) | Remove known source BioProjects and rebuild >=2-project support before recovery | MethodS5 taxonomy/source exclusion / TableS1; inputs to Table1 and Supplementary FigureS1 |
| [analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_null_models.py](../analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_null_models.py) | Compute random-pair and degree-preserving null recovery distributions | Supplementary FigureS1 / Table1 / TableS2 / MethodS5 |
| [analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_report.py](../analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_report.py) | Independently recalculate null moments, quantiles, empirical P and swap checks; write report | Supplementary FigureS1 / Table1 / TableS2 / MethodS5 |
| [analysis/external_validation/literature_networks/stage4a_recurrence/scripts/04a_build_literature_recurrence.py](../analysis/external_validation/literature_networks/stage4a_recurrence/scripts/04a_build_literature_recurrence.py) | Collapse multiple networks to publication-level unique pairs and study membership | Supplementary FigureS2 / TableS4 / MethodS7 |
| [analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_independent_validate.py](../analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_independent_validate.py) | Independently verify combined source exclusions and per-pair support | Supplementary FigureS2 / TableS4 / MethodS7 |
| [analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_recurrence_source_excluded_recovery.py](../analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_recurrence_source_excluded_recovery.py) | Exclude the union of a pair's reporting-study sources and recalculate recovery | Supplementary FigureS2 / TableS4 / MethodS7 |
| [analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_independent_recalculation.py](../analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_independent_recalculation.py) | Independently recalculate ranks, Spearman, rank tests and exact tables | Supplementary FigureS2 / TableS4 / MethodS7 |
| [analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py](../analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py) | Build eligible-project denominators and test recurrence against support count/fraction/rank | Supplementary FigureS2 / TableS4 / MethodS7 |
| [analysis/external_validation/shared_protocol/chen_wallen_null_models_final_v1.py](../analysis/external_validation/shared_protocol/chen_wallen_null_models_final_v1.py) | Provide the accepted encoded graph, rewire and null-summary protocol functions | Supplementary FigureS1 / Table1 / TableS2 / MethodS5 |
| [analysis/external_validation/stage5B/01_sample_rescue/audit_local_rescue.py](../analysis/external_validation/stage5B/01_sample_rescue/audit_local_rescue.py) | Match original Run/Profile values to local abundance columns used by complete-project comparisons | Figure4B / MethodS8; exact Run/subgroup input preparation |
| [analysis/external_validation/stage5B/01_sample_rescue/build_rescue.py](../analysis/external_validation/stage5B/01_sample_rescue/build_rescue.py) | Recover exact sample/Run membership and subgroup evidence for published projects | Figure4B / MethodS8; exact Run/subgroup input preparation |
| [analysis/external_validation/stage5B/01_sample_rescue/compile_rescue_metadata.py](../analysis/external_validation/stage5B/01_sample_rescue/compile_rescue_metadata.py) | Cross-check cached accession/BioSample/experiment membership and public metadata | Figure4B / MethodS8; exact Run/subgroup input preparation |
| [analysis/external_validation/stage5B/01_sample_rescue/rescue_io.py](../analysis/external_validation/stage5B/01_sample_rescue/rescue_io.py) | Shared evidence registration, path, TSV/JSON, checksum and opt-in cached retrieval helpers | Figure4B / MethodS8; exact Run/subgroup input preparation |
| [analysis/external_validation/stage5B/02_prepare_fastspar_input/finalize_stage5B1A.py](../analysis/external_validation/stage5B/02_prepare_fastspar_input/finalize_stage5B1A.py) | Independently read back prepared abundance/count matrices and freeze verified column identities | Figure4B / MethodS8; exact Run/subgroup input preparation |
| [analysis/external_validation/stage5B/02_prepare_fastspar_input/prepare_inputs.py](../analysis/external_validation/stage5B/02_prepare_fastspar_input/prepare_inputs.py) | Prepare exact Run-selected Profile/count inputs and mapping records; source of retained provenance | Figure4B / MethodS8; exact Run/subgroup input preparation |

## Shared runtime and tools

| Script | Function | Paper/method scope |
|---|---|---|
| [lib/mackdb_paths.py](../lib/mackdb_paths.py) | Resolve opt-in analysis data/code paths from safe JSON configuration; no scientific parameters | Required shared runtime, input validation or release-integrity support |
| [tools/check_dependencies.py](../tools/check_dependencies.py) | Probe installed Python/R packages and configured external tool paths without installation | Required shared runtime, input validation or release-integrity support |
| [tools/epi_project.py](../tools/epi_project.py) | Validate and freeze one project, invoke unchanged EPI worker, verify success and2-of3 outputs | Required shared runtime, input validation or release-integrity support |
| [tools/network_matrix_io.py](../tools/network_matrix_io.py) | Validate exact labelled matrix ID/order/NA contracts and adapt TSV network matrices to CSV | Required shared runtime, input validation or release-integrity support |
| [tools/verify_paper_inputs.py](../tools/verify_paper_inputs.py) | 按安全相对路径逐项校验已解压论文数据包大小和SHA | Required shared runtime, input validation or release-integrity support |
| [tools/verify_release_candidate.py](../tools/verify_release_candidate.py) | Verify current script/file hashes, immutable historical provenance, shared EPI layout, accepted scopes and input boundaries | Required shared runtime, input validation or release-integrity support |

## Examples

| Script | Function | Paper/method scope |
|---|---|---|
| [examples/real_bioproject/PRJEB81106/run_workflow.py](../examples/real_bioproject/PRJEB81106/run_workflow.py) | Portable complete33Run matrix workflow; formal unchanged scientific programs; fresh DKI training | Real-project demonstration; does not replace cohort-wide paper statistical inputs |

## Tests

| Script | Function | Paper/method scope |
|---|---|---|
| [tests/run_small_tests.py](../tests/run_small_tests.py) | Run nine bounded synthetic scientific test groups with explicit result logs | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |
| [tests/test_portable_project_workflow.py](../tests/test_portable_project_workflow.py) | Portable workflow I/O, branch dependency and failure rejection tests | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |
| [tests/test_single_run_workflow.py](../tests/test_single_run_workflow.py) | Small input/error/timeout/cleanup safety contracts for the portable upstream wrapper | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |
| [tests/unit/dki.R](../tests/unit/dki.R) | Check unified-R deletion baseline, Bray-Curtis, median, sampleSD/Top-K and legacy counterexample | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |
| [tests/unit/epi.R](../tests/unit/epi.R) | Check prevalence, D1/D2/Q, ECDF ties/fallback, strict threshold and legitimate zero candidates | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |
| [tests/unit/zipi.R](../tests/unit/zipi.R) | Check Zi/Pi values, sampleSD, boundary equality, valid NA and role selection | Synthetic scientific or workflow-boundary regression; not an additional manuscript result |

## Structure and history

- Both EPI ranks call one R core and one Python worker; rank-specific launchers select inputs and isolated runtime outputs.
- Figure4B subgroup source is under `analysis/paper_network_comparison/original/matched_subgroups/`. The different `stage4b_recurrence_recovery` task builds literature-support inputs for SupplementaryS2.
- Stage5A directed-complementarity, the04c recovery/logistic extension and three old Query/EPI diagnostics are retired from the default tree; original bytes remain in Git history. SupplementaryS1/S2 statistics and validators remain.
- Original84/93/151 inventories are immutable historical evidence, not active source counts. See [provenance](provenance/README.md).
