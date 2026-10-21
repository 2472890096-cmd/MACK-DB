# Script reference —82 entries

The [complete TSV](SCRIPT_REFERENCE.tsv) records every current Python/R/Shell entry, inputs/outputs, CLI, dependencies, source identity and execution scope. The default runnable entry is the [complete-project matrix example](../examples/real_bioproject/PRJEB81106/README.md); the [core guide](CORE_WORKFLOW.md) explains the scientific branches. Sources and source guard remain distinct from external installation recipes.

PASS describes a bounded stated route, not all cohorts or every historical branch. External final-results modules may remain BLOCKED for held inputs/helpers. Do not run a frozen-path program with `--help` unless it actually declares that CLI; some historical modules enter main regardless.

## abundance

| Script | Function | Status |
|---|---|---|
| [pipeline/abundance/count_sparcc.R](../pipeline/abundance/count_sparcc.R) | Multiply archived abundances by the accepted fixed scale for the FastSpar count input | PASS |
| [pipeline/abundance/shape_sparcc_format.py](../pipeline/abundance/shape_sparcc_format.py) | Merge per-Run Profile records into taxon-by-sample matrices with original prevalence/order behavior | ARCHIVED_REPLAY |
| [pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py](../pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py) | Original species-label token genus aggregation; sums within sample and writes4decimals | PASS |

## associations

| Script | Function | Status |
|---|---|---|
| [pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py](../pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py) | Apply the accepted final fastCCLasso absolute-correlation0.3 and P threshold to saved outputs | PASS |
| [pipeline/associations/fastcclasso/run_fastcclasso_parallel.R](../pipeline/associations/fastcclasso/run_fastcclasso_parallel.R) | Infer species associations using the separately supplied official fastCCLasso core and bootstrap protocol | PASS |
| [pipeline/associations/fastcclasso/run_fastcclasso_parallel_genus.R](../pipeline/associations/fastcclasso/run_fastcclasso_parallel_genus.R) | Infer genus associations with the same official core and rank-specific input/seed layout | ARCHIVED_REPLAY |
| [pipeline/associations/fastspar/corr_inter.R](../pipeline/associations/fastspar/corr_inter.R) | Convert selected correlation/P pairs to undirected project link files | PASS |
| [pipeline/associations/fastspar/run_fastspar.sh](../pipeline/associations/fastspar/run_fastspar.sh) | Run fixed main estimation,10000 bootstrap fits and10000 P-value permutations using external FastSpar | PASS |
| [pipeline/associations/fastspar/sparcc_corr_pvalue.R](../pipeline/associations/fastspar/sparcc_corr_pvalue.R) | Map correlation/P matrix IDs to taxon labels and select final significant pairs | PASS |

## dki

| Script | Function | Status |
|---|---|---|
| [pipeline/dki/DKI_formula_fixed.py](../pipeline/dki/DKI_formula_fixed.py) | Train and predict cNODE2 compositions with accepted model/ODE implementation; not scoring replay | PASS |
| [pipeline/dki/calculate_structural_keystoneness_rk4_fast.R](../pipeline/dki/calculate_structural_keystoneness_rk4_fast.R) | Calculate unified-R structural Ks from matched saved predictions and perturbations | PASS |
| [pipeline/dki/generate_dki_input_from_abundance.py](../pipeline/dki/generate_dki_input_from_abundance.py) | Generate normalized Ptrain, focal-deletion Ptest and1-based sample/taxon perturbation indices | PASS |
| [pipeline/dki/map_keystoneness_genus_names.py](../pipeline/dki/map_keystoneness_genus_names.py) | Map DKI genus score IDs to exact archived genus names | ARCHIVED_REPLAY |
| [pipeline/dki/map_keystoneness_species_names.py](../pipeline/dki/map_keystoneness_species_names.py) | Map DKI species score IDs to exact archived species names | PASS |
| [pipeline/dki/prepare_dki_keys_for_mysql.py](../pipeline/dki/prepare_dki_keys_for_mysql.py) | Export local named Top-K and mean+2sampleSD candidates; performs no database write | PASS |
| [pipeline/dki/r_scoring_revision/score_R_streaming.R](../pipeline/dki/r_scoring_revision/score_R_streaming.R) | Replay accepted R scoring in chunks after a validated legacy-reproduction gate; specialized audit adapter | NOT_VERIFIED |

## environment

| Script | Function | Status |
|---|---|---|
| [pipeline/environment/fastcclasso/get_link_info_fastcclasso.py](../pipeline/environment/fastcclasso/get_link_info_fastcclasso.py) | Join environment support to exact fastCCLasso project coefficients/P/sample counts with fail-fast checks | NOT_VERIFIED |
| [pipeline/environment/fastcclasso/link_to_project_fastcclasso.py](../pipeline/environment/fastcclasso/link_to_project_fastcclasso.py) | Build environment pair membership and filter minimum independent BioProject support | NOT_VERIFIED |
| [pipeline/environment/get_link_info5.py](../pipeline/environment/get_link_info5.py) | Annotate FastSpar environment pairs with rounded project coefficients/P and sample counts; silent-error risks documented | BLOCKED |
| [pipeline/environment/link_to_project.py](../pipeline/environment/link_to_project.py) | Build project-support belonging files and optional minimum-support filtering | NOT_VERIFIED |

## epi_v3

| Script | Function | Status |
|---|---|---|
| [pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R](../pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R) | Compute unchanged EPIv3 D1/D2/discrete-ECDF Q and per-metric candidate keys | PASS |
| [pipeline/epi_v3/genus/scripts/capture_v3.R](../pipeline/epi_v3/genus/scripts/capture_v3.R) | Source the pinned EPI core and save score/candidate objects for diagnostics | NOT_VERIFIED |
| [pipeline/epi_v3/genus/scripts/genus_workflow.py](../pipeline/epi_v3/genus/scripts/genus_workflow.py) | Preserve genus production pilot, resource-aware scheduling, resume and verified final export orchestration | NOT_VERIFIED |
| [pipeline/epi_v3/genus/scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R](../pipeline/epi_v3/genus/scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R) | Compute unchanged EPIv3 D1/D2/discrete-ECDF Q and per-metric candidate keys | ARCHIVED_REPLAY |
| [pipeline/epi_v3/genus/scripts/prepare_inputs.py](../pipeline/epi_v3/genus/scripts/prepare_inputs.py) | Validate genus matrices and pin the unmodified species worker with exact production input/pilot gates | NOT_VERIFIED |
| [pipeline/epi_v3/genus/scripts/run_one_project.sh](../pipeline/epi_v3/genus/scripts/run_one_project.sh) | Run a single isolated EPI worker with the preserved production argument contract | NOT_VERIFIED |
| [pipeline/epi_v3/genus/scripts/species_production_workflow.py](../pipeline/epi_v3/genus/scripts/species_production_workflow.py) | Preserve the species production worker/scheduler implementation as a genus compatibility dependency | NOT_VERIFIED |
| [pipeline/epi_v3/genus/scripts/timed_rscript.sh](../pipeline/epi_v3/genus/scripts/timed_rscript.sh) | Record R execution resource usage for the genus production adapter | NOT_VERIFIED |
| [pipeline/epi_v3/species/scripts/run_one_project.sh](../pipeline/epi_v3/species/scripts/run_one_project.sh) | Run a single isolated EPI worker with the preserved production argument contract | PASS |
| [pipeline/epi_v3/species/scripts/workflow.py](../pipeline/epi_v3/species/scripts/workflow.py) | Preserve EPI production worker, integrity/resume, status and2-of3 postprocessing | PASS |

## examples

| Script | Function | Status |
|---|---|---|
| [examples/real_bioproject/PRJEB81106/run_stage.py](../examples/real_bioproject/PRJEB81106/run_stage.py) | Orchestrate fixed PRJEB81106 matrix branches with input hashes, new outputs,600s cap and archive-only genus gate | PASS |

## external_validation

| Script | Function | Status |
|---|---|---|
| [analysis/external_validation/literature_networks/scripts/01_parse_external_networks.py](../analysis/external_validation/literature_networks/scripts/01_parse_external_networks.py) | Parse published Zhou, Guo and Lacruz edges and sample provenance without inventing taxonomy matches | BLOCKED |
| [analysis/external_validation/literature_networks/scripts/02_harmonize_external_taxonomy.py](../analysis/external_validation/literature_networks/scripts/02_harmonize_external_taxonomy.py) | Resolve external labels through GTDB/NCBI and documented genome mappings; preserve ambiguity | BLOCKED |
| [analysis/external_validation/literature_networks/stage2b/02b_taxonomy_rescue.py](../analysis/external_validation/literature_networks/stage2b/02b_taxonomy_rescue.py) | Rescue Zhou labels and Lacruz SGB mappings using recorded evidence | BLOCKED |
| [analysis/external_validation/literature_networks/stage2c/li2024/scripts/02c_li2024_taxonomy.py](../analysis/external_validation/literature_networks/stage2c/li2024/scripts/02c_li2024_taxonomy.py) | Map Li2024 taxa, remove mapping self-loops and deduplicate undirected edges | BLOCKED |
| [analysis/external_validation/literature_networks/stage3a/scripts/03a_baseline_overlap.py](../analysis/external_validation/literature_networks/stage3a/scripts/03a_baseline_overlap.py) | Compute baseline external-pair recovery and sign agreement | BLOCKED |
| [analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_provenance.py](../analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_provenance.py) | Build study-specific source-project exclusion configuration from accession evidence | BLOCKED |
| [analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_source_exclusion.py](../analysis/external_validation/literature_networks/stage3b_source_exclusion/scripts/03b_source_exclusion.py) | Remove known source BioProjects and rebuild >=2-project support before recovery | BLOCKED |
| [analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_null_models.py](../analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_null_models.py) | Compute random-pair and degree-preserving null recovery distributions | BLOCKED |
| [analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_report.py](../analysis/external_validation/literature_networks/stage3c_null_model/scripts/03c_report.py) | Independently recalculate null moments, quantiles, empirical P and swap checks; write report | BLOCKED |
| [analysis/external_validation/literature_networks/stage4a_recurrence/scripts/04a_build_literature_recurrence.py](../analysis/external_validation/literature_networks/stage4a_recurrence/scripts/04a_build_literature_recurrence.py) | Collapse multiple networks to publication-level unique pairs and study membership | BLOCKED |
| [analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_independent_validate.py](../analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_independent_validate.py) | Independently verify combined source exclusions and per-pair support | BLOCKED |
| [analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_recurrence_source_excluded_recovery.py](../analysis/external_validation/literature_networks/stage4b_recurrence_recovery/scripts/04b_recurrence_source_excluded_recovery.py) | Exclude the union of a pair's reporting-study sources and recalculate recovery | BLOCKED |
| [analysis/external_validation/literature_networks/stage4c_recurrence_statistics/scripts/04c_independent_recalculation.py](../analysis/external_validation/literature_networks/stage4c_recurrence_statistics/scripts/04c_independent_recalculation.py) | Independently recompute logistic and exact contingency statistics | BLOCKED |
| [analysis/external_validation/literature_networks/stage4c_recurrence_statistics/scripts/04c_recurrence_statistics.py](../analysis/external_validation/literature_networks/stage4c_recurrence_statistics/scripts/04c_recurrence_statistics.py) | Fit recurrence/recovery logistic models, Fisher/trend tests and BH correction | BLOCKED |
| [analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_independent_recalculation.py](../analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_independent_recalculation.py) | Independently recalculate ranks, Spearman, rank tests and exact tables | BLOCKED |
| [analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py](../analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py) | Build eligible-project denominators and test recurrence against support count/fraction/rank | BLOCKED |
| [analysis/external_validation/literature_networks/stage5a_cross_study_complementarity/logs/build_stage5a.py](../analysis/external_validation/literature_networks/stage5a_cross_study_complementarity/logs/build_stage5a.py) | Preserve historical study-complementarity processing rules and provenance used by the held final720-order builder | BLOCKED |
| [analysis/external_validation/shared_protocol/chen_wallen_null_models_final_v1.py](../analysis/external_validation/shared_protocol/chen_wallen_null_models_final_v1.py) | Provide the accepted encoded graph, rewire and null-summary protocol functions | BLOCKED |
| [analysis/external_validation/stage5B/01_sample_rescue/audit_local_rescue.py](../analysis/external_validation/stage5B/01_sample_rescue/audit_local_rescue.py) | Match original Run/Profile values to local abundance columns used by complete-project comparisons | BLOCKED |
| [analysis/external_validation/stage5B/01_sample_rescue/build_rescue.py](../analysis/external_validation/stage5B/01_sample_rescue/build_rescue.py) | Recover exact sample/Run membership and subgroup evidence for published projects | BLOCKED |
| [analysis/external_validation/stage5B/01_sample_rescue/compile_rescue_metadata.py](../analysis/external_validation/stage5B/01_sample_rescue/compile_rescue_metadata.py) | Cross-check cached accession/BioSample/experiment membership and public metadata | BLOCKED |
| [analysis/external_validation/stage5B/01_sample_rescue/rescue_io.py](../analysis/external_validation/stage5B/01_sample_rescue/rescue_io.py) | Shared evidence registration, path, TSV/JSON, checksum and opt-in cached retrieval helpers | BLOCKED |
| [analysis/external_validation/stage5B/02_prepare_fastspar_input/finalize_stage5B1A.py](../analysis/external_validation/stage5B/02_prepare_fastspar_input/finalize_stage5B1A.py) | Independently read back prepared abundance/count matrices and freeze verified column identities | BLOCKED |
| [analysis/external_validation/stage5B/02_prepare_fastspar_input/prepare_inputs.py](../analysis/external_validation/stage5B/02_prepare_fastspar_input/prepare_inputs.py) | Prepare exact Run-selected Profile/count inputs and mapping records; source of retained provenance | BLOCKED |
| [analysis/external_validation/stage5B/04_network_comparison/logs/build_reference.py](../analysis/external_validation/stage5B/04_network_comparison/logs/build_reference.py) | Freeze GTDB R226 vocabulary and reference fingerprints for external comparisons | BLOCKED |
| [analysis/external_validation/stage5B/04_network_comparison/logs/independent_recalculation.py](../analysis/external_validation/stage5B/04_network_comparison/logs/independent_recalculation.py) | Independently verify historical subgroup overlap, signs, ranks and saved random graphs | BLOCKED |
| [analysis/external_validation/stage5B/04_network_comparison/logs/run_degree_null.py](../analysis/external_validation/stage5B/04_network_comparison/logs/run_degree_null.py) | Provide the rewire function used by final complete-project nulls; historical main uses older input scope | BLOCKED |
| [analysis/external_validation/stage5B/04_network_comparison/logs/stage5b2_analysis.py](../analysis/external_validation/stage5B/04_network_comparison/logs/stage5b2_analysis.py) | Restore numeric taxonomy, common-universe mappings and historical network comparison reference outputs | BLOCKED |

## keystone_concordance

| Script | Function | Status |
|---|---|---|
| [analysis/keystone_concordance/concordance_release.py](../analysis/keystone_concordance/concordance_release.py) | Compute Broad/Strict support, exclusive method combinations and unique taxa after accepted four-project R replacement | BLOCKED |

## lib

| Script | Function | Status |
|---|---|---|
| [lib/mackdb_paths.py](../lib/mackdb_paths.py) | Resolve opt-in analysis data/code paths from safe JSON configuration; no scientific parameters | NOT_VERIFIED |

## metakssd_profile

| Script | Function | Status |
|---|---|---|
| [pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh](../pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh) | Apply the original adaptive species prevalence loop and delegate matrix generation | NOT_VERIFIED |
| [pipeline/metakssd_profile/historical_linked_wrappers/query.sh](../pipeline/metakssd_profile/historical_linked_wrappers/query.sh) | Query a saved sample Sketch against a separately supplied compatible MarkerDB | ARCHIVED_REPLAY |
| [pipeline/metakssd_profile/profile.R](../pipeline/metakssd_profile/profile.R) | Filter query shared-k-mer support and calculate Poisson-corrected abundance percentages rounded to4decimals | ARCHIVED_REPLAY |
| [pipeline/metakssd_profile/run_single_run.py](../pipeline/metakssd_profile/run_single_run.py) | Bounded portable orchestration of fixed ERR14036375 pairedFASTQ to counted Sketch, pinned R226Query and originalRProfile; independent full-output acceptance | PASS |

## portable_project

| Script | Function | Status |
|---|---|---|
| [examples/real_bioproject/PRJEB81106/run_workflow.py](../examples/real_bioproject/PRJEB81106/run_workflow.py) | Portable complete33Run matrix workflow; formal unchanged scientific programs; fresh DKI training | PASS |

## test

| Script | Function | Status |
|---|---|---|
| [tests/test_portable_project_workflow.py](../tests/test_portable_project_workflow.py) | Portable workflow I/O, branch dependency and failure rejection tests | PASS |

## tests

| Script | Function | Status |
|---|---|---|
| [tests/epi_v3/scripts_used/capture_v3_scores.R](../tests/epi_v3/scripts_used/capture_v3_scores.R) | Capture EPI core scores for the independent numerical diagnostics | ARCHIVED_REPLAY |
| [tests/epi_v3/scripts_used/compare_saved_diagnostic_scores.R](../tests/epi_v3/scripts_used/compare_saved_diagnostic_scores.R) | Compare captured three-project EPIv3 scores and candidates to frozen references | ARCHIVED_REPLAY |
| [tests/query/query_reproduction.sh](../tests/query/query_reproduction.sh) | Run the bounded saved-Sketch query reproduction helper; no FASTQ generation | ARCHIVED_REPLAY |
| [tests/run_small_tests.py](../tests/run_small_tests.py) | Run nine bounded synthetic scientific test groups with explicit result logs | PASS |
| [tests/unit/dki.R](../tests/unit/dki.R) | Check unified-R deletion baseline, Bray-Curtis, median, sampleSD/Top-K and legacy counterexample | PASS |
| [tests/unit/epi.R](../tests/unit/epi.R) | Check prevalence, D1/D2/Q, ECDF ties/fallback, strict threshold and legitimate zero candidates | PASS |
| [tests/unit/zipi.R](../tests/unit/zipi.R) | Check Zi/Pi values, sampleSD, boundary equality, valid NA and role selection | PASS |
| [tests/test_single_run_workflow.py](../tests/test_single_run_workflow.py) | Small input/error/timeout/cleanup safety contracts for the portable upstream wrapper | PASS |

## tools

| Script | Function | Status |
|---|---|---|
| [tools/check_dependencies.py](../tools/check_dependencies.py) | Probe installed Python/R packages and configured external tool paths without installation | NOT_VERIFIED |
| [tools/epi_project.py](../tools/epi_project.py) | Validate and freeze one project, invoke unchanged EPI worker, verify success and2-of3 outputs | PASS |
| [tools/network_matrix_io.py](../tools/network_matrix_io.py) | Validate exact labelled matrix ID/order/NA contracts and adapt TSV network matrices to CSV | PASS |
| [tools/verify_release_candidate.py](../tools/verify_release_candidate.py) | Verify all77 script hashes, historical manifests, compatibility copies and approved data boundaries | PASS |

## zipi

| Script | Function | Status |
|---|---|---|
| [pipeline/zipi/FastSpar_signed_spinglass_network.R](../pipeline/zipi/FastSpar_signed_spinglass_network.R) | Build the filtered signed FastSpar network, optimize modules, calculate Zi/Pi and roles, then render | PASS |
| [pipeline/zipi/extract_signed_zipi_keys_by_environment.py](../pipeline/zipi/extract_signed_zipi_keys_by_environment.py) | Extract Connector/Module hub/Network hub rows while preserving valid zero-candidate projects | NOT_VERIFIED |
| [pipeline/zipi/prepare_signed_keys_for_mysql.py](../pipeline/zipi/prepare_signed_keys_for_mysql.py) | Export local named Zi-Pi candidate tables; performs no database write | NOT_VERIFIED |
