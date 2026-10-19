# Portable entry points and IO

Commands below use actual implemented arguments and run from the repository root. `PRJEXAMPLE` and `data/...` are input placeholders, not claims that real project data are bundled. The accepted nine-group synthetic suite was rerun in this private checkout. Full real-data Figure3 inputs are intentionally held pending upload authorization; prior v0.2 validation does not imply they are bundled here. Inference/training/full network commands are documented for later authorized runs, not executed here.

Use an existing compatible Python/R environment. Output directories should be new and separate from inputs. Historical batch launchers remain parameter references; do not use them as reviewer quick starts.

## Query, Profile, abundance and counts

```bash
mkdir -p work/query work/profile work/abundance
bash pipeline/metakssd_profile/historical_linked_wrappers/query.sh \
  data/reference/marker_db data/sketch/SAMPLE_sketch work/query/SAMPLE_all.txt
Rscript pipeline/metakssd_profile/profile.R \
  -l data/sample_ids.txt -i work/query -s 12 -o work/profile
bash pipeline/metakssd_profile/historical_linked_wrappers/batch_matrix.sh \
  data/project_ids.txt data/project_samples work/profile work/abundance
```

`data/project_samples/{project}_uniq.txt` contains sample IDs matching Profile filenames. The unchanged adaptive f loop delegates to `pipeline/abundance/shape_sparcc_format.py`, and writes `work/abundance/filtered_{f}/abundance|species`. The path migration does not introduce sorting or renumbering beyond the original dictionary enumeration. Historical Python ordering across versions is still a reproduction limitation.

```bash
Rscript pipeline/abundance/count_sparcc.R \
  -i data/species/abundance/PRJEXAMPLE_abundance.matrix \
  -n 1000000 -o work/PRJEXAMPLE_count.matrix
```

Genus abundance must already exist with a matching genus mapping. This repository does not invent a missing canonical species→genus generator.

## FastSpar and final postprocessing

This inference command performs the historical10000-bootstrap protocol; it is **not** the small test:

```bash
MACKDB_THREADS=1 bash pipeline/associations/fastspar/run_fastspar.sh \
  work/PRJEXAMPLE_count.matrix work/PRJEXAMPLE
```

Executables may be selected with `FASTSPAR_BIN`, `FASTSPAR_BOOTSTRAP_BIN`, `FASTSPAR_PVALUES_BIN`, `PARALLEL_BIN`; default lookup uses PATH. Resource concurrency is configurable, scientific iteration counts are fixed to the reviewed script. Successful completion removes only its own bootstrap directories. Existing directories are not silently overwritten by `mkdir`.

```bash
python tools/network_matrix_io.py \
  --corr work/PRJEXAMPLE_correlation.tsv --pvalue work/PRJEXAMPLE_pvalues.tsv \
  --mapping data/species/species/PRJEXAMPLE_species.matrix --output work/PRJEXAMPLE_csv
Rscript pipeline/associations/fastspar/sparcc_corr_pvalue.R \
  -c work/PRJEXAMPLE_csv/correlation.csv -p work/PRJEXAMPLE_csv/pvalues.csv \
  -s data/species/species/PRJEXAMPLE_species.matrix \
  --corr_num 0.3 --pvalue_num 0.0001 -o work/PRJEXAMPLE_corr_pvalue
Rscript pipeline/associations/fastspar/corr_inter.R \
  -f work/PRJEXAMPLE_corr_pvalue.csv -n work/PRJEXAMPLE
```

The source postprocessor distinguishes CSV matrices with headers from its historical tabular no-header branch. The new IO adapter accepts labelled TSV matrices only when row and column IDs exactly match the mapping order, checks upper-triangle NA-mask agreement, and emits CSV without changing numeric strings. It rejects mismatches instead of silently reordering/transposing. Synthetic execution compares direct CSV and adapted TSV results. An unexpected upstream matrix format is an input error, not permission to guess its axes.

## fastCCLasso

The external official core is required. Species example preserves the original initial correlation filter and explicit final refilter:

```bash
Rscript pipeline/associations/fastcclasso/run_fastcclasso_parallel.R \
  --project PRJEXAMPLE --abundance_dir data/species/abundance \
  --species_dir data/species/species --output work/fastcclasso/species \
  --fast_core external/fastCCLasso/RCode/fastCCLasso.R \
  --nboot 100 --kmax 20 --pvalue 0.0001 --cor_cutoff 0.1 \
  --workers 1 --seed 20260804 --write_full_matrices T --write_pair_csv T
python pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py \
  --level species --work-root work/fastcclasso/species
```

The refilter writes to the explicitly chosen working root in place, just as its original scientific postprocessing did; never point it at production results. Genus uses `run_fastcclasso_parallel_genus.R`, `--species_dir data/genus/genus`, `--abundance_dir data/genus/abundance`, seed20260806, genus output root and `--level genus`. No P-value formula is recreated in a wrapper. The configured default core is not bundled; use the official pinned retrieval and SHA check in EXTERNAL_DEPENDENCIES.md. Full inference was not executed in v0.2.

## Environment networks

```bash
python pipeline/environment/link_to_project.py \
  -l data/env/human_gut.list -d work/project_links -o work/environment \
  -n human_gut -f 2 -s T
python pipeline/environment/get_link_info5.py \
  -i work/environment/human_gut_filtered.belonging -o work/environment -n human_gut.info \
  -s data/projects.source -c work/corr_pv_csv -a data/species/abundance
```

For fastCCLasso use `pipeline/environment/fastcclasso/link_to_project_fastcclasso.py` with the same `-l/-d/-o/-n/-f 2/-s T` contract, then:

```bash
python pipeline/environment/fastcclasso/get_link_info_fastcclasso.py \
  -i work/environment/human_gut_filtered.belonging \
  -d work/fastcclasso/species/cor-pv_tsvformat -a data/species/abundance \
  -o work/environment -n human_gut_fastcclasso.info -e human_gut --rank species
```

Both support filters write `human_gut.belonging` (unfiltered) and `human_gut_filtered.belonging` (f=2). The attribute step must read the filtered file as above. The historical FastSpar filter console message mislabels the extension as .txt; the actual source writes .belonging. Do not infer a completed environment network from an empty/missing file alone.

## Zi-Pi

The full signed spin-glass/plotting entry retains the original computation. It now accepts data roots and taxonomy path:

```bash
Rscript pipeline/zipi/FastSpar_signed_spinglass_network.R \
  -p PRJEXAMPLE -l species -C data/networks/PRJEXAMPLE_corr.tsv \
  -P data/networks/PRJEXAMPLE_p.tsv --species-base data/species \
  --genus-base data/genus --taxonomy-file data/reference/bac120_taxonomy_r226_species_uniq.tsv \
  -f 0.3 --spins 25 --gamma 1 --gamma-minus 1 --seed 12345 \
  -a work/zipi/plots_a -b work/zipi/plots_b -c work/zipi/plots_c \
  -g work/zipi/gexf -d work/zipi/species_nodes
python pipeline/zipi/extract_signed_zipi_keys_by_environment.py \
  --level species --env-dir data/env --species-node-dir work/zipi/species_nodes \
  --species-out-dir work/zipi/species_keys --write-empty
python pipeline/zipi/prepare_signed_keys_for_mysql.py \
  --level species --species-dir work/zipi/species_keys --output-dir work/candidates
```

These last Python programs export local TSVs; they do not write MySQL. Full optimization and plots are not executed by small tests. A fallback basic node table after role failure is not equivalent to complete valid Zi/Pi output.

## EPI v3: recommended isolated single-project entry

```bash
python tools/epi_project.py --level species --project PRJEXAMPLE \
  --abundance data/species/abundance/PRJEXAMPLE_abundance.matrix \
  --mapping data/species/species/PRJEXAMPLE_species.matrix \
  --run-root work/epi_species/PRJEXAMPLE --rscript "$(command -v Rscript)"
```

For genus, change level, input paths and run root. The adapter checks matrix orientation/dimensions, exact ordered IDs/names and values. It copies the unchanged core and preserved worker into the new run root, freezes input/config hashes, then invokes `worker` and `postprocess`. It does not call the all-cohort `master`. Both ranks use exactly mode prevalence,.25–.75,pQ.1.

Outputs: three `keys/{project}_EPI-{D1,D2,Q}.key`, summary, status, raw logs, committed hashes and `final_2of3/{project}_EPI-2of3.key`. Zero-byte key files are accepted only with successful exit+summary+hash evidence. Resume the exact same command with `--resume`; it verifies the frozen setup and inputs and avoids rerunning successful projects. The production worker's per-file hard-link commits are unchanged; status/committed hashes determine complete success after an interruption.

The original workflow names remain in both rank paths. Their `ROOT` can be supplied by `MACKDB_EPI_RUN_ROOT`; config-relative core paths and PATH-resolved Rscript work without Dell paths. The separate full-genus adapter (`genus/scripts/prepare_inputs.py`) still requires the exact1575 formal list and historical pilot gates; it now initializes only an explicit isolated run root. It is not used as a toy-test shortcut and was not run on the full cohort.

For direct score diagnostics use the retained capture wrapper with `EPI_CAPTURE_CORE` set to the repository core and the same R CLI flags. `compare_saved_diagnostic_scores.R` now takes `REPRODUCED_ROOT AUDIT_DIR REFERENCE_VALIDATION_DIR`; its three named-project assertions are deliberately retained.

## DKI: accepted R scoring entry

Start from saved predictions if available; no training is required to rescore:

```bash
Rscript pipeline/dki/calculate_structural_keystoneness_rk4_fast.R \
  --project-list data/dki/projects.txt --data-base data/dki \
  --pred-base work/dki/species --threads 1
python pipeline/dki/map_keystoneness_species_names.py \
  --project-list data/dki/projects.txt --dki-base data/dki \
  --pred-base work/dki/species --species-dir data/species/species \
  --status-file work/dki/species/name_mapping_status.tsv
python pipeline/dki/prepare_dki_keys_for_mysql.py \
  --level species --species-root work/dki/species --env-dir data/env \
  --output-dir work/dki/candidates --sd-mode sample
```

`data/dki/{project}/data/` contains Ptrain/Ptest/species_id/sample_id CSVs; `work/dki/species/{project}/results/` contains qtrn/qtst. Check `reports/ks_rk4_fast_status.tsv`: `DONE`, `SKIP`, `NO_VALID_ROWS`, `MISSING_INPUT` and `FAILED` are distinct. The historical scorer can return process exit0 while reporting a failed project; its status table must be checked. This behavior is retained and tested where applicable, not concealed by the IO migration.

Genus mapping uses `map_keystoneness_genus_names.py --project-list ... --pred-base ... --genus-dir data/genus/genus --status-file ...`, followed by selection `--level genus --genus-root ...`. Neither export script writes a database.

`r_scoring_revision/score_R_streaming.R PROJECT DATA_DIR RESULT_DIR OUTPUT_DIR TOTAL_CASES` evaluates expressions from the adjacent accepted R scorer. It remains the specialized historical audit adapter and requires `OUTPUT_DIR/LOWMEM_REPRODUCTION.json` with PASS plus `taxon_name_map.tsv`. It is **not** a Lowmem/full-qtrn scoring implementation. Its gate must not be fabricated for real projects. The test uses a clearly labelled synthetic gate after checking the legacy counterexample. The regular R scorer above is the default published scoring entry.

For new DKI input, `generate_dki_input_from_abundance.py --project ... --abundance-dir ... --species-dir ... --out-base ...` implements the reviewed perturbation mapping. Training is intentionally excluded from quick tests. If separately authorized later, the actual accepted model command is:

```bash
python pipeline/dki/DKI_formula_fixed.py \
  --train-file data/dki/PRJEXAMPLE/data/Ptrain.csv \
  --test-file data/dki/PRJEXAMPLE/data/Ptest.csv \
  --output-dir work/dki/species/PRJEXAMPLE/results --project PRJEXAMPLE \
  --max-epochs 1000 --early-stop-patience 80 --lr 0.005 --batch-size 32 \
  --activation relu --bottleneck auto --normalize-z true --time-mode endpoint \
  --ode-method rk4 --ode-step-size 0.05 --device cpu --seed 20260629
```

Select the existing torch/torchdiffeq environment explicitly. The original default3000/200 does not describe the accepted production override1000/80.

## Figure3 — final code retained; real inputs held

The synthetic Figure3 test is part of `python tests/run_small_tests.py`. The final unified-R real-data program is retained unchanged, but full real candidate data is not committed. Availability and SHA references are in docs/FIGURE3_INPUT_AVAILABILITY.tsv. Do not run the following until separately authorized inputs have been supplied under the ignored data/ directory:

```bash
python analysis/figures/figure3/figure3_release.py \
  --base-analysis data/figure3/base_analysis \
  --scopes-dir frozen_inputs/figure3/scopes \
  --r-candidates data/figure3/r_replacement/FOUR_PROJECT_R_CANDIDATES_FINAL.tsv \
  --impact-audit data/figure3/impact_audit --output work/figure3
```

base_analysis must include inputs/key_{species,genus}{,_epi,_dki}.tsv plus FIGURE3_ALL_FOUR_SCOPES_22_VALUES.tsv; impact_audit must include FIGURE3_POTENTIAL_IMPACT.tsv. The accepted program checks the historical baseline, then replaces all four Species DKI projects with the61 accepted R candidates. Its intended final science is unified R, not mixed Lowmem. Do not pass already-replaced candidates as the historical baseline. Rendering is separate and was not run.

## Supporting analyses

`SUPPLEMENTARY_USAGE.tsv` gives the fixed45 roles and adoption state. The literature-network paths use `lib/mackdb_paths.py`; optionally set `MACKDB_PATHS_CONFIG` to a local JSON derived from `config/paths.example.json`. Same-directory helpers and selected legacy shared-protocol functions remain in the repository. Missing upstream taxonomy helper code, external datasets and frozen manifests are listed in DEPENDENCY_RESOLUTION.tsv. Their absence must stop full analysis; do not replace them with guessed data or change assertions. Stage5B and historical wrappers are not part of the default pipeline.
