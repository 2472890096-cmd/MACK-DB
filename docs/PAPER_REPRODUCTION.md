# Reproduce the paper

This guide is bound to the author-confirmed final `MACKDB_Manuscript.docx` (the requested manuscript9), supplementary information8 and revised workbook. Their [identities](paper/MANUSCRIPT_BINDING.tsv) were recorded on12November2026. Current scope: Figures1–4/Table1, MethodsS1–S8, FiguresS1/S2 and TablesS1–S4. The manuscript's older code reference is9f4fce7; this branch has additional analysis adapters and evidence. Record the actual checkout with `git rev-parse HEAD`.

Start from the **fixed processed inputs actually used for the paper**. Neither raw FASTQ nor fresh whole-cohort inference is required for these statistical checks. The separate [real matrix example](../examples/real_bioproject/PRJEB81106/README.md) demonstrates inference, including fresh DKI training. This guide does not substitute saved predictions for training.

## Access and preparation

[DATA_PACKAGES](paper/DATA_PACKAGES.md) lists local prepared archives, sizes, SHA256 and missing download identifiers. [The per-file manifest](paper/DATA_PACKAGE_MANIFEST.tsv) records logical paths and access holds. Obtain the approved data and any [held original sources](paper/SOURCE_ACCESS.tsv) before running the relevant commands. **There is no public URL for those packages yet.** A path on the author's server is not a download location. Sources/fonts are not hidden inside data packages.

Use Python3.11+ with NumPy, pandas, SciPy, Matplotlib, Pillow, openpyxl and statsmodels for the relevant modules; the exact tested environment and module dependencies are documented alongside each entry. Set `MACKDB_PAPER_PYTHON` to that interpreter. The existing [environment instructions](../environments/portable-matrix/README.md) retain SciPy1.15.3 from PyPI, version-object R checks and the fixed fastCCLasso core. No installer ran during this task.

Before using paper commands, select the installed analysis interpreter and check its dependencies (this check installs nothing):

```bash
export MACKDB_PAPER_PYTHON="$MACKDB_PYTHON"
"$MACKDB_PAPER_PYTHON" -c 'import sys; assert sys.version_info >= (3,11); import numpy,pandas,scipy,matplotlib,PIL,openpyxl,statsmodels; print(sys.executable,sys.version)'
```

Use another explicitly configured interpreter if the matrix orchestrator environment lacks the paper-only packages. The optional DKI training interpreter is not the paper interpreter.

Extract approved packages below ignored `data/`. Verify them before computation, for example:

```bash
"$MACKDB_PAPER_PYTHON" tools/verify_paper_inputs.py   --manifest analysis/associations/INPUT_DATA_MANIFEST.tsv --root data
"$MACKDB_PAPER_PYTHON" tools/verify_paper_inputs.py   --manifest docs/paper/FIGURE3_INPUT_MANIFEST.tsv --root data
"$MACKDB_PAPER_PYTHON" tools/verify_paper_inputs.py   --manifest analysis/paper_support/INPUT_ARCHIVE_MANIFEST.tsv --root data/paper-inputs
```

Every analysis uses a **new output directory**, stops on validation errors and records input/source hashes. Successful execution is separate from comparison with archived results.

## Figure 2: environment associations

**Input:** all45 environments × species/genus × FastSpar/fastCCLasso,180 `*.info` files. Required columns include taxon endpoints, `aver_corr`, `aver_pv`, `proj_count`, `proj_ID`, `proj_corr`, `proj_pv` and `proj_sam`. These are processed associations with method-specific project support, not raw sequencing data.

**Source:** the exact `plot_species_genus_coabundance_complete.py` is separately supplied and SHA-pinned. The adapter calls its unchanged statistical functions. Configure `FIGURE2_STATISTICS_SOURCE` and the final renderer `FIGURE2_SELECTION_SOURCE`; see [module instructions](../analysis/associations/README.md) for source hashes, dependencies and the complete verification/rendering commands.

```bash
"$MACKDB_PAPER_PYTHON" analysis/associations/reproduce_environment_statistics.py   --statistics-source "$FIGURE2_STATISTICS_SOURCE"   --species-fastspar data/figure2/species/FastSpar/f2   --genus-fastspar data/figure2/genus/FastSpar/f2   --species-cclasso data/figure2/species/CCLasso/f2   --genus-cclasso data/figure2/genus/CCLasso/f2   --output work/paper/figure2-statistics

"$MACKDB_PAPER_PYTHON" analysis/associations/verify_environment_statistics.py   --summary work/paper/figure2-statistics/environment_method_summary.tsv   --reference-summary data/figure2/reference/environment_method_summary.tsv   --plot-values data/figure2/reference/panel_plot_values.tsv   --selection-source "$FIGURE2_SELECTION_SOURCE"   --output work/paper/figure2-validation
```

**Outputs:**90-row environment summary, support-count distribution, global totals, input hashes, and all numerical panel comparisons. A/B use the same top15 environments by **combined species+genus union**. C uses the combined-union selection; D/E select by combined common-pair counts and order by mean sign agreement. Selection is taken from the actual final renderer, not invented independently.

**2E definition:** each method uses the unweighted signed mean across its own supporting BioProjects. The two support sets may differ. Shared unordered pairs supply one saved mean from each method. Spearman drops missing means, **retains zero**, and returns missing for fewer than3 pairs or a constant vector. Sign agreement separately removes zero means; it must not provide the Spearman mask. The replay starts from saved `aver_corr`; an independent scan checked all8,594,529 project lists/means. Maximum mean discrepancies from stored project coefficients were below8.34e-11, consistent with stored precision. It did not newly estimate project correlations. Do not round rho before comparison.

**Checks:** species(FastSpar,fastCCLasso,common,union)=(4661190,2287364,2059581,4888973); genus=(999851,646124,519093,1126882). Weighted Jaccard totals round to0.421/0.461. All90 rows reproduce byte-for-byte; all30 E points and all120 plot records match;13 missing rho values remain missing.

**Drawing:** the separately supplied final renderer can be run through `analysis/associations/render_environment_figure.py --renderer-source ... --summary ... --font-dir ... --output ...`; the [full command](../analysis/associations/README.md#redraw-from-a-saved-summary) uses actual implemented arguments. Four Arial font files remain external. This task actually redrew Figure2 in isolation and regenerated its plot-value table byte-for-byte. Images/fonts are not committed.

## Figure 3: candidate-method concordance

**Input:** six final candidate TSVs with `name`, `bioproject`, `environment`; accepted61-row Species DKI R reference; three project scope lists; expected aggregate tables. No production database query is used. The main1,571-project list includes valid zero-candidate projects. It must not be inferred from whether a candidate row exists. `PRJEB51353` is absent from all three scopes.

The preferred entry accepts an **already corrected final snapshot** and asserts its four-project membership without replacing anything:

```bash
"$MACKDB_PAPER_PYTHON" analysis/keystone_concordance/from_final_candidates.py   --candidates-dir data/figure3/final_candidates   --r-candidates data/figure3/r_replacement/FOUR_PROJECT_R_CANDIDATES_FINAL.tsv   --scopes-dir frozen_inputs/figure3/scopes   --expected-dir frozen_inputs/figure3/expected   --output work/paper/figure3
```

It calls the unchanged accepted `concordance_release.calc()`. The historical `concordance_release.py` remains available for pre-correction snapshots; its explicit61-row replacement and historical baseline assertions must not be applied to already-final inputs. Both paths were executed and agreed in this task; the new entry rejects the old snapshot.

**Outputs:**66 values (main22 plus two sensitivity22s), every Broad/Strict combination and sources, full Strict All3 support and per-project candidate coverage. Broad merges method evidence within an environment; Strict requires the same taxon in the same project. At-least-two includes All3. Candidate coverage is not relabeled as independent scientific-success evidence.

| Panel | Final main values |
|---|---|
| A species, Z only/E only/D only/ZE/ZD/ED/All3 |4575,3604,4786,879,1141,356,434|
| B genus, same order |2108,1260,1826,330,621,161,273|
| C Broad species, ≥2/All3 unique |2491,416|
| C Broad genus |982,243|
| D Strict species |970,25|
| D Strict genus |399,19|

**Checks:**66/66 aggregate values; all138 Strict All3 support rows across three scopes; main25 species/19 genus groups (27/19 project–taxon supports); human-gut9 species/8 genera (11/8 supports). All match accepted references. The old-mixed and final Species DKI snapshots differ in the intended31 removed/31 added taxa; no scoring is repeated.

**Artwork boundary:** the final Arial Figure3 script updates an existing registered SVG. It does **not** calculate UpSet intersections. Its existing SVG,22-value registration, protected files, renderer and fonts remain a separate legacy rendering dependency. Statistics above were rerun; that SVG chain was inspected but not executed or advertised as a TSV-to-UpSet generator. See [source access](paper/SOURCE_ACCESS.tsv).

## Figure 4: fixed-order and corresponding-project comparisons

See [the runnable module guide](../analysis/paper_network_comparison/README.md) for the exact A/B/C/D commands, JSON template, original-source SHA pins, path rebasing, input layout and40-value validator. Existing source directories still say `Figure7`. No new FastSpar inference or simulations are launched by these entries.

| Panel | Processed starting input and entry | Check |
|---|---|---|
| A | Six mapped publication pair sets, per-pair original project-support lookup,9-source exclusion; `fixed_order_recovery.py` | Fixed Chen,Wallen,Li,Guo,Lacruz-Pleguezuelos,Zhou order; cumulative2045,3284,3550,3712,4805,6067; final6067/11900 |
| B | Five groups' saved FastSpar matrices, exact Run/species mappings and literature tables; `run_original_statistics.py --mode b` | Guo Common species31/130; genus groups28/101 |
| C | Complete BioProject matrices, existing FastSpar matrices, mapped published sets/common taxa; `--mode c` | Guo shared38/201; Wallen523/1413 and185/1413 |
| D | Each side's complete saved species-edge set and production-label grouping; `--mode d` | Guo shared34/149; Wallen307/796 and118/796 |

Copy the real configuration template before editing local roots:

```bash
cp -n config/paper_network_comparison.example.json config/paper_network_comparison.local.json
```

After the four computations described in that guide:

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_network_comparison/validate_statistics.py   --a work/paper/Figure4A --b work/paper/Figure4B   --c work/paper/Figure4C --d work/paper/Figure4D   --output work/paper/Figure4_validation
```

**A:** nine exclusions are PRJNA319574,PRJNA942468,PRJNA834801,PRJEB32762,PRJNA751448,PRJEB42304,PRJEB7774,PRJEB11419,PRJEB41311; require≥2 remaining support projects. The independent fixed-order check validates exclusions for the processed literature-pair lookup. It does not rescan the whole database to independently reconstruct the archived full D0 size652972. MHO is excluded from the primary publication union. All six steps match historical order22. The old720-order analysis remains history, not the current default computation.

**B:** the public [Figure4B_sequencing_runs.tsv](../analysis/paper_network_comparison/reference/Figure4B_sequencing_runs.tsv) records91 Guo Control,106 MECFS,197 pooled Common,490 Wallen PD and234 NHC Runs:1118 group–Run rows,921 distinct Runs. Common pools runs before FastSpar. MACK-DB Guo uses `abs(r)>=0.3 && P<=0.0001`; Wallen uses `abs(r)>0.2 && P<0.05`. Literature-side filters and mapped-common-species denominators are recorded in the original functions/module guide. Complete-project Wallen has725 Runs; it is not a newly pooled724-Run B inference.

**C/D:** complete-project conditions distinguish Guo aggregate, Wallen published-matched and Wallen MACK-DB-standard (`abs(r)>=0.3 && P<=0.0001`) thresholds. The orange full-project network-size bar is not the overlap denominator. D removes within-group loops and deduplicates unordered genus-group pairs separately on both side-sets, then intersects. It uses the frozen production-label token rule, not a fresh genus network or conversion of shared species pairs alone. No new genus-impact audit is performed.

**Checks:** all40 final values and numerators/denominators match, maximum difference0. B/C/D also pass against separately copied/rebased input roots. The final four-panel renderer consumes52 saved records and displays40; its protected resources/fonts/layout still require separate packaging. It was not run here. Numeric reconstruction must not be described as new Figure4 artwork rendering.

## Supplementary Figure S1, Table 1 and Table S2

[The null replay](../analysis/paper_support/README.md#null-replay-s1-table-1-and-table-s2) reads **11 literature sets ×2 models ×1000 existing outcomes**, configs/seeds and observed summaries. It reuses the original `summarize()` and checks the source SHA before execution.

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_support/replay.py nulls   --literature-root data/paper-inputs/literature   --chen-wallen-root data/paper-inputs/chen_wallen   --output work/paper/nulls
```

**Outputs/checks:**22 recomputed summaries and excess-recovery display values, per-iteration seed/config checks, `RESULT.json`, `CHECKS.tsv` and hashes. Sample SD/linear quantiles and upper-tail `P=(1+n_ge_observed)/1001` are preserved. Compare S1 intervals as observed−null97.5% to observed−null2.5%. All22 summaries/S1 records agree (maximum serialization-level difference1.25e-14); Table1's11 rows and TableS2's22 rows agree at displayed precision. Eight sets are significant under both nulls; Guo Control/MECFS and exploratory MHO only under NullB.

This is **stored-simulation summary replay**. The [original simulation entries](../analysis/paper_support/README.md#rerunning-simulations) remain available, but full graph/source-exclusion inputs and resource approval are needed for new simulations. Two historical seed families have different master seeds and overlapping integer ranges; no new cross-family independence guarantee is asserted. Final Arial S1 artwork was not executed here.

## Supplementary Figure S2 and Table S4

The [support replay](../analysis/paper_support/README.md#support-replay-s2-and-table-s4) uses per-pair publication membership/source-excluded support,859 frozen tested-species sets, mapping decisions and the original653071-pair human-gut network.

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_support/replay.py support   --literature-root data/paper-inputs/literature   --rank-reference data/paper-inputs/rank_reference/human_gut.info   --output work/paper/support
```

**Outputs/checks:**6067 primary/6474 sensitivity pair rows; group summaries, rank/fraction statistics, original tests and input/source hashes. Primary report-count groups1–5 have5114,749,144,55,5 pairs. Supporting-project medians are10,35,117,147,176. Count rho=0.35635554668211467; fraction rho=0.26727784815927857. Pair/summary/test fields match the archived results; current TableS4A–D agrees at displayed precision.

Support counts exclude the pair's known source projects. Fractions divide by remaining projects whose actual FastSpar input contains both taxa. Only859/875 human-gut projects have complete tested-set evidence; the other16 are omitted, not treated as absence. Ranks use **all653071 original pairs before source exclusion**, descending counts, average ties and `100*(N-rank)/(N-1)`. The current MethodS7 and workbook now agree on this baseline. Final S2 artwork contains count/fraction only; rank results remain for TableS4. Its legacy Arial renderer was not rerun.

## Other tables and presentation

Figure1 is workflow/site presentation, not an unimplemented cohort statistic. TableS1 contains source taxonomy-mapping decisions; regenerating them requires the cited original third-party inputs and held mapping helpers, not a new inferred mapping. TableS3 is curated evidence for the17 human-gut candidates; code does not independently establish experimental validity. The current final supplementary files are recorded by hash, not silently bundled or substituted. No missing TableS7/S8 generator is asserted.

## Recorded validation

| Route | Actual scope in this task | Status |
|---|---|---|
| PRJEB81106 | Full dependency preflight; both-rank prepare+EPI; all1254 genus cells; valid zero candidates | PASS; fresh DKI not repeated; prior16-stage evidence retained |
| Figure2 |180 real association tables→90 summary rows;1260 numeric checks; final renderer regenerated120 plot records | PASS; summary/plot TSV byte identity |
| Figure3 | Old+replacement and final-only inputs;66 aggregate values;138 support records | PASS; no candidate inference/scoring |
| Figure4 | Fixed-order processed-support check; original B/C/D stats; relocation and40-value check | PASS; no full-D0 reconstruction/new networks/artwork |
| S1/Table1/S2 |22000 saved random outcomes→original summaries/P | PASS; no new random simulation |
| S2/TableS4 | Original653071 pairs+frozen support/tested sets→pair records/tests | PASS; no new859-matrix scan/artwork |
| Lightweight regression |9 scientific groups,12 upstream safety,7 portable matrix I/O plus new adapter safety checks | PASS; not full-cohort execution |

[VALIDATION_RECORDS.tsv](paper/VALIDATION_RECORDS.tsv) records input/output identities, scope and execution evidence. The complete detailed logs remain in the author's external publication-preparation directory; they are not GitHub payload. [Remaining work](FOLLOW_UP_ITEMS.md) concerns source access, reviewed data deposition, licenses and freezing the author-approved commit. None is replaced by a successful server-side command.
