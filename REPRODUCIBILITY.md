# Reproducibility scope

There are two independently documented routes:

1. [PRJEB81106](examples/real_bioproject/PRJEB81106/README.md) starts from the bundled **33 archived Runs**, species53×33 and genus38×33 matrices. It can run FastSpar→Zi-Pi, fastCCLasso, EPI v3 and fresh CPU DKI training/prediction/unified-R scoring. EPI and DKI are independent matrix branches. See [installation](environments/portable-matrix/README.md) for the pinned SciPy/PyPI, R version-comparison and fastCCLasso instructions.
2. [Paper reproduction](docs/PAPER_REPRODUCTION.md) starts from the **processed inputs actually used for the final figures/tables**. It recomputes association summaries, candidate-set concordance, literature/network comparisons and summaries of saved random outcomes. This route does not repeat metagenomic profiling, whole-cohort network inference or DKI training.

The separate [ERR14036375](examples/real_fastq/ERR14036375/README.md) paired-FASTQ route is an accepted one-Run upstream reproduction through MetaKSSD2.21, L3K11, GTDB R226 Query (`-p 1 -m 2`) and Species Profile (`sharekmer=12`). Its final Profile matched the archived file byte-for-byte. It does not establish fresh upstream execution for all33 project Runs.

## Scientific versions and scope

The paper describes **264,767 shotgun metagenomic sequencing runs**,1,575 BioProjects and45 environments. The main candidate comparison uses the fixed1,571-project scope;1,570 and1,568 are sensitivity scopes. Eligibility is not reconstructed from candidate counts. Valid zero-candidate projects remain included.

Scientific definitions remain EPI v3 and DKI **MACKDB_DKI_R_20261012_v1**. Genus abundance uses the author-adopted [species-label token rule](docs/GENUS_AGGREGATION_STATUS.md). Figure4 genus groups coarsen existing species association sets; they are not new genus FastSpar networks.

Current numbering is main Figures1–4/Table1, MethodsS1–S8, FiguresS1/S2 and TablesS1–S4. Historical Stage numbers are development stages, not manuscript figure labels. Legacy paths named Figure5/6/7 retain historical names. MethodS7 concerns publication recurrence/support; S8 concerns matched-subgroup and complete-project comparisons. Figure4A uses a fixed six-study order, not a fresh720-order experiment.

## What execution evidence establishes

The historical PRJEB81106 full run completed16 stages across both ranks with fresh model training. An independent new-environment run at `9fbb1c6647b2f3b6fe7e2a60b9528a6c51446a82` completed after disclosed R-version/C++ library handling. Commit `9f4fce74c0d04716903c25caff612ffd1661c699` corrected installation-source/version checks without changing scientific code. Neither statement promises identical floating-point model predictions on every system or a clean installation pass on every Linux distribution.

The paper checks on this branch use isolated output directories and real frozen inputs. Their outputs are compared to independently archived summaries and final plot-value tables. [The per-result record](docs/PAPER_REPRODUCTION.md#recorded-validation) distinguishes actual statistics, stored-result replay, artwork and prior evidence. Nine lightweight scientific groups test specified formulas/boundaries, not all1,575 projects or every historical production environment.

## Sources, data access and exact version

The [data-package manifest](docs/paper/DATA_PACKAGE_MANIFEST.tsv) records proposed archive files, hashes and access status. The required Figure2/4 statistics/selection sources are included with source identities and repository-relative defaults. The Figure4B subgroup sources now share the network-comparison directory; their scientific functions retain the accepted definitions. The processed input packages still have no public download link, so a reader must obtain those separately. Other historical taxonomy helpers and optional final-artwork sources remain outside this migration. No unavailable DOI, source license or public-data permission is implied.

On23November2026 the author confirmed that the final Figure2/4 axis-label and legend edits were manual presentation changes, with no change to the calculations. The original source labels are preserved; pixel-identical final typography is not the statistical reproduction criterion.

Record `git rev-parse HEAD` with the input SHA256 manifests when using this version. The author must approve and archive the final submission commit; branch names are mutable. No formal Release or manuscript modification is performed by these checks.
