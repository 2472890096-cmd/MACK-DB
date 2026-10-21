# Association inference

Two independent abundance-based routes are retained. **PASS** is the Stage2 single-project species route; genus is **ARCHIVED_REPLAY**. Neither is a whole-cohort/clean-install claim.

FastSpar1.0.0 requires external `fastspar`, `fastspar_bootstrap`, `fastspar_pvalues` and GNUparallel. Use the count matrix and `fastspar/run_fastspar.sh COUNT OUTPUT_PREFIX`, setting executable paths via FASTSPAR_BIN/FASTSPAR_BOOTSTRAP_BIN/FASTSPAR_PVALUES_BIN/PARALLEL_BIN and resource-only MACKDB_THREADS. Main iterations50/exclusion20/internal0.2;10000bootstraps with five fit iterations;10000permutations. Convert labelled square TSV matrices with `tools/network_matrix_io.py`, then apply `sparcc_corr_pvalue.R` and `corr_inter.R`. Final |r|≥0.3 and P≤1e-4. Outputs are corr/cov/P, named selected pairs and per-project links. Zi-Pi consumes the network; environment integration consumes links plus original coefficients/P.

fastCCLasso wrappers require a separately obtained official Rcore (LGPL; not bundled), abundance/map directories, nboot100,kmax20,isCntFALSE,pseudocount1e-9, rank-specific seeds20260804/20260806. Initial0.1 filtering is followed by the accepted0.3 refilter. The official P derives from the t-tail of averaged original+bootstrap estimates as min(u,1-u), without multiplication by2; it is not conventional doubled two-sided P or sign-frequency. Archive percentage inputs are transposed and pseudocount-added as the unchanged runner does. Final0edges with passstatus is legitimate.

Stage2 full FastSpar matrices matched the independent species reference and production genus reference exactly; full fastCC matrices matched production with zero retained edges. Tool licenses remain external and are not reassigned.

Read [PIPELINE](../../PIPELINE.md) for exact commands and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for every helper/CLI. [Input/output contracts](../../docs/PIPELINE_INPUT_OUTPUT_CONTRACTS.tsv) and [status evidence](../../docs/REPRODUCIBILITY_STATUS_MATRIX.tsv) define the acceptance scope.
