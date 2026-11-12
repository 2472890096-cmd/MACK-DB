# Replay the stored null and literature-support analyses

This module reproduces **Supplementary Figure S1, Table 1, Table S2** and
**Supplementary Figure S2 / Table S4** from existing processed research inputs.
It neither reruns the 22,000 randomized networks nor infers associations. The
unchanged scientific functions come from the retained
[Chen/Wallen protocol](../external_validation/shared_protocol/chen_wallen_null_models_final_v1.py)
and [support analysis](../external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py).
Only named function definitions are loaded; their historical initialization and
`main()` functions do not run. These two existing scientific sources remain
byte-identical to the installation-fix baseline.

## Dependencies and inputs

Use Python with NumPy, pandas, SciPy and statsmodels. The recorded replay used
NumPy 2.1.3, pandas 2.2.3, SciPy 1.15.3 and statsmodels 0.14.4; it does not claim
that the earlier production used this same complete environment. Set
`MACKDB_PAPER_PYTHON` to an interpreter containing these dependencies. This is
separate from the optional DKI model interpreter.

[INPUT_ARCHIVE_MANIFEST.tsv](INPUT_ARCHIVE_MANIFEST.tsv) records the required
processed-data layout, bytes and SHA256. A local candidate archive has been
assembled, but **no public download URL or data DOI exists yet**. Distribution
requires author/source-rights review. Do not substitute participant metadata,
raw reads or an unrelated current database export. Create `data/paper-inputs/` and extract the archive there (it does not contain a `paper-inputs/` prefix). The resulting layout is:

```text
paper-inputs/
  literature/                 # frozen Stage3C, Stage4B/4D and S1 values
  chen_wallen/                # four original 1000-replicate tables + summary
  rank_reference/human_gut.info # original 653071-pair FastSpar reference
```

Each output directory must be new. From the repository root:

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_support/replay.py nulls \
  --literature-root data/paper-inputs/literature \
  --chen-wallen-root data/paper-inputs/chen_wallen \
  --output work/paper-null-replay

"$MACKDB_PAPER_PYTHON" analysis/paper_support/replay.py support \
  --literature-root data/paper-inputs/literature \
  --rank-reference data/paper-inputs/rank_reference/human_gut.info \
  --output work/paper-support-replay
```

These commands are actual implemented entry points; without the external
archive they stop at the first missing input. Existing output paths are refused.
`INPUT_SHA256.tsv`, `CHECKS.tsv`, reused function hashes and `RESULT.json` document
each execution. Nonzero exit and `FAIL` must not be relabeled as success.

## Null replay: S1, Table 1 and Table S2

The input contains **11 association sets × 2 null models × 1000 iterations**.
The original `summarize()` calculates sample SD, NumPy linear quantiles,
enrichment, the original zero-SD `z=inf` rule and upper-tail empirical
`P=(1+n_ge_observed)/1001`. It compares every rebuilt summary to the saved
scientific summaries and the 22 S1 display records. S1 intervals use
`observed_rate - null_q975` and `observed_rate - null_q025`; they are not
observational confidence intervals. The aggregate reference is
[NULL_SUMMARY_22.tsv](reference/NULL_SUMMARY_22.tsv).

The Chen/Wallen master seed is 20260910; the nine-network extension uses
20260917. Within each historical family, the recorded offsets and all 1000
iteration seeds are checked. These two jobs have overlapping integer seed
ranges: `SEED_SCOPE.json` records 22000 family/seed keys and 18028 globally
distinct integer seeds. No new independence guarantee is asserted and no saved
randomization is changed.

### Rerunning simulations

Saved iterations are sufficient to recompute the table/S1 summary. **Rerunning
the simulations** needs additional source-excluded edge sets, mapping inputs,
prepared graphs and historical protocol evidence used by the
[original Stage3C entry](../external_validation/literature_networks/stage3c_null_model/scripts/03c_null_models.py).
That expensive entry is deliberately not invoked by this lightweight replay;
its full rerun dependency package is not bundled here.

## Support replay: S2 and Table S4

The support branch starts from the frozen Stage4B pair-level recovery/support
tables, **859 actual tested-species sets**, the exact database-label decoder and
the entire original human-gut network. It checks the network SHA256, scans all
**653071** unique pairs and reconstructs ranks using the original function:
descending **pre-exclusion** support counts, average tie rank, percentile
`100*(N-rank)/(N-1)`.

For each recovered pair the numerator is frozen remaining support; the
denominator is the number of projects testing both taxa minus that pair's known
source projects. The 16 unavailable projects among 875 listed human-gut projects
are excluded, not counted as absence. This uses archived tested-set evidence;
it does not newly inspect or reconstruct all 859 original matrices.

The branch rebuilds **6067 primary** and **6474 sensitivity** pair rows,
support/percentile/fraction summaries and original statistical tests, and checks
them against archived outputs. Reference aggregates:
[primary summary](reference/SUPPORT_PRIMARY_SUMMARY.tsv) and
[primary Spearman](reference/SUPPORT_PRIMARY_SPEARMAN.tsv).
Key primary checks are counts **5114, 749, 144, 55, 5** for one to five reports;
support medians **10, 35, 117, 147, 176**; count rho
**0.35635554668211467**, fraction rho **0.26727784815927857**.

The final S2 artwork shows only support count and support fraction; the rank
analysis remains available for Table S4. Graphical redrawing from the saved
values is a separate step. Historical Arial renderers require licensed fonts
and legacy directory resources and are not represented as portable renderers.

## Actual verification

On 2026-11-12, final real-data null replay took **1.52 s / 201268 KiB maximum RSS**;
support replay took **17.92 s / 453804 KiB**. All 22 S1 records matched within
**1.25e-14 absolute numerical difference** (serialization/reduction precision).
All regenerated support pair/summary/statistics fields matched the frozen
records exactly. Current Table 1's 11 rows, Table S2's 22 rows and all sections
A–D of Table S4 also matched their displayed precision. No fresh randomization or
plot rendering was performed. These results establish server-local statistical
replay, not public input availability or a full-cohort independent rerun.

Small adapter tests:

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_support/test_replay.py
```

Five tests exercise the original sample-SD/upper-tail/zero-SD behavior, changed
scientific-source refusal, duplicate-key refusal, output escape refusal and
refusal to overwrite an existing result.
