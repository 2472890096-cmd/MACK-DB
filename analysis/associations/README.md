# Environment association statistics (Figure 2)

This entry starts with the paper's processed, environment-level association
tables. It does not infer networks, download reads, or estimate project-level
correlations. The fixed input has 45 environments at each of two taxonomic
levels for each of FastSpar and fastCCLasso: 180 `*.info` files in four
directories. Keep the project-support columns (`proj_ID`, `proj_corr`,
`proj_pv`, `proj_sam`) alongside the means (`aver_corr`, `aver_pv`) and
`proj_count`; these are needed to inspect method-specific support and rounding.

## Input and source availability

The two original scientific sources are now included **byte-for-byte** under
[`original/`](original/), following the author's 23 November 2026 authorization:

| Source | SHA256 |
|---|---|
| [Statistics producer](original/plot_species_genus_coabundance_complete.py) | `48920f30e40f7212fb73c3394b6718a5634acd6d5d5f4554fa5d6f97f29dfb9c` |
| [Final renderer and environment selection](original/plot_coabundance_panel_ABCDE_ARIAL.py) | `50dda8e6c2169aa7fd864a8d68f230ac77535529d4b79201a09ff1929956c14a` |

Use the adapters below, which default to these repository-local files and retain
all SHA256 checks. The originals preserve historical path defaults for provenance;
do not run them directly against those paths. The large association inputs remain
a separate proposed supporting-data archive with no public download identifier yet.
See [paper reproduction](../../docs/PAPER_REPRODUCTION.md) for data access status.

The local review archive is
`Figure2_processed_environment_associations_20261112.tar.gz` (507,859,601 bytes;
SHA256 `4093c398bc8ae79cdfdab3e64ebfb19c41a65aff0b093b5ef5ccb3b312fcd0f2`).
It is prepared for author review and has **not** been uploaded or licensed for
redistribution. The [input metadata manifest](INPUT_DATA_MANIFEST.tsv) records
182 research files (180 association inputs and two small references), with
relative paths, sizes and SHA256. Once legitimately obtained and extracted,
check it before running:

```bash
"$MACKDB_PAPER_PYTHON" tools/verify_paper_inputs.py \
  --manifest analysis/associations/INPUT_DATA_MANIFEST.tsv \
  --root "$FIGURE2_DATA_ROOT"
```

The four input directories within that root are
`figure2/{species,genus}/{FastSpar,CCLasso}/f2/`; references are under
`figure2/reference/`.

Use Python3.11+ via `MACKDB_PAPER_PYTHON`. Required Python packages are NumPy, pandas, SciPy and Matplotlib. The final
renderer import used for selection also requires Pillow. The numerical
validation does not load Arial font files or call the rendering entry point.

## Recompute from actual association tables

After obtaining the listed data, set the four input-directory paths for your
machine. Source paths need no configuration: the adapters locate `original/`
relative to their own files. Explicit `--statistics-source`, `--selection-source`
and `--renderer-source` overrides remain available, and must match the same
frozen source SHA256.

Run from the repository root, using a new output directory:

```bash
"$MACKDB_PAPER_PYTHON" analysis/associations/reproduce_environment_statistics.py \
  --species-fastspar "$FIGURE2_SPECIES_FASTSPAR" \
  --genus-fastspar "$FIGURE2_GENUS_FASTSPAR" \
  --species-cclasso "$FIGURE2_SPECIES_CCLASSO" \
  --genus-cclasso "$FIGURE2_GENUS_CCLASSO" \
  --output work/paper/figure2-statistics
```

Outputs are `environment_method_summary.tsv` (90 rows),
`support_count_distribution.tsv`, `global_rank_summary.tsv`,
`INPUT_SHA256.tsv` and `EXECUTION.json`. The implementation directly calls the
archived producer's `build_summaries()`, `global_summary()` and `disp_env()`.

In Figure 2E each method's already-saved **signed** `aver_corr` is joined on
the common unordered pair. The archived `safe_spearman()` drops missing means,
retains zero means and returns NA for fewer than three usable pairs or a
constant vector. The sign-agreement branch separately excludes zero means;
its mask must not be substituted for the Spearman mask. Different methods
can have different supporting BioProjects for the same association.

## Check the frozen statistics and final panel selection

`FIGURE2_REFERENCE_SUMMARY` denotes the original 90-row summary;
`FIGURE2_PLOT_VALUES` denotes the final 120-row `panel_plot_values.tsv`.

```bash
"$MACKDB_PAPER_PYTHON" analysis/associations/verify_environment_statistics.py \
  --summary work/paper/figure2-statistics/environment_method_summary.tsv \
  --reference-summary "$FIGURE2_REFERENCE_SUMMARY" \
  --plot-values "$FIGURE2_PLOT_VALUES" \
  --output work/paper/figure2-validation
```

The checker validates all 90 environment/rank records, the actual top-15
selection and all 180 numerical marks (three segments per A/B bar and the
C/D/E points), including the 30 E points. Output files are
`VALUE_COMPARISONS.tsv` and `VALIDATION.json`; inspect the status and mismatch
counts. This is statistical reproduction and numeric panel checking, not
an assertion that the final Arial artwork was regenerated.

From the fixed inputs, global environment-level totals are:

| Level | FastSpar | fastCCLasso | Common | Union |
|---|---:|---:|---:|---:|
| Species | 4,661,190 | 2,287,364 | 2,059,581 | 4,888,973 |
| Genus | 999,851 | 646,124 | 519,093 | 1,126,882 |

Means remain at their historical stored precision. The saved FastSpar
**project coefficients** are typically four decimals; their means can have
more digits (for example `0.55145`). fastCCLasso coefficients and saved
means retain about ten significant digits. Do not round the summary's
Spearman values before comparison.

## Redraw from a saved summary

This separate entry only redraws the fixed summary; it does not recompute
association statistics. The original final renderer is included; its four Arial TTF
files must be obtained separately with appropriate font permissions.
`FIGURE2_FONT_DIR` must contain `Arial.TTF`, `Arialbd.TTF`, `Ariali.TTF`
and `Arialbi.TTF`.

```bash
"$MACKDB_PAPER_PYTHON" analysis/associations/render_environment_figure.py \
  --summary work/paper/figure2-statistics/environment_method_summary.tsv \
  --font-dir "$FIGURE2_FONT_DIR" \
  --output work/paper/figure2-render
```

The adapter configures only the output directory, font directory and the
original `--summary` argument. It retains the original layout/selection and
numeric-coordinate assertions. Outputs are the four final image formats,
`panel_plot_values.tsv`, `render_audit.json` and `RENDER_EXECUTION.json`.
Compare the generated panel values independently from successful rendering.

## Manuscript presentation edits

The author confirmed that the manuscript's final axis wording and small legend
changes were applied manually and did not change the original calculations.
The included original sources therefore retain their original labels, including
Figure 2E's `mean r` rather than the manuscript's `average r`. Statistical values
and all five data-plot regions matched the current manuscript; exact final
annotation/legend pixels are not promised by the original renderer. These
presentation differences do not change source identity, formulas or thresholds.
