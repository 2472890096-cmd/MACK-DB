# Figure 4: fixed-order recovery and Guo/Wallen comparisons

These entries recheck existing processed data. They never infer a network, run
FastSpar/bootstrap, train a model, or regenerate the historical 720 orders.
Figure 4 is still named **Figure7** in the input archive. Scientific functions
and thresholds are unchanged; source-layout adaptations are recorded separately.

## What is available

| Entry | Starting input | Actual validation |
|---|---|---|
| `fixed_order_recovery.py` | Six frozen publication pair sets, the original-database publication-pair support lookup, nine-source exclusion audit | All six fixed-order steps recomputed; 6,067/11,900 |
| `run_original_statistics.py --mode b` | Five saved FastSpar matrix outputs, exact Run/taxon mappings, frozen literature workbooks and mapping decisions | Original tracked Stage5B functions executed; five species comparison summaries reproduced |
| `run_original_statistics.py --mode c` | Existing complete-project matrices, published mapped networks and common analysis taxa | Original SHA-pinned producer executed; complete-project species statistics reproduced |
| `run_original_statistics.py --mode d` | Existing complete-side species edge sets and frozen species-label mapping | Original SHA-pinned producer executed; genus-group statistics reproduced |
| `validate_statistics.py` | Those four freshly computed outputs | All 40 final display values and their relevant numerators/denominators match, maximum difference 0 |

The Figure 4B producers and helpers now live under `original/matched_subgroups/`,
alongside the C/D sources under `original/complete_project/` and
`original/genus_groups/`. The B files were relocated from a historical `logs/`
directory: they are scientific code, not log files. Their data-root initialization
now uses the existing `stage5b` path configuration rather than source-file parents;
their scientific functions and constants are unchanged. The C/D producers and
their separate `common.py` helpers remain byte-identical.
The author approved their integration on 23 November 2026 and confirmed that
manual axis-label changes in the manuscript did not alter their calculations.
The adapter selects these bundled sources by default and checks the approved
source SHA256 before execution. `--source` remains an optional override accepting
only the same pinned bytes for the current release.

The actual input archive still has **no public download identifier**. Code is
included, but the checkout is not a self-contained Figure 4 data distribution.
Original source files retain historical path defaults/import-time IO; use the
adapter below rather than running `original/*.py` files directly.

### Source layout and historical Stage names

| Location | Role |
|---|---|
| `original/matched_subgroups/stage5b2_analysis.py` | Figure 4B: five matched-Run groups, common taxa and threshold-specific overlap |
| `original/matched_subgroups/build_reference.py` | Build the frozen R226 vocabulary used by those comparisons |
| `original/matched_subgroups/run_degree_null.py` | Optional formal 1,000-replicate degree-preserving null generation; not run by the statistics adapter |
| `original/matched_subgroups/independent_recalculation.py` | Independent readback of subgroup statistics and existing null graph archives |
| `original/complete_project/` | Figure 4C complete-project species comparison |
| `original/genus_groups/` | Figure 4D and B genus-group comparison |

The four matched-subgroup sources use `resolve_path('stage5b',
'04_network_comparison')` for their historical data workspace. Configure `stage5b`
through `MACKDB_PATHS_CONFIG` only when deliberately using those historical
preparation/validation entries. The recommended `run_original_statistics.py`
adapter instead supplies isolated input/output roots through its explicit JSON.
Neither route infers a scientific data directory from the moved source location.

The separate `external_validation/literature_networks/stage4b_recurrence_recovery/`
module performs publication-pair source exclusion and recovery. Its **Stage 4B**
is a historical processing label, not the current **Figure 4B** panel; it remains
part of the literature-support chain for Supplementary Figure S2/Table S4.

## Inputs and paths

From the repository root, copy the [configuration template](../../config/paper_network_comparison.example.json)
and point `base` and `project_inputs` to the unpacked input roots. Paths in JSON
are relative to that JSON file, independent of the current working directory.

```bash
cp -n config/paper_network_comparison.example.json config/paper_network_comparison.local.json
```

The required layout preserves these original directory names under `base`:

```
9_16_external_validation/               # publication inputs/mapping/provenance
9_22_stage5B_same_data_cross_method/    # matched Run inputs and saved outputs
10_13_group_meeting_followup/          # complete-project and genus pair inputs
figure56/                             # frozen Wallen mapping/component evidence
```

`project_inputs` contains the two projects' original species abundance/count/
taxon mapping files, saved FastSpar matrices and required production-command
evidence. The exact executed-input manifest, including file hashes, is recorded
in the accompanying publication-preparation audit; publication and download
packaging are pending. This is not a request for all cohort matrices.

The optional `path_rebase` dictionary maps each **original frozen metadata
prefix** to its local unpacked root. It is applied in memory to path fields,
historical manifest keys and `comparison_datasets.json`; expected SHA256 values
and scientific data fields are unchanged. It also avoids reusing old server
paths accidentally. Supply the actual archived prefixes from the eventual
data-package metadata, replacing the illustrative prefixes in the template.
An optional `comparison_datasets` key selects an explicit JSON descriptor; its
`*_path` values may be absolute or relative to that descriptor. Never edit the
original frozen metadata in place.

The adapter compiles the pinned original function ASTs unchanged, supplies
explicit path/state bindings, and excludes top-level historical IO. It requires
Python 3.11+ (`hashlib.file_digest`), NumPy, SciPy and openpyxl. Use a new output
directory for each command; an existing destination is rejected.

## Run statistics

The commands below use the bundled accepted producers and their adjacent
helpers. They still require the separately provided real input data. Set
`PAPER_DATA` to the unpacked `analysis_inputs/10_13_group_meeting_followup`
directory before running them.

```bash
module=analysis/paper_network_comparison
cfg=config/paper_network_comparison.local.json
"$MACKDB_PAPER_PYTHON" "$module/fixed_order_recovery.py" \
  --pairsets-dir "$PAPER_DATA/01_stage6A_sequential_recovery/data" \
  --support-lookup "$PAPER_DATA/01_stage6A_sequential_recovery/data/stage6A_published_pair_support_lookup.tsv" \
  --source-exclusion "$PAPER_DATA/01_stage6A_sequential_recovery/stage6A_source_exclusion_audit.tsv" \
  --archive "$PAPER_DATA/01_stage6A_sequential_recovery/stage6A_all_720_orders.tsv" \
  --output work/paper/Figure4A

"$MACKDB_PAPER_PYTHON" "$module/run_original_statistics.py" --mode b --config "$cfg" \
  --output work/paper/Figure4B
"$MACKDB_PAPER_PYTHON" "$module/run_original_statistics.py" --mode c --config "$cfg" \
  --output work/paper/Figure4C
"$MACKDB_PAPER_PYTHON" "$module/run_original_statistics.py" --mode d --config "$cfg" \
  --output work/paper/Figure4D

"$MACKDB_PAPER_PYTHON" "$module/validate_statistics.py" \
  --a work/paper/Figure4A --b work/paper/Figure4B \
  --c work/paper/Figure4C --d work/paper/Figure4D \
  --output work/paper/Figure4_validation
```

Run A independently of B/C/D. B/C/D each start from existing processed inputs;
the adapter does not rerun prior network inference or depend on new FastSpar.
The legacy B producer prints `NULL_PENDING`; this means the edge-comparison
stage ended and **no new null model ran**. `ADAPTER_EXECUTION.json` records that
boundary. In the historical analysis, B/C/D input edge files were frozen at
different production steps; their source identities are retained.

Figure 4A's transparent set checker is an independent, narrowly scoped check,
not a rerun of the original 720-order production script. It recalculates the
nine-project support exclusion from the saved per-pair original support IDs,
requires at least two remaining supporting projects, and applies the fixed order
**Chen → Wallen → Li → Guo → Lacruz-Pleguezuelos → Zhou**. `Lacruz` is the stored
study ID. The lookup covers literature pairs present in the original audited
f≥2 database. It does not re-scan the full database or independently establish
the archived full D0 size **652,972**. Absence from this lookup is not asserted to
mean zero underlying support. The six cumulative counts are
**2,045, 3,284, 3,550, 3,712, 4,805, 6,067**, with cumulative publication-union
sizes **2,733, 4,396, 4,871, 5,142, 6,392, 11,900**.

## B: sample selection and thresholds

[Figure4B_sequencing_runs.tsv](reference/Figure4B_sequencing_runs.tsv) contains
the requested public Run identifiers and retained ordering: 91 Guo Control,
106 Guo MECFS, their 197-Run Common union, 490 Wallen PD and 234 Wallen NHC.
There are **1,118 network–Run rows and 921 distinct Runs**. The Common rows
intentionally repeat the disjoint Guo subgroup Runs. Participant identifiers,
BioSample/clinical metadata and private paths are excluded from this list.

- MACK-DB Guo: `abs(r) >= 0.3 AND P <= 0.0001`.
- MACK-DB Wallen primary: `abs(r) > 0.2 AND P < 0.05`.
- Published Guo: supplied final single-fit FlashWeave-S edge list, adjusted
  P<0.05 as documented in the source Methods; no new weight cutoff.
- Published Wallen: source tables already use uncorrected P<0.05; retain
  `abs(r)>0.2`. Raw signed r/P values and frozen mapping decisions are reused.
- Guo Common check: species **31/130**; genus groups **28/101**.

All B overlaps use the same mapped common-analysis species universe before
comparison. Complete-project C uses 197 Guo and 725 Wallen Runs; the five B
groups do not imply a pooled 724-Run Wallen inference.

## C/D: correct universes and genus grouping

C's three conditions, in order, are Guo aggregate published-network union,
Wallen published-matched threshold, and Wallen MACK-DB-standard threshold.
The MACK-DB-standard threshold is `abs(r)>=0.3 AND P<=0.0001`; published-matched
is `abs(r)>0.2 AND P<0.05`. The published Wallen set is unchanged between those
two conditions. The plotted counts are:

| Condition | Published | Complete project, size only | Comparable project | Shared |
|---|---:|---:|---:|---:|
| Guo species | 201 | 1003 | 73 | 38 |
| Wallen species, published-matched | 1413 | 7525 | 733 | 523 |
| Wallen species, MACK-DB standard | 1413 | 1246 | 201 | 185 |
| Guo genus groups | 149 | 615 | 57 | 34 |
| Wallen genus groups, published-matched | 796 | 3722 | 440 | 307 |
| Wallen genus groups, MACK-DB standard | 796 | 699 | 128 | 118 |

The orange complete-project bar is an input-scope network size. Overlap uses
the published **common-analysis** denominator (201/1413 species or 149/796
genus groups), not that orange size. Jaccard uses the union within the common
analysis universe. For D, orange uses `ALL_MAPPED_SPECIES`; other bars use
`COMMON_SPECIES`.

Genus groups are obtained by applying the frozen production-label token rule
separately to **both complete side-sets** of existing species associations,
removing intra-group loops and deduplicating unordered genus-group pairs, then
intersecting. This folds some GTDB suffixes. It is neither a newly inferred genus
FastSpar network nor conversion of shared species edges only; no coefficient
averaging or additional sign-concordance condition is applied.

## Validation and drawing boundary

On 2026-11-12, B/C/D were executed in isolated destinations using original
scientific functions, then repeated against separately copied, relocated input
roots with path rebasing. Both validations matched **40/40** final values with
maximum absolute difference **0**. The original adapter safety tests passed. Current tests also check all three bundled producer paths, both C/D helper hashes, and the configured data roots of all four relocated B sources:

```bash
"$MACKDB_PAPER_PYTHON" analysis/paper_network_comparison/test_io_adapter.py
```

The first relocation attempt exposed an IO prefix-remapping bug when the new
root was nested under an old prefix; the adapter now makes translation
idempotent and excludes its output root. A second preparation gap was a missing
historical comparison manifest; it was added to the isolated input set. Those
failed preparation logs are preserved externally. No scientific function,
expected result or source input was changed to obtain PASS.

This validates statistics from existing data. It does not assert that the
unpublished input archive is already available to another user. The
final renderer, `plot_figure7_four_panels.py`, reads existing 52-row plotting
input/metadata and emits 40 displayed entries. It still needs its original
layout, protected reference files and Arial fonts, and has no portable CLI.
The original four-panel renderer was executed separately on 16 November 2026
and reproduced the archived PNG byte-for-byte; it is not part of these six
required source integrations. Its B-axis wording differs from the paper's
manually edited title, which the author has confirmed is a presentation-only
change. Neither calculation nor saved reference values were changed. The small
[40-value reference](reference/Figure4_expected_values.tsv) contains aggregates
only, without private source paths. It is a checking reference, not a
replacement for statistical input data.

See the [paper reproduction guide](../../docs/PAPER_REPRODUCTION.md) for current
download/permission gaps and the distinction between fresh statistics and TSV
redrawing.

## Source integration check

On 23 November 2026 the four original C/D files were copied without byte changes.
The adapter's default-source paths were then used for isolated C/D statistical
replays against the already prepared relocated input roots. Original SHA guards,
AST functions, thresholds and input hashes remain unchanged. This is a bounded
path-integration check, not another network inference or model run.

The 23 November default-source runs reproduced the complete-project 5-row summary
and genus-group 20-row summary exactly against the relocated reference outputs.
The combined 40-value check passed again with maximum absolute difference 0;
no `--source` argument or external source path was used.
