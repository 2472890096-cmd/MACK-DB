# Figure 4: fixed-order recovery and Guo/Wallen comparisons

These entries recheck existing processed data. They never infer a network, run
FastSpar/bootstrap, train a model, or regenerate the historical 720 orders.
Figure 4 is still named **Figure7** in the source archive. Original scientific
producer files remain unchanged.

## What is available

| Entry | Starting input | Actual validation |
|---|---|---|
| `fixed_order_recovery.py` | Six frozen publication pair sets, the original-database publication-pair support lookup, nine-source exclusion audit | All six fixed-order steps recomputed; 6,067/11,900 |
| `run_original_statistics.py --mode b` | Five saved FastSpar matrix outputs, exact Run/taxon mappings, frozen literature workbooks and mapping decisions | Original tracked Stage5B functions executed; five species comparison summaries reproduced |
| `run_original_statistics.py --mode c` | Existing complete-project matrices, published mapped networks and common analysis taxa | Original SHA-pinned external producer executed; complete-project species statistics reproduced |
| `run_original_statistics.py --mode d` | Existing complete-side species edge sets and frozen species-label mapping | Original SHA-pinned external producer executed; genus-group statistics reproduced |
| `validate_statistics.py` | Those four freshly computed outputs | All 40 final display values and their relevant numerators/denominators match, maximum difference 0 |

The actual input archive and held C/D source files have **no public download
identifier yet**. This checkout therefore supplies runnable entry points and
small references, not a self-contained Figure 4 data distribution. The external
sources retain a per-file publication-rights hold. Do not substitute new scripts
with the same filename: source SHA256 is checked before execution.

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

The commands below are executable interfaces. They require the separately
provided real input data and, for C/D, the held original source plus its adjacent
`common.py`. Set `$PAPER_DATA` to the unpacked
`analysis_inputs/10_13_group_meeting_followup` directory and `$PAPER_SOURCES` to
the separately authorized original source directory before running them.

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
  --source "$PAPER_SOURCES/complete_project/compare_networks.py" \
  --output work/paper/Figure4C
"$MACKDB_PAPER_PYTHON" "$module/run_original_statistics.py" --mode d --config "$cfg" \
  --source "$PAPER_SOURCES/genus_groups/genus_overlap.py" \
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
maximum absolute difference **0**. Six adapter safety tests also passed:

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
unpublished input/source archive is already available to another user. The
final renderer, `plot_figure7_four_panels.py`, reads existing 52-row plotting
input/metadata and emits 40 displayed entries. It still needs its original
layout, protected reference files and Arial fonts, and has no portable CLI.
It was **not executed** here. The small
[40-value reference](reference/Figure4_expected_values.tsv) contains aggregates
only, without private source paths. It is a checking reference, not a
replacement for statistical input data.

See the [paper reproduction guide](../../docs/PAPER_REPRODUCTION.md) for current
download/permission gaps and the distinction between fresh statistics and TSV
redrawing.
