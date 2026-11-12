# MACK-DB

MACK-DB (Microbial Association and Candidate Keystone Database) integrates microbial associations and candidate keystone taxa from **264,767 shotgun metagenomic sequencing runs, 1,575 BioProjects and 45 environments**. Records retain their supporting BioProjects and inference methods. The main three-method candidate comparison uses a fixed **1,571-project** set, including successfully analyzed projects with no candidates.

The workflow includes MetaKSSD profiling, species/genus abundance matrices, FastSpar and fastCCLasso associations, Zi-Pi, EPI v3, unified-R DKI scoring, environment integration and external literature comparisons. **Zi-Pi uses a network; EPI and DKI start independently from abundance matrices.**

Choose one of the two reproduction routes below. Running a complete example and reproducing the paper-wide statistics are separate tasks, with separate input requirements and validation records.

## 1. Run a complete real BioProject: PRJEB81106

### Install

Follow the [Linux environment instructions](environments/portable-matrix/README.md). They preserve **SciPy 1.15.3 from official PyPI**, R package-version-object comparisons, and the fixed fastCCLasso commit/core checksum. Configure `MACKDB_PYTHON` as the Python≥3.11 orchestrator and `MACKDB_DKI_PYTHON` as the separate model interpreter. Select existing R/FastSpar/GNU parallel paths and obtain the pinned external fastCCLasso core as documented. The workflow does not install software or download reads automatically.

### Inputs

The repository [bundles](examples/real_bioproject/PRJEB81106/bundled/README.md) the **33 archived Run samples**, a **53×33 species matrix**, a **38×33 genus matrix**, exact taxon/Run mappings, input checksums, minimal taxonomy and small reference outputs. This route **starts from abundance matrices**. It needs neither FASTQ downloads nor saved DKI predictions. The archived 33 Runs define this example; current BioProject membership must not replace them.

### Run

Run from the repository root. Copy the template once, then edit/check local tool paths; preserve an existing local configuration.

```bash
cp -n config/portable_matrix.example.json config/portable_matrix.local.json
```

The template intentionally leaves interpreter selection to the exported environment variables unless you add explicit interpreter paths. After setup:

```bash
example=examples/real_bioproject/PRJEB81106/run_workflow.py
config=config/portable_matrix.local.json

"$MACKDB_PYTHON" "$example" --config "$config" --steps all --preflight

"$MACKDB_PYTHON" "$example" --config "$config" --steps all \
  --levels species genus --output work/PRJEB81106_complete --threads 4
```

Use a **new output directory**. The default limit is 600 seconds per stage; explicitly add `--timeout-seconds 3600` on slower servers if needed. The recorded independent run used 3600 seconds. An execution limit does not change scientific parameters. `--steps all` includes **new CPU cNODE2 training, prediction and unified-R scoring**; `--steps epi` selects an independent matrix-only branch.

### Check

Read `work/PRJEB81106_complete/WORKFLOW_STATUS.json` and each `species/STAGE/SUCCESS.json` / `genus/STAGE/SUCCESS.json`. Completion requires the recorded output checks and hashes, not merely an existing filename.

- Existing evidence covers **16 completed stages**, eight per rank.
- Archived FastSpar references contain **25 species** and **3 genus** selected associations. Compare numerical outputs separately from successful execution.
- EPI has legitimate zero candidates at both ranks. A valid zero-node genus Zi-Pi network is also a successful result.
- DKI must record new training, model/prediction output and unified-R scoring. The recorded example produced three candidates per rank; fresh prediction/checkpoint bytes are not promised identical across numerical environments. Saved-prediction replay is a separate historical entry.

The [complete example guide](examples/real_bioproject/PRJEB81106/README.md) documents parameters, outputs, expected results, failure handling and resume checks. The separate [ERR14036375 paired-FASTQ example](examples/real_fastq/ERR14036375/README.md) provides the accepted one-Run MetaKSSD→GTDB R226 Query→Species Profile upstream route; it does not replace the project matrix or establish fresh upstream execution for all 33 Runs.

## 2. Reproduce the paper

Use [Reproduce the paper](docs/PAPER_REPRODUCTION.md) for **Figures 2–4, Supplementary Figures S1/S2, Table 1 and related supplementary tables**. It lists processed inputs, source versions, commands, outputs, check values and download/permission status. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) explains how these checks differ from the real example and from fresh cohort-wide inference.

Real-input statistical checks on this branch reproduce the Figure 2 environment summary and displayed numbers, Figure 3 main/sensitivity concordance, Figure 4 comparisons, and supplementary null/support statistics. The example preflight and both-rank matrix/EPI checks also pass. Numerical checks and final artwork rendering have separate records. **Some required sources and real input packages remain held or unpublished**; these server-side passes do not make the checkout alone a complete public paper-reproduction package. See [input packages and access](docs/paper/DATA_PACKAGES.md).

## Validation and exact revision

| Route | Recorded execution / revision | Meaning for this checkout |
|---|---|---|
| PRJEB81106, original portable run | Both ranks, 16 stages; 188.03s with eight FastSpar workers | Existing-server execution, including fresh DKI training |
| PRJEB81106, independent clone/new analysis environments | Commit `9fbb1c6647b2f3b6fe7e2a60b9528a6c51446a82`; 16 stages, 288.16s, four workers | Completed after disclosed R-version and local C++ library handling; not an unmodified-installation pass |
| Installation-fix baseline | [`9f4fce74c0d04716903c25caff612ffd1661c699`](https://github.com/2472890096-cmd/MACK-DB/tree/9f4fce74c0d04716903c25caff612ffd1661c699) | Same scientific sources/inputs; SciPy source and R version-check instructions corrected. This was not a newly repeated full installation/run |
| Current-manuscript processed-input checks, 12 November 2026 | See [per-result execution records](docs/PAPER_REPRODUCTION.md#recorded-validation) | Real statistics were recalculated from existing associations, candidates and stored simulations; no cohort-wide inference or new random simulations. Source/data access limitations remain |

Get the exact checked-out version with `git rev-parse HEAD`. Before submission, freeze the author-approved commit/release and record that actual identifier and its matching validation evidence in the manuscript/data archive. A branch name is not an immutable version; no new final commit, tag or DOI is implied by this documentation revision.

[TESTING](TESTING.md) describes nine small scientific groups, twelve upstream safety cases and seven matrix-workflow I/O cases. Their passes do not substitute for full paper-wide input availability or execution. The known CXXABI library issue in one independent environment remains a platform-specific installation limitation, separate from the two corrected installation instructions.

## Code and data guide

| Location | Purpose |
|---|---|
| `pipeline/` | Profile/matrices, associations, Zi-Pi, EPI v3, unified-R DKI and environment support |
| `analysis/keystone_concordance/` | Broad/Strict three-method candidate statistics |
| `analysis/external_validation/` | Literature mapping, source exclusion, null, recurrence and support analyses |
| `examples/`, `tests/` | Real examples, synthetic inputs and small checks |
| `config/`, `lib/`, `tools/` | Nonsecret configuration, path handling and required adapters |

Scientific versions remain **EPI v3** and **DKI `MACKDB_DKI_R_20261012_v1`**. Genus abundance uses the author's [original species-label token grouping](docs/GENUS_AGGREGATION_STATUS.md); species-association regrouping for Figure 4 is a distinct manuscript operation. See [PIPELINE](PIPELINE.md), [core entries](docs/CORE_WORKFLOW.md), [core script reference](docs/SCRIPT_REFERENCE.md) and [current source manifest](docs/RELEASE_SCRIPT_MANIFEST.tsv), [parameters](PARAMETERS.md), [installation](INSTALL.md), [external dependencies](EXTERNAL_DEPENDENCIES.md), [third-party software](THIRD_PARTY_SOFTWARE.md), [data availability](DATA_AVAILABILITY_GUIDE.md) and [remaining release work](docs/FOLLOW_UP_ITEMS.md).

Large reads, whole-cohort matrices/networks, Sketches, MarkerDBs, predictions and database backups are not GitHub payloads. Supporting-data download/reviewer-access identifiers and source/data licensing must be supplied for the final release; no unavailable data DOI or license is invented here.
