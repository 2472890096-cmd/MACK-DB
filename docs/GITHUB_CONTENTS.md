# Repository organization

The public code tree serves two routes: [run the real examples](../examples/README.md) and [recompute paper statistics](PAPER_REPRODUCTION.md). Helpers and tests are dependencies, not additional workflows a reader must run separately.

## Keep in GitHub

| Location | Required content |
|---|---|
| `pipeline/` | Accepted scientific implementations and their runtime wrappers |
| `analysis/` | Actual statistical calculations, necessary validation and optional rendering adapters |
| `examples/` | One complete-project matrix entry and one single-Run upstream entry, with small inputs/references |
| `environments/`, `config/` | Dependencies, pinned external versions and portable configuration templates |
| `lib/`, `tools/` | Imports, matrix adapters and integrity/input checks used by the entries |
| `tests/` | Small scientific tests, synthetic fixtures and expected values |
| `docs/`, `frozen_inputs/` | Detailed instructions, file hashes, accepted project scopes and small numeric references |

The required Figure2/4 originals are kept in `analysis/associations/original/` and `analysis/paper_network_comparison/original/`. Run their documented adapters; the original files retain historical paths and may perform IO when run directly.

The [versioned source manifest](RELEASE_SCRIPT_MANIFEST.tsv) lists active scripts. [The core guide](CORE_WORKFLOW.md) explains scientific call order. EPI species/genus launchers use one shared core/worker; generated run snapshots are kept only under each isolated output root. Immutable source-provenance files required by the integrity guard remain under `docs/provenance/`; they are maintenance evidence, not user entry points.

## Keep outside GitHub

Large processed association/candidate inputs belong in separately reviewed [data packages](paper/DATA_PACKAGES.md), with access instructions and SHA256 in GitHub. Raw FASTQ/SRA, Sketches, MarkerDBs, models/predictions, database backups, credentials, caches, server logs and redundant historical reports are not code payloads. Third-party or individually held files are not bundled without authorization.

## Retired optional material

On24November2026 six optional historical/extension scripts were additionally removed from the default tree; their bytes remain in Git history. The Supplementary S1/S2 statistics and Figure4B calculations remain. Figure4B sources were moved from a misleading `logs/` directory into `analysis/paper_network_comparison/original/matched_subgroups/`.

The former `run_stage.py` replay entry, its private-input templates, superseded FASTQ/matrix plans and three unreferenced phase-one snapshots were removed from this candidate on16November2026 after checking callers. The current `run_workflow.py`, bundled data, scientific implementations, required provenance guards and paper-reproduction functions remain. Original bytes are preserved in Git history and a separate local backup. No scientific formula, threshold or accepted reference was changed.

The README is the short starting page. Historical timing details, all-script inventories and per-panel evidence stay in the linked specialist guides.
