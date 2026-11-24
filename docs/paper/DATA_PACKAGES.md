# Paper input packages and access

These local review candidates contain processed research results. **None has been uploaded or assigned a public data URL/DOI.** Their availability is not established by the code repository being public. Review rights, third-party redistribution, metadata and licenses before deposition. No raw reads, models, Sketches, MarkerDBs, production database backup or credentials are proposed for these packages.

| Package | Prepared contents | Raw bytes | Compressed bytes | Current state |
|---|---|---:|---:|---|
| `Figure2_processed_environment_associations_20261112.tar.gz` |180 associations plus90-row summary and120-row display reference;182 files|1724945736|507859601|Extracted and SHA-verified; data-release review pending|
| `MACKDB_PAPER_FIGURE3_FINAL_INPUTS_20261112.tar.gz` |Six FINAL unified-R candidates,61-row assertion,scopes,aggregate/support references;15 files|2748328|480745|Extracted and SHA-verified; held row-level data not authorized for Git upload|
| `MACKDB_PAPER_EXTERNAL_NULL_SUPPORT_INPUTS_20261112.tar.gz` |22 stored1000-iteration outcomes,configs,source-filtered pair records,tested sets,653071-pair rank reference,summary oracles;53 files|242351521|49184060|Local prepared candidate; internal provenance fields need review; no new simulations|
| `MACKDB_PAPER_FIGURE4_REVIEW_INPUTS_20261112.tar.gz` |149 data/provenance files for Figure4 statistics;4 source/provenance scripts deliberately excluded|360433995|276474987|**Partial review package**; third-party workbooks and metadata need review; remaining source/data access must be checked individually|

Archive SHA256 values:

```text
4093c398bc8ae79cdfdab3e64ebfb19c41a65aff0b093b5ef5ccb3b312fcd0f2  Figure2_processed_environment_associations_20261112.tar.gz
075d197b784022235e297e08860b812a42f1343626136660bbf08c575222655e  MACKDB_PAPER_FIGURE3_FINAL_INPUTS_20261112.tar.gz
9264fed4a6408c42d141a5ca3ff46b01bb744768f3850ca5c9dee3d525312d23  MACKDB_PAPER_EXTERNAL_NULL_SUPPORT_INPUTS_20261112.tar.gz
0e16e76c999841fa082dfe7b14772b869b0b6bf7ec33235782a1e1537ebff46d  MACKDB_PAPER_FIGURE4_REVIEW_INPUTS_20261112.tar.gz
```

[DATA_PACKAGE_MANIFEST.tsv](DATA_PACKAGE_MANIFEST.tsv) lists399 logical data paths, sizes, SHA256 and access statuses. This count describes these four packages, not399 files to commit to Git. The Figure4 executed-input audit has153 items; four code/provenance sources are access dependencies, not data payload. No claim of independent anonymity review is made for the unpublished large packages.

## Layout and checking

Figure2 extracts with `figure2/` below the selected data root; see [its input manifest](../../analysis/associations/INPUT_DATA_MANIFEST.tsv). Figure3 extracts with `figure3/`; see [FIGURE3_INPUT_MANIFEST.tsv](FIGURE3_INPUT_MANIFEST.tsv). The external package has `literature/`, `chen_wallen/` and `rank_reference/` directly at its root; extract it **into** `data/paper-inputs/`; see [its manifest](../../analysis/paper_support/INPUT_ARCHIVE_MANIFEST.tsv). These layouts match the implemented replay commands in [PAPER_REPRODUCTION](../PAPER_REPRODUCTION.md).

Figure4 extracts below `figure4/analysis_inputs` and `figure4/project_inputs`. Its runtime configuration must map the original frozen path prefixes to those local roots. The local149-file archive preserves original provenance bytes; it is **not a sanitized ready-to-publish replacement** for all153 dependencies. The validated relocated144-file input setup and complete local source manifest are retained outside Git. Before distribution, review the three publication workbooks, sample/participant-named mapping records and frozen metadata, and resolve the remaining source/data dependencies individually. The required Figure4 statistics sources are now included in the code tree; they are not hidden inside the data archive. Public Run accession identifiers are separately provided in [Figure4B_sequencing_runs.tsv](../../analysis/paper_network_comparison/reference/Figure4B_sequencing_runs.tsv).

Small code-repository references and synthetic tests do not replace these real input packages. Original third-party tables should be obtained from their cited publishers when redistribution permission is absent; a freely downloadable paper is not automatically a data redistribution license. Arial fonts and held original source code remain outside all four packages.

For an approved local download, set `PAPER_PACKAGE_DIR` to the directory containing the exact checked archives. Run from the repository root; use a new empty data destination to avoid replacing earlier files:

```bash
mkdir -p data/paper-inputs
tar -xzf "$PAPER_PACKAGE_DIR/Figure2_processed_environment_associations_20261112.tar.gz" -C data
tar -xzf "$PAPER_PACKAGE_DIR/MACKDB_PAPER_FIGURE3_FINAL_INPUTS_20261112.tar.gz" -C data
tar -xzf "$PAPER_PACKAGE_DIR/MACKDB_PAPER_EXTERNAL_NULL_SUPPORT_INPUTS_20261112.tar.gz" -C data/paper-inputs
```

The Figure4 partial review archive is not advertised as a ready-to-run download: metadata review, third-party redistribution and finalized input-prefix mapping must first be resolved. The required B/C/D statistics sources are now bundled with the code. No example here downloads an unpublished URL.

## Submission next step

The author must approve the data release, resolve source/font rights, provide reviewer access or deposit the reviewed processed inputs, and record actual persistent download links. Then record the matching immutable Git commit. No GigaDB/Zenodo identifier or formal code Release is assumed here.
