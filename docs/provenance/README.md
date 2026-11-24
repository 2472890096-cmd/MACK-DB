# Source provenance and current layout

ACCEPTED_84_SCRIPTS_ORIGINAL.tsv is the exact accepted v0.2 list (35main/45supplementary/4test), retained as immutable historical evidence. BASELINE_93_SCRIPTS.tsv and BASELINE_151_TRACKED_FILES.tsv retain the original Git commit, hashes and review classes. These three files are unchanged; their old counts are not current release-file expectations.

ORIGINAL_93_DISPOSITION.tsv accounts for every original script: the earlier17 optional removals plus six author-approved historical/extension removals on24November2026; retained sources, path adaptations and shared implementations have explicit current destinations. Deduplicated aliases resolve to a shared implementation rather than an absent scientific method.

SOURCE_RESTRUCTURE_20261124.tsv records old/current paths and SHA256, plus preserved Python function-body fingerprints for approved path-only migrations. The EPI R core remains byte-identical. Current script/file hashes are in docs/RELEASE_SCRIPT_MANIFEST.tsv and docs/RELEASE_FILE_MANIFEST.tsv. The guard checks those hashes, the shared EPI layout, historical source identities, retained function bodies, tests, scopes and frozen-input boundaries; it is not disabled by the restructuring.

Removed scripts remain available in Git history; this working tree does not carry another history-code folder. Original production sources and scientific results outside the repository are untouched. Intentional future changes must update the versioned current manifests and their declared provenance rather than bypassing the guard.
