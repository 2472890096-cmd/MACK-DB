# Immutable source provenance

ACCEPTED_84_SCRIPTS_ORIGINAL.tsv is the exact accepted v0.2 list (35main/45supplementary/4test), retained as historical evidence. BASELINE_93_SCRIPTS.tsv and BASELINE_151_TRACKED_FILES.tsv record the original Git commit, hash and review class. These old counts are not current release-file expectations.

ORIGINAL_93_DISPOSITION.tsv accounts for every original script, including17 removals from the default tree and current paths. Current hashes are checked by docs/RELEASE_SCRIPT_MANIFEST.tsv and docs/RELEASE_FILE_MANIFEST.tsv. The guard pins the immutable historical manifests separately and continues to verify aliases, tests, scopes and frozen-input boundaries.

The original Git repository, commit objects and external historical evidence are preserved. The removed script files are not distributed in this working tree; their original bytes remain in Git history and the separate local phase-one archive. No original source is overwritten.
