# Runnable examples

- **Recommended complete-project start:** [PRJEB81106 matrices](real_bioproject/PRJEB81106/README.md). Bundled original 53×33 species and 38×33 genus matrices, names, Run mapping and small references. `run_workflow.py` branches to FastSpar→Zi-Pi, independent fastCCLasso/EPI, and new DKI input→training→unified-R scoring. No FASTQ download or archived predictions are needed for this example.
- **Real upstream start:** [ERR14036375 paired FASTQ](real_fastq/ERR14036375/README.md). Official MetaKSSD2.21→count Sketch→GTDB R226 Query→Species Profile was executed and Profile matches the archive byte-for-byte. Raw data and reference databases are obtained separately.
- [Archived 33-Run identities](prjeb81106/archived_runs.txt) fix the project example cohort. Current BioProject membership is larger and must not replace these Runs.
- [Earlier matrix plan](matrix/README.md) and [earlier full-project FASTQ plan](fastq/README.md) are historical planning documents. The old `run_stage.py` is retained for archived-input/scoring replay compatibility; it is not the default portable route.

A single-Run Profile is not a substitute for the complete project matrix required by the downstream methods. The matrix example uses the author's adopted species-label genus-token aggregation, which folds some GTDB suffix genera. It does not claim fresh upstream processing for all33Runs, full-cohort reconstruction, or a successful environment install on every machine. Dependencies and actual evidence are listed in each guide.
