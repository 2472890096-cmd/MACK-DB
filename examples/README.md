# Real examples

| Start from | Example | Included input | What runs |
|---|---|---|---|
| Complete abundance matrices | [PRJEB81106](real_bioproject/PRJEB81106/README.md) | Archived33-Run species53×33 and genus38×33 matrices, mappings and small references | FastSpar→Zi-Pi; independent fastCCLasso/EPI; fresh DKI training, prediction and unified-R scoring |
| One paired-FASTQ Run | [ERR14036375](real_fastq/ERR14036375/README.md) | Download instructions, input pins and small expected Profile/Query references | Official MetaKSSD→GTDB R226 Query→Species Profile |

The matrix workflow needs no FASTQ download or saved model/predictions. The upstream example requires separately obtained reads and reference resources. A single Run cannot replace the complete matrix required by downstream project analyses.

Each guide gives installation, commands, input/output formats, expected results and validation limits. Neither example repeats all1,575 BioProjects. Paper-wide statistics use the separate [paper reproduction route](../docs/PAPER_REPRODUCTION.md).

Superseded planning files and the earlier private-path/scoring-replay example were retired from this candidate. Their source remains in Git history; they are not required by these two entries.
