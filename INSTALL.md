# Requirements and environment setup

The supported portable example targets **Linux x86_64**, Python ≥3.11 for orchestration and existing R/external scientific tools. See the [portable-matrix setup guide](environments/portable-matrix/README.md) for environment YAMLs, an observed R environment and supplemental packages, FastSpar 1.0.0, the pinned external fastCCLasso core, and separate Python/torch training dependencies. No dependency is installed automatically. Recipes are provided for reproducibility; clean installation on a second server has not yet been tested.

```bash
python tools/verify_release_candidate.py
python examples/real_bioproject/PRJEB81106/run_workflow.py --steps all --preflight
```

Preflight checks only dependencies required by selected steps. `--steps epi --preflight` does not require FastSpar or PyTorch. Full `all` requires Python/numpy/pandas, R packages, the three FastSpar binaries, GNU parallel, the pinned fastCCLasso R file, and a separate compatible DKI Python with torch/torchdiffeq. Paths use PATH, explicit CLI arguments or [a nonsecret JSON template](config/portable_matrix.example.json). Use ignored local configuration; never add credentials.

The [bundled PRJEB81106 inputs](examples/real_bioproject/PRJEB81106/bundled/README.md) provide all real matrix/mapping and small taxonomy inputs for the default example. They do not include raw reads, model files or predictions. The original species-to-genus converter is included byte-identically. EPI/DKI cores and the preceding 79 accepted scripts remain unchanged.

R execution evidence includes R4.5.3, vegan2.7-5, microeco2.0.0 and igraph2.1.4; [environment records](environments/SOFTWARE_AND_PACKAGE_VERSIONS.tsv) distinguish historical production, independent reproduction and current observations. Do not infer a complete historical lockfile from these records, or switch microeco versions merely to pass tests.

The separate [ERR14036375 FASTQ example](examples/real_fastq/ERR14036375/README.md) requires the configured official MetaKSSD2.21, pinned L3K11/R226 files, GNU gzip/zcat, `/usr/bin/time`, and R tidyverse/optparse. Its binary, raw FASTQs, shuffle and MarkerDB remain external. SRA Toolkit/AWS are not dependencies of the accepted paired-FASTQ route.

The candidate contains 82 scripts: the preceding79 plus the original converter, portable complete-project orchestrator and its I/O tests. Original84 provenance remains immutable. Fourteen unrelated final-statistics sources remain held for individual rights confirmation; the full manuscript input/archive and project LICENSE remain separate release decisions. See [data](DATA_AVAILABILITY_GUIDE.md), [testing](TESTING.md), and [workflow](PIPELINE.md).
