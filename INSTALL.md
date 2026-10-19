# Environment and dependency setup — v0.2

No installation or upgrade was performed. [EXTERNAL_DEPENDENCIES.md](EXTERNAL_DEPENDENCIES.md) supplies verified official acquisition routes, version pins, licenses and configuration examples for a future isolated environment; those installation/build commands were not executed.

## Code-only private checkout

Python>=3.10 is required by Figure3 int.bit_count(). `python tools/verify_private_checkout.py` checks source/scopes/reference integrity without running science. Run the nine synthetic groups using the explicit existing R interpreter in TESTING.md. Real Figure3 candidate inputs are held pending authorization, so full-data Figure3 reproduction is not a self-contained entry in this commit.

## Small-test environment

The archived v0.1 executed suite used Python 3.13.5, R 4.5.3, vegan 2.7.5, microeco 2.0.0, igraph 2.1.4, tidyverse 2.0.0, optparse 1.7.5, dplyr 1.2.0, and data.table 1.17.8. The role test explicitly requires microeco 2.0.0 because its actual package implementation defines the tested equality/NA boundaries. Python >=3.10 is needed for Figure3 `int.bit_count()` and current annotations. OS support: Linux, Bash, `/usr/bin/time`, POSIX locks/hard links for the EPI worker.

```bash
python --version
Rscript --version
python tools/check_dependencies.py --rscript "$(command -v Rscript)"
python tests/run_small_tests.py --rscript "$(command -v Rscript)"
```

The test driver reports BLOCKED when an explicitly probed dependency is unavailable. Other execution errors remain FAIL, with logs. Do not silently substitute a different package implementation or upgrade a production environment to make a test pass. Historical environment evidence is summarized in `environments/SOFTWARE_AND_PACKAGE_VERSIONS.tsv`; it is not a claim that every historical production run used these versions. Archived version evidence remains in `environments/SOFTWARE_AND_PACKAGE_VERSIONS.tsv`.

## Other main-pipeline requirements

- Abundance/DKI input/name adapters: Python numpy and pandas. The tested base environment has numpy 2.1.3 and pandas 2.2.3.
- DKI model training/prediction, **not part of this test run**: torch, torchdiffeq, numpy. An existing separate environment was probed: Python 3.10.19, torch 2.7.1+cu118, torchdiffeq 0.2.5, numpy 2.2.6. The default Python has no torch/torchdiffeq; select the appropriate existing interpreter explicitly. No training or model prediction was run.
- Full Zi-Pi program additionally uses htmlwidgets, networkD3, gridExtra, optparse, tidyverse and plotting dependencies; only role calculations were tested. `corr_inter.R` additionally loads VennDiagram.
- MetaKSSD composite needs an external executable and the exact GTDB R226 marker Sketch DB plus a sample Sketch. A MetaKSSD executable is present on the current host, but these scientific inputs were not imported into this candidate.
- FastSpar requires `fastspar`, `fastspar_bootstrap`, `fastspar_pvalues`, GNU parallel. The three FastSpar executables are absent from the current PATH. No executable was installed or searched for across unrelated historical directories.
- fastCCLasso wrappers need the official `fastCCLasso.R` core, supplied explicitly through `--fast_core` or `FASTCCLASSO_CORE`. It is not vendored; its official retrieval command, current checkout/source SHA and LGPL license are now documented in EXTERNAL_DEPENDENCIES.md. The wrappers use optparse/parallel; **future.apply is not a wrapper dependency** even though the optional probe reports it missing.

Copy `config/runtime.env.example` to a private, ignored local configuration and set only paths to software/data you already possess. Do not put passwords or credentials into these files. `config/paths.example.json` maps optional supplementary inputs; relative values resolve from the repository root.

```bash
export MACKDB_RSCRIPT="$(command -v Rscript)"
export MACKDB_PYTHON="$(command -v python)"
export MACKDB_THREADS=1
# Set METAKSSD_BIN, FASTSPAR_BIN, FASTSPAR_BOOTSTRAP_BIN,
# FASTSPAR_PVALUES_BIN, FASTCCLASSO_CORE only when available locally.
```

External inference software and large upstream data must be obtained under their own terms. Exact acquisition/configuration is documented, but fresh installation and a clean container have not been tested. Only scope lists and aggregate Figure3 references are bundled; row-level inputs are held. Production-version gaps are explicit in EXTERNAL_DEPENDENCIES.md; no unsupported promise of a single fully locked historical environment is made.
