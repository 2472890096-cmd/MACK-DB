# Portable Linux matrix environment

Target: **Linux x86_64**, Bash, writable local work space, `/usr/bin/time`, GNU core utilities and gzip. This directory separates the Python orchestrator, R analysis runtime and FastSpar tools. It does not promise support for every operating system or CPU architecture. No installer, solver, package download, neural training or association inference was run while preparing these files.

The real PRJEB81106 matrix example uses Python **>=3.11**, including compatibility with the earlier `run_stage.py` (`hashlib.file_digest`). Its [bundled inputs](../../examples/real_bioproject/PRJEB81106/bundled/README.md) now contain the complete 53-taxon/33-sample species matrix, archived 38-taxon/33-sample genus matrix, ordered Run/sample mappings, minimal matching taxonomy references and small expected outputs. No separate archived matrix input download is needed. The independent one-Run FASTQ example does not replace these project matrices. EPI and DKI start from abundance independently of FastSpar; Zi-Pi starts from a network. The current portable entry is [run_workflow.py](../../examples/real_bioproject/PRJEB81106/run_workflow.py); see [the matrix example](../../examples/real_bioproject/PRJEB81106/README.md).

## 1. Python and FastSpar

The following are **installation instructions for a new, isolated environment**, not commands already executed during this release preparation. Conda must already be installed. The YAML files pin direct packages but are not complete solver locks; stop and report a solver error rather than silently changing scientific tool versions.

```bash
conda env create --file environments/portable-matrix/matrix-python.yml
conda env create --file environments/portable-matrix/fastspar.yml
conda activate mackdb-matrix-python
export MACKDB_PYTHON="$(command -v python)"
```

Select the FastSpar tool environment and retain its explicit binary directory:

```bash
conda activate mackdb-fastspar-1.0.0
export FASTSPAR_DIR="$CONDA_PREFIX/bin"
export PATH="$FASTSPAR_DIR:$PATH"
"$FASTSPAR_DIR/fastspar" --version
"$FASTSPAR_DIR/fastspar_bootstrap" --version
"$FASTSPAR_DIR/fastspar_pvalues" --version
parallel --version
```

All three FastSpar programs must report **1.0.0**. Its official [versioned instructions](https://raw.githubusercontent.com/scwatts/fastspar/v1.0.0/README.md) document Conda and source builds. The formal example retains 10,000 bootstrap/permutation replicates. GNU parallel must be on PATH during the example. The exact local binary/build and timing record is separate from this untested installation recipe.

## 2. R runtime

The R runtime observed for the accepted small tests and bounded matrix runs is **R 4.5.3**, with **microeco 2.0.0**, **igraph 2.1.4**, **vegan 2.7-5**, tidyverse 2.0.0, optparse 1.7.5, dplyr 1.2.0 and data.table 1.17.8. These are recorded runtime versions, not proof that every historical production project used one identical environment. [Version evidence](../SOFTWARE_AND_PACKAGE_VERSIONS.tsv) preserves that distinction.

- [r-direct-requirements.tsv](r-direct-requirements.tsv): 14 packages directly required by the retained matrix workflow and its plot/export chain. Plot imports must still load because the unchanged Zi-Pi entry imports them at startup.
- [r-observed-linux64-explicit.txt](r-observed-linux64-explicit.txt): actual 305-package Conda metadata snapshot with package URLs and MD5 values. This is a reproducible reference for the **R tool environment**, not a minimal package list or a tested fresh installation. It contains Python 3.14.3; keep the separate Python 3.11 orchestrator above.
- [r-observed-packages.tsv](r-observed-packages.tsv): observed 196 R packages and provenance metadata; library paths and credentials are omitted.
- [r-required-supplemental.tsv](r-required-supplemental.tsv): eight required R packages outside the Conda package snapshot, including their version/source metadata.
- [r-supplemental-packages.tsv](r-supplemental-packages.tsv): all 43 non-Conda, non-base packages observed locally. Most development packages here are **not** matrix workflow dependencies.
- [r-optional-export.tsv](r-optional-export.tsv): rgexf 0.16.3 for GEXF export; the unchanged Zi-Pi entry catches an export error and writes GraphML instead. Pandoc is also needed for self-contained HTML output and is present in the observed Conda reference.

To reconstruct the observed R reference in a new environment:

```bash
conda create --name mackdb-matrix-r --file environments/portable-matrix/r-observed-linux64-explicit.txt
conda activate mackdb-matrix-r
export MACKDB_RSCRIPT="$CONDA_PREFIX/bin/Rscript"
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
"$MACKDB_RSCRIPT" --version
```

Then install only the required non-Conda supplements into that **new** environment. The bootstrap helper `remotes` is an installation tool. The [official remotes documentation](https://cran.r-project.org/web/packages/remotes/refman/remotes.html) documents versioned CRAN and commit-specific GitHub installation. This source installation needs an appropriate compiler/system libraries from the reference environment and access to upstream archives. Network retrieval and availability have not been tested here.

```bash
"$MACKDB_RSCRIPT" --vanilla -e 'install.packages("https://cran.r-project.org/src/contrib/remotes_2.5.0.tar.gz", repos=NULL, type="source", lib=Sys.getenv("R_LIBS_USER"))'
"$MACKDB_RSCRIPT" --vanilla - <<'RS'
lib <- Sys.getenv("R_LIBS_USER")
stopifnot(nzchar(lib), getRversion() == "4.5.3")
pins <- c(permute="0.9-10", ape="5.8-1", plyr="1.8.9", reshape2="1.4.5", vegan="2.7-5", ggrepel="0.9.7", microeco="2.0.0")
for (p in names(pins)) {
  remotes::install_version(p, version=pins[[p]], repos="https://cloud.r-project.org", lib=lib, dependencies=NA, upgrade="never")
  stopifnot(as.character(utils::packageVersion(p)) == pins[[p]])
}
remotes::install_github("mattflor/chorddiag@1688d72cd93071abb373e054190363bdfb3af2af", lib=lib, dependencies=NA, upgrade="never")
stopifnot(as.character(utils::packageVersion("chorddiag")) == "0.1.3")
RS
```

If a package reports a missing dependency or upstream archive, stop and preserve the message. Do not install a different microeco or igraph to make a numerical oracle pass. No R package source has been vendored here. Final plot output and scientific node-table validation remain distinct.

## 3. fastCCLasso official core

The upstream LGPL-3 source is obtained separately and is **not included** in this repository. The following retrieves the exact reviewed commit and verifies the core, without relying on a production directory. Only execute acquisition commands when preparing your own environment:

```bash
mkdir -p external
git clone https://github.com/ShenZhang-Statistics/fastCCLasso.git external/fastCCLasso
git -C external/fastCCLasso checkout --detach 7833b61b02723a34e307e6325ab9b48564debcf3
export FASTCCLASSO_CORE="$PWD/external/fastCCLasso/RCode/fastCCLasso.R"
printf '%s  %s\n' '41004c99db932346c651447ef561d7965fddee0dc3252aad63ac6681b57a8638' "$FASTCCLASSO_CORE" | sha256sum --check -
```

Pass `--fastcclasso-core "$FASTCCLASSO_CORE"` to the real-matrix example. The unchanged formal wrapper uses nboot=100 and kmax=20. Its actual P definition is the smaller Student-t tail of the averaged original plus bootstrap correlations, **without doubling**; see [dependency details](../../EXTERNAL_DEPENDENCIES.md). A successful source hash check is not an inference run.

## 4. Fresh DKI training versus the earlier replay

The new `run_workflow.py --steps all` entry includes **fresh perturbation generation, training, neural prediction and unified-R scoring**. Its bundle contains no saved `qtrn`/`qtst` or checkpoint, so these cannot silently replace training. Training requires torch and torchdiffeq even when a GPU is absent; the example explicitly uses CPU. EPI and Zi-Pi can be requested without the training environment.

The earlier `run_stage.py --stage dki-replay` is a separate archived-prediction replay path: it needs matched saved Ptrain/Ptest/indices/qtrn/qtst, which are not supplied by the new bundle. Replay does not train or generate predictions and must not be used as the evidence for the new training path.

[dki-training-optional.yml](dki-training-optional.yml) describes the separate model environment. It is optional only when selecting steps that exclude `dki-train`/`dki-score`; it is **required for `--steps all`**. The archived model runtime was Python 3.10.19, torch 2.7.1+cu118 and torchdiffeq 0.2.5; these are reproduction records, not a newly tested universal installation. The [official PyTorch previous-version page](https://pytorch.org/get-started/previous-versions/) supplies 2.7.1 CPU and CUDA wheels; [torchdiffeq](https://github.com/rtqichen/torchdiffeq) documents its PyTorch dependency.

Optional CPU environment preparation, **not executed here**:

```bash
conda env create --file environments/portable-matrix/dki-training-optional.yml
conda activate mackdb-dki-training
python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cpu
python -m pip install torchdiffeq==0.2.5
python -c 'import torch,torchdiffeq; print(torch.__version__,torchdiffeq.__version__)'
export MACKDB_DKI_PYTHON="$(command -v python)"
```

This CPU wheel is a documented alternative build to the archived cu118 wheel, not a claim of bitwise execution equivalence. A new model's actual execution and scientific checks must be read from its `dki-train` and `dki-score` status/log files; installation instructions alone establish no training PASS. Do not describe saved-prediction replay as full training reproducibility. Pass this model interpreter via `--dki-python`; the matrix orchestrator remains Python>=3.11.

## 5. Check before running

After activating/configuring the correct environments, perform read-only dependency probes from the repository root. The output is written to a user-selected work directory, not to source directories:

```bash
export PATH="$FASTSPAR_DIR:$PATH"
"$MACKDB_PYTHON" tools/check_dependencies.py --rscript "$MACKDB_RSCRIPT" --output work/environment/dependency_probe.json
"$MACKDB_RSCRIPT" --vanilla -e 'p<-c("vegan","microeco","igraph","tidyverse","optparse","dplyr","data.table","htmlwidgets","networkD3","gridExtra","ggrepel","RColorBrewer","VennDiagram","chorddiag");for(x in p){stopifnot(requireNamespace(x,quietly=TRUE));cat(x,as.character(packageVersion(x)),"\n")}'
"$MACKDB_PYTHON" tests/run_small_tests.py --rscript "$MACKDB_RSCRIPT"
```

The final command is a bounded scientific test run, not an environment probe. On a new server, record its actual result; an environment recipe alone cannot be marked PASS. The existing nine-group results are local accepted evidence only. No automatic software installation is performed by the scientific runner.

## 6. Complete-project entry with bundled matrices

Run from the repository root after configuring the absolute executable variables above. `--preflight` checks the input hashes, scientific source hashes and tools required by the selected steps without starting those analyses:

```bash
EXAMPLE=examples/real_bioproject/PRJEB81106/run_workflow.py
export PATH="$FASTSPAR_DIR:$PATH"
"$MACKDB_PYTHON" "$EXAMPLE" --steps all --preflight \
  --rscript "$MACKDB_RSCRIPT" --dki-python "$MACKDB_DKI_PYTHON" \
  --fastspar-bin-dir "$FASTSPAR_DIR" --fastcclasso-core "$FASTCCLASSO_CORE"
```

The following explicitly requests all branches for both species and genus, including **new CPU training and prediction** for each rank, and writes only a new work directory:

```bash
"$MACKDB_PYTHON" "$EXAMPLE" --steps all --levels species genus \
  --output work/PRJEB81106_complete --threads 4 --timeout-seconds 600 \
  --rscript "$MACKDB_RSCRIPT" --dki-python "$MACKDB_DKI_PYTHON" \
  --fastspar-bin-dir "$FASTSPAR_DIR" --fastcclasso-core "$FASTCCLASSO_CORE"
```

Each stage is bounded at 600 seconds in this command. A timeout stops its child process group, preserves failure logs and does not automatically retry, reduce the formal scientific parameters or switch to saved predictions. This is the bounded local acceptance configuration, not a guarantee that every server can finish each stage within 600 seconds. An existing output root is rejected unless `--resume` is explicitly selected; resume validates completed artifacts and input/code/tool identities before reuse. Use a new output root for an intentionally separate run.

For an initial abundance-only branch without installing FastSpar or torch:

```bash
"$MACKDB_PYTHON" "$EXAMPLE" --steps epi --preflight --rscript "$MACKDB_RSCRIPT"
"$MACKDB_PYTHON" "$EXAMPLE" --steps epi --output work/PRJEB81106_epi --rscript "$MACKDB_RSCRIPT"
```

`--steps dki-score` automatically includes `prepare`, `dki-input` and **fresh** `dki-train`; `--steps zipi` includes `prepare` and FastSpar. `--steps all` adds fastCCLasso as an independent network-inference branch. The entry also accepts [config/portable_matrix.example.json](../../config/portable_matrix.example.json) through `--config`; copy/edit a local configuration if preferred. The template omits interpreter keys so the exported `MACKDB_PYTHON`, `MACKDB_RSCRIPT` and `MACKDB_DKI_PYTHON` remain active. Adding those keys as absolute paths overrides the environment; do not set both to an ambiguous `python3` after switching environments. Always launch the orchestrator with `"$MACKDB_PYTHON"`, not the currently active model Python. Bare executable names resolve on PATH, configuration paths resolve relative to the JSON file, and explicit CLI tool flags override them. No hidden Dell path is needed. Default invocation is a prerequisite check only, not a scientific run.

Successful stages produce `SUCCESS.json`, command logs and artifact checksums; failures preserve `FAILED.json` and `WORKFLOW_STATUS.json`. Expected archival results are references for comparison, not inputs that substitute for a scientific computation. Fresh stochastic model results are not required to equal historical predictions merely because the seed is fixed. Keep actual runtime completion, output checks and historical-reference equality as separate findings.

Documentation and official-source checks: 2026-10-21. Unverified by this environment-documentation work: fresh Conda solve/download, reconstruction from these exact explicit URLs, source compilation and clean-server execution. Actual bounded training/inference acceptance belongs to the workflow execution records, not this environment recipe. The bundled matrices remove the earlier missing-example-input requirement; public distribution still needs the author-approved source/data license and retained GTDB attribution.
