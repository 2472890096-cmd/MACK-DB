# External dependencies: acquisition, configuration and evidence

No software was installed/upgraded, no external scientific program was run and no AWS connection was made. The instructions below are for a future isolated reviewer environment. Source and license pages were checked on 2026-10-17. A documented install route is not a clean-install PASS.

## Version evidence

| Component | Confirmed historical production | Independent reproduction | Current local evidence |
|---|---|---|---|
| MetaKSSD | Exact production runtime version not recovered | 2.21, three-project archive | Local checkout c56806aed7b29b218c88c341dc6472545e39d8bb; src/metakssd.c reports 2.21; executable on PATH; not executed this round |
| FastSpar, bootstrap, pvalues | Do not infer a uniform production version | 1.0.0, three-project archive | Not found on current PATH; official v1.0.0 release exists |
| fastCCLasso | Logs identify official RCode/fastCCLasso.R and kmax20/nboot100 | Same actual code definition is the accepted provenance | Current local unmodified core SHA256 41004c99db932346c651447ef561d7965fddee0dc3252aad63ac6681b57a8638, checkout 7833b61b02723a34e307e6325ab9b48564debcf3; checkout date/commit is not a historical run snapshot |
| Zi-Pi R / igraph / microeco | R4.5.3, igraph2.1.4, microeco2.0.0 printed in formal project logs | Archived reproduction environment separately recorded | r_network metadata query gives the same versions |
| EPI / DKI dependencies | Some historical environment snapshots absent | EPI R4.5.3/vegan2.7-5; DKI Python3.10.19/torch2.7.1+cu118/torchdiffeq0.2.5 in recorded reproduction | Current package metadata is not retrospective proof |

Evidence: environments/SOFTWARE_AND_PACKAGE_VERSIONS.tsv preserves historical labels with internal source paths reduced to archive basenames. The source hashes/commits in this document were verified in accepted v0.2; full local probe/audit files are not committed. Frozen official-core source mapping is supported by the previously accepted inventory, not a claim that current Git HEAD proves every historical run.

## MetaKSSD

Official source: [yhg926/MetaKSSD](https://github.com/yhg926/MetaKSSD). The documented build uses make; the locally reviewed Makefile uses GCC GNU11, zlib, OpenMP and libm. The local LICENSE.txt is Apache-2.0. Obtain prerequisites in a new environment without changing a production installation.

From the MACK-DB repository root, these are optional acquisition/build commands, **not executed here**:

```bash
mkdir -p external
git clone https://github.com/yhg926/MetaKSSD.git external/MetaKSSD
git -C external/MetaKSSD checkout --detach c56806aed7b29b218c88c341dc6472545e39d8bb
make -C external/MetaKSSD
export METAKSSD_BIN="$PWD/external/MetaKSSD/bin/metakssd"
"$METAKSSD_BIN" --version
```

The pin identifies the current inspected 2.21 source, not a recovered production binary. Verify the executable version and archive its SHA before use. The official README links [GTDB R226 MarkerDB](https://zenodo.org/records/16317275); database/sample Sketch inputs are separate and not bundled. Follow PIPELINE.md for the existing composite query and MACK-DB Profile sharekmer12; do not replace its Profile stage with the upstream generic default18.

## FastSpar 1.0.0

Official [v1.0.0 release](https://github.com/scwatts/fastspar/releases/tag/v1.0.0) and [versioned build instructions](https://raw.githubusercontent.com/scwatts/fastspar/v1.0.0/README.md). Required source-build tools include C++11/OpenMP, Fortran, Armadillo, LAPACK/OpenBLAS, GSL, make and autotools. License: GPL-3.0. Use a private prefix with the prerequisites already available:

```bash
mkdir -p external
git clone --branch v1.0.0 --depth 1 https://github.com/scwatts/fastspar.git external/fastspar
MACKDB_FASTSPAR_PREFIX="$PWD/external/fastspar-install"
(
  cd external/fastspar
  ./autogen.sh
  ./configure --prefix="$MACKDB_FASTSPAR_PREFIX"
  make
  make install
)
export FASTSPAR_BIN="$MACKDB_FASTSPAR_PREFIX/bin/fastspar"
export FASTSPAR_BOOTSTRAP_BIN="$MACKDB_FASTSPAR_PREFIX/bin/fastspar_bootstrap"
export FASTSPAR_PVALUES_BIN="$MACKDB_FASTSPAR_PREFIX/bin/fastspar_pvalues"
export PARALLEL_BIN="$(command -v parallel)"
"$FASTSPAR_BIN" --version
"$FASTSPAR_BOOTSTRAP_BIN" --version
"$FASTSPAR_PVALUES_BIN" --version
```

The official alternative is conda; in a new environment request `fastspar=1.0.0` from conda-forge/bioconda and verify all three binaries. Availability/solver success has not been tested here. Full inference is the real `pipeline/associations/fastspar/run_fastspar.sh COUNT_MATRIX OUTPUT_PREFIX` entry in PIPELINE.md, with 10000 bootstraps/permutations and final abs(r)>=0.3, P<=1e-4. It was not executed. GNU parallel is a separate prerequisite.

## fastCCLasso official R core

Official [ShenZhang-Statistics/fastCCLasso](https://github.com/ShenZhang-Statistics/fastCCLasso), actual path **RCode/fastCCLasso.R** (case matters). The repository's [GNU LGPL-3.txt](https://github.com/ShenZhang-Statistics/fastCCLasso/blob/main/GNU%20LGPL-3.txt) is LGPLv3. The inspected core is 32612 bytes; it is not vendored into MACK-DB.

```bash
mkdir -p external
git clone https://github.com/ShenZhang-Statistics/fastCCLasso.git external/fastCCLasso
git -C external/fastCCLasso checkout --detach 7833b61b02723a34e307e6325ab9b48564debcf3
export FASTCCLASSO_CORE="$PWD/external/fastCCLasso/RCode/fastCCLasso.R"
printf '%s  %s\n' '41004c99db932346c651447ef561d7965fddee0dc3252aad63ac6681b57a8638' "$FASTCCLASSO_CORE" | sha256sum -c -
Rscript --vanilla -e 'invisible(parse(file=Sys.getenv("FASTCCLASSO_CORE"))); cat("PARSE_ONLY_PASS\n")'
```

This parse check does not run inference. The official fastCCLasso path uses base/stats R functions; wrappers additionally use optparse and parallel. No requirement for future.apply is introduced. Pass the checked core through `--fast_core "$FASTCCLASSO_CORE"` in the species/genus commands in PIPELINE.md. These retain nboot100/kmax20, accepted seeds and the final refilter.

The actual P definition is **not a bootstrap sign-frequency test**. Let rbar be the mean of the 100 bootstrap estimated correlations **plus the original fitted correlation** (101 columns). The code computes `t = rbar * sqrt((n-2)/(1-rbar^2))`, `u = pt(t, df=n-2)`, then `P = min(u, 1-u)` elementwise, without multiplying by 2. It symmetrizes the matrix, with diagonal 1. Bootstrap estimates the mean correlation; Student-t tails provide the returned P. Do not describe it as a conventional doubled two-sided P. The v0.1 phrase “bootstrap sign tails” is corrected only in documentation; no code or parameters changed. Direct source evidence: boot_fastCCLasso in the pinned, hashed official file; wrappers read result$p_vals unchanged. No inference/100-bootstrap validation was rerun.

## Configuration and read-only checks

Keep external software under ignored external/ or another explicit user-owned path. Do not commit installations or credentials. The configuration template contains paths only. Export variables when invoking child processes:

```bash
export MACKDB_PYTHON="$(command -v python3)"
export MACKDB_RSCRIPT="$(command -v Rscript)"
export MACKDB_THREADS=1
python tools/check_dependencies.py --rscript "$MACKDB_RSCRIPT"
```

Use a compatible existing R environment, not whichever Rscript happens to be first on PATH. The dependency checker reads package/executable availability; its optional missing-package entries are not automatically scientific blockers. No fresh-machine or clean-install execution was authorized. Figure3 synthetic set tests are independent of these inference dependencies. Real candidate inputs are separately held pending upload authorization.
