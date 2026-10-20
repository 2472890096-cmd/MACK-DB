# MACK-DB

MACK-DB brings together microbial association networks and candidate keystone taxa across BioProjects and environments at species and genus levels. Its scientific workflow connects metagenomic abundance profiling, association inference, environment-level network integration and three complementary approaches to identifying candidate keystone taxa: Zi-Pi, EPI and DKI.

This repository contains the research code supporting the database and its associated analyses, including external validation and manuscript statistics. The database scope comprises 1,575 BioProjects; individual analyses use their documented input and project scopes. Large scientific inputs and the deployed website are maintained separately.

## Analysis workflow

| Module | Role in MACK-DB | Code |
| --- | --- | --- |
| MetaKSSD and Profile | Query sample Sketches against the marker reference, estimate abundance and assemble project abundance/count matrices. | [Profiling](pipeline/metakssd_profile/), [matrix preparation](pipeline/abundance/) |
| FastSpar | Infer associations from count matrices, calculate P values and apply the reviewed correlation/significance filters. | [FastSpar workflow](pipeline/associations/fastspar/) |
| fastCCLasso | Infer compositional associations using the external official R core and the reviewed bootstrap, P-value and final-filter implementation. | [fastCCLasso workflow](pipeline/associations/fastcclasso/) |
| Environment networks | Combine project-level association evidence and retain links according to the environment support rules. | [Network integration](pipeline/environment/) |
| Zi-Pi | Calculate within-module connectivity and participation coefficients, then identify Connector, Module hub and Network hub candidates. | [Zi-Pi workflow](pipeline/zipi/) |
| EPI v3 | Calculate D1, D2 and the corrected Q/ECDF metric with 25%–75% prevalence filtering; retain candidates supported by at least two metrics. Species and genus share the verified core. | [EPI v3 workflow](pipeline/epi_v3/), [single-project entry](tools/epi_project.py) |
| DKI, unified R scoring | Prepare perturbations, train/predict with cNODE2 when required, and score saved predictions using the focal-removed, renormalized `q_null` baseline. The accepted scoring version is **MACKDB_DKI_R_20261012_v1**. | [DKI workflow](pipeline/dki/) |
| External validation | Harmonize literature-network taxonomy, compare recovered associations, evaluate source exclusion and null models, and analyze recurrence and cross-study support. | [Validation analyses](analysis/external_validation/), [adoption and usage records](SUPPLEMENTARY_USAGE.tsv) |
| Manuscript statistics | Summarize candidate-method overlap and support at environment and BioProject levels using fixed analysis scopes. | [Figure analysis code](analysis/figures/) |

The default DKI scoring entry uses the unified R definition. The historical full-`qtrn` Lowmem baseline is not the default workflow. Scientific formulas, thresholds, seeds and their evidence levels are documented in [PARAMETERS.md](PARAMETERS.md). In particular, the fastCCLasso P-value definition follows the reviewed official implementation; see [external dependency and P-value notes](EXTERNAL_DEPENDENCIES.md).

## Repository structure

```text
pipeline/          Profiling, abundance, associations, environment networks, Zi-Pi, EPI and DKI
analysis/          External validation and manuscript analysis code
tools/             Single-project adapters, input checks and dependency checks
lib/               Shared path/configuration helpers
config/            Path-only configuration templates
tests/             Lightweight scientific tests and synthetic fixtures
frozen_inputs/     Approved analysis scope lists and aggregate reference tables
environments/      Software-version records with evidence levels
docs/              Script provenance, input availability and scope limitations
```

The reviewed collection contains **84 scientific scripts**: 35 main, 45 supplementary and 4 retained tests, with required local helpers. Their provenance is recorded in the [accepted-script manifest](docs/ACCEPTED_84_SCRIPTS.tsv). The 35 main files retain portable paths. Some supplementary or historical launchers retain original layout references and batch defaults; use the documented entry points rather than historical whole-cohort commands. See [source and dependency boundaries](docs/PRIVATE_SNAPSHOT_SCOPE.md).

## Getting started

1. Check the Python/R requirements and existing environment using [INSTALL.md](INSTALL.md).
2. Obtain and configure MetaKSSD, FastSpar or the official fastCCLasso core as needed using [EXTERNAL_DEPENDENCIES.md](EXTERNAL_DEPENDENCIES.md). Third-party code and licenses are described in [THIRD_PARTY_SOFTWARE.md](THIRD_PARTY_SOFTWARE.md).
3. Select the required inputs and isolated output directory using [PIPELINE.md](PIPELINE.md), which documents actual commands for each module. Path-only templates are provided in [config/](config/).
4. Run the lightweight checks before working with real inputs:

```bash
python tools/verify_private_checkout.py
python tests/run_small_tests.py --rscript "$(command -v Rscript)"
```

Select an existing compatible Python/R environment as described in the installation guide. The test suite uses synthetic inputs and covers Profile, EPI core/workers/2-of-3, Zi-Pi roles, unified-R DKI scoring/selection, candidate-set statistics and FastSpar postprocessing. **All nine groups passed in this checkout on 2026-10-19.** Test commands, assertions and limits are documented in [TESTING.md](TESTING.md); small input/output examples are in [tests/fixtures/](tests/fixtures/).

## Reproducibility

Reproducibility is documented at three levels: small scientific tests, independent project-level reproduction, and manuscript analysis from saved results. The nine lightweight groups check numerical and selection behavior; they do not perform full association inference, network optimization, neural-network training or a cohort rerun. Prior independent reproduction covered three human-gut projects and 111 samples, starting from existing Sketches. Software records distinguish confirmed production versions, independent reproduction versions and current local evidence.

Paper analyses and figures require the corresponding validated results, input permissions and documented project scopes. Available code includes association/network processing, keystone-candidate selection, external validation and cross-method statistics. Figure 3 uses EPI v3 and unified-R DKI with a 1,571-project main scope and 1,570/1,568-project sensitivity scopes. Its frozen-input provenance, authorization boundaries and real-data reproduction instructions are maintained in the [specialist input guide](frozen_inputs/figure3/README.md).

**Complete manuscript reproduction is not yet available from this checkout.** Figure 2 final assembly, Supplementary Table S3 literature assembly, S7/S8 generators and other unresolved final figure/table source bindings remain follow-up items. Canonical genus-abundance generation and website source also await final binding or export. MicrobiomeNet/Stage5B adoption is not represented as a default main analysis while its final manuscript scope remains unresolved. Some external validation helpers and datasets are unavailable; their gaps are listed in [DEPENDENCY_RESOLUTION.tsv](DEPENDENCY_RESOLUTION.tsv).

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the verified scope and limitations, [DATA_AVAILABILITY_GUIDE.md](DATA_AVAILABILITY_GUIDE.md) for included and separately required data, and [follow-up items](docs/FOLLOW_UP_ITEMS.md) for unresolved publication materials. Full fresh-environment installation and independent recomputation of all 1,575 projects have not been established by the lightweight tests.

## License and data use

The project-wide license awaits author confirmation. Upstream notices and third-party license obligations remain applicable; this checkout does not grant a new license for third-party software or research data. Final licensing, data-upload permissions and long-term archival deposition are tracked in the [follow-up list](docs/FOLLOW_UP_ITEMS.md).
