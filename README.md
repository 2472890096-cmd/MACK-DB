# MACK-DB

MACK-DB integrates microbial associations and candidate keystone taxa from **264,767 shotgun metagenomic sequencing runs, 1,575 BioProjects and 45 environments**.

The workflow covers MetaKSSD/GTDB R226 profiling, abundance matrices, FastSpar and fastCCLasso, Zi-Pi, EPI v3, unified-R DKI, and literature validation. **Zi-Pi uses networks; EPI and DKI start independently from abundance matrices.**

## Run a real example

Follow the [Linux installation guide](environments/portable-matrix/README.md), then configure your local tool paths. The repository includes PRJEB81106's **33 archived Runs**, species **53×33** and genus **38×33** abundance matrices, mappings and small reference outputs.

From the repository root:

```bash
cp -n config/portable_matrix.example.json config/portable_matrix.local.json
# Edit the local configuration after installing the documented dependencies.
example=examples/real_bioproject/PRJEB81106/run_workflow.py
config=config/portable_matrix.local.json

"$MACKDB_PYTHON" "$example" --config "$config" --steps all --preflight
"$MACKDB_PYTHON" "$example" --config "$config" --steps all \
  --levels species genus --output work/PRJEB81106_complete --threads 4
```

Check `WORKFLOW_STATUS.json` and each stage's `SUCCESS.json`. Prior full runs completed **16 stages**, including fresh DKI training/prediction. FastSpar references contain **25 species and 3 genus associations**; EPI zero candidates are valid. See the [example guide](examples/real_bioproject/PRJEB81106/README.md) for parameters, expected outputs and numerical comparison limits.

For the upstream FASTQ→Sketch→Query→Profile route, use the separately validated [ERR14036375 example](examples/real_fastq/ERR14036375/README.md).

## Reproduce the paper

[Paper reproduction](docs/PAPER_REPRODUCTION.md) provides inputs, commands, outputs and checks for **Figures 2–4, Supplementary Figures S1/S2 and related tables**. The main candidate comparison uses **1,571 BioProjects**, with Strict All3 support for **25 species and 19 genus groups**.

Real-input statistical checks passed on the server. The six required Figure 2/4 sources are included and used by default. Paper-wide processed-data packages still need separate access; the checkout alone does not yet supply those inputs. See [data access](docs/paper/DATA_PACKAGES.md) and [validation scope](REPRODUCIBILITY.md). Large research data stay outside GitHub.

## Repository guide

| Directory | Contents |
|---|---|
| `pipeline/` | Profiling, abundance, associations, keystone methods and environment integration |
| `analysis/` | Paper statistics and external validation |
| `examples/` | The two real-data starting points |
| `environments/`, `config/` | Installation requirements and nonsecret configuration templates |
| `lib/`, `tools/` | Required shared code and input/version checks |
| `tests/` | Small scientific tests, synthetic inputs and expected results |
| `docs/`, `frozen_inputs/` | Detailed instructions, input manifests, fixed scopes and small references |

[Workflow](PIPELINE.md) · [Parameters](PARAMETERS.md) · [Tests](TESTING.md) · [Dependencies and third-party terms](THIRD_PARTY_SOFTWARE.md) · [File organization](docs/GITHUB_CONTENTS.md)
