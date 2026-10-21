# Scientific tests and source integrity

Run from the repository root in an already available environment:

```bash
python tools/verify_release_candidate.py
python tests/run_small_tests.py --rscript "$(command -v Rscript)"
```

On 20 October 2026, all nine groups actually passed in this candidate: PROFILE, EPI_CORE, EPI_WORKERS, EPI_2OF3, ZIPI_ROLES, DKI_R_SCORING, DKI_SELECTION, FIGURE3_SETS and FASTSPAR_POSTPROCESS. Current timings are in [tests/expected/PHASE1_TEST_RESULTS_20261020.tsv](tests/expected/PHASE1_TEST_RESULTS_20261020.tsv). The historical SMALL_TEST_RESULTS.tsv remains unchanged. Full command/output logs are retained locally outside the candidate publication files and under ignored `tests/results/`.

Only the concordance module import path changed in the scientific test driver. Assertions, numerical expected values, independent oracles and R unit-test sources were not changed. The legacy FIGURE3_SETS label now exercises `analysis/keystone_concordance/concordance_release.py`. The complete concordance program itself is byte-identical to the accepted source.

Coverage: Profile X3/sharekmer12 and X5 Poisson/precision boundaries; EPI prevalence, discrete ECDF ties/fallback, D1/D2/Q, sample-SD thresholds, valid zero candidates, resume and 2-of-3; fixed-module Zi/Pi formulas and roles; unified-R DKI q_null removal/normalization, Bray–Curtis, median and selection; Broad/Strict mutually exclusive combinations and unique taxa; FastSpar postprocessing and labelled-matrix alignment. The old full-qtrn Lowmem baseline fails the R-definition counterexample.

Not covered by those nine synthetic groups: real FASTQ retrieval/sketch generation, complete FastSpar 10000-bootstrap inference, fastCCLasso fitting, signed module optimization/plotting, DKI neural training/prediction, complete-project S9 analysis or any1575-project rerun. The accepted real single-Run path is documented separately below. Static syntax/help checks are recorded separately from scientific execution. Missing dependencies must be BLOCKED; assertions must never be weakened to obtain PASS.

The new version guard verifies immutable original84/93/151 provenance, all current script hashes, reviewed original93 dispositions, file hashes, EPI aliases, fixed scopes, approved six frozen inputs, the 64 held-input boundary and synthetic fixture hashes. It does not silently disable the old count check: the original84 snapshot remains independently pinned while the current version has its own manifest. Generated test outputs are ignored.


## Stage2 real matrix acceptance

The PRJEB81106 matrix example was additionally exercised in isolated local directories at both ranks: formal10000FastSparbootstrap/permutations,100fastCCLassobootstraps,smallZiPi,EPI,andDKIinput/scoringreplay. This is separate from the nine synthetic groups, whose assertions/oracles were unchanged. The new input contract was tested with minimal branch-specific JSON and rejects wrongSHA or genus computation without the explicit historical-input flag. Initial speciesZiPi validation used the wrongcolumncase; actualz/p columns were then compared read-only with no scientific rerun or changedexpectedvalues. Real-data sharing remains separately controlled; genus follows the author-adopted name-token grouping. DownstreamPASS is not a claim of whole-cohort end-to-end reconstruction or unrestricted publicfixture access. See examples/real_bioproject/PRJEB81106/README.md and the separate localStage2reports.

## Real one-Run upstream acceptance and portable wrapper safety

The2026-10-21 ERR14036375 experiment actually downloaded and verified only its two approved mates, generated a new official2.21 count Sketch, queried the pinnedR226MarkerDB and ran the unchanged RProfile-s12. All31,236key/count mappings,12scientific Query rows and all5Profile rows matched; Profile bytes are identical. [Small evidence summary](examples/real_fastq/ERR14036375/accepted_evidence_summary.json) contains the prior measurements and scope. This is neither33Run fresh upstream acceptance nor complete DKI training.

```bash
python tests/test_single_run_workflow.py
```

These additional small regressions test paired-record validation, corruption/mismatches, relative paths, rejected parameter overrides, nonzero exits, bounded child-session timeout, safe owned cleanup and nonexecuting SRA plans. They perform no download or neural training and are **safety/IO tests**, not proof of real scientific output. A relocated-copy run using the existing verified real FASTQs is the separate actual CLI/science acceptance. The original nine-group scientific test driver, its expected values and core algorithms are unchanged.

Actual2026-10-21 integration results are in [UPSTREAM_INTEGRATION_20261021.tsv](tests/expected/UPSTREAM_INTEGRATION_20261021.tsv). Original nine groups all passed again,12new safety cases passed, and the [relocated cached-input CLI](examples/real_fastq/ERR14036375/portable_wrapper_regression.json) matched the real reference. This is same-server path relocation, not a fresh second-machine installation test.

## Portable complete-project matrix workflow

```bash
python tests/test_portable_project_workflow.py
python examples/real_bioproject/PRJEB81106/run_workflow.py --steps all --preflight
python examples/real_bioproject/PRJEB81106/run_workflow.py --steps all --output work/PRJEB81106
```

Seven new I/O tests check independent branch planning, fresh-training dependency, bundle schema/hashes, altered-artifact refusal and nonfinite prediction rejection. They are separate from the unchanged nine scientific groups and twelve upstream safety cases. The real33Run example also actually ran all16stages at both ranks with formal parameters, including fresh training/prediction (15.05/11.45s) and unified-R selection (3candidates each). The [execution table](tests/expected/PORTABLE_MATRIX_EXECUTION_20261021.tsv) records real stage timings; complete local logs/results remain outside publication files. A subsequent code relocation and resume check tests path portability without another training run. These checks do not establish a fresh environment install on every server.
