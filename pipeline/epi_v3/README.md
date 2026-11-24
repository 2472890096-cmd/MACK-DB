# EPI v3: one implementation for both taxonomic levels

EPI starts from the abundance matrix and exact taxon mapping, independently of network inference. Use [tools/epi_project.py](../../tools/epi_project.py) for a bounded isolated project; both species and genus now share [one R core](core/get_keystone_EPI_single_safe_qecdf_v3_fast.R) and [one Python worker](scripts/workflow.py). Rank-specific launchers/configuration select the inputs and output root; they do not implement different mathematics.

The accepted definitions remain inclusive25–75% prevalence, focal removal and renormalization, D1/D2 and the corrected unique off-diagonal sample-pair Q/ECDF with pQ0.1. Each metric uses strict score>mean+2sampleSD and the final candidate needs at least2of3. A valid zero-candidate result requires complete success/status/summary evidence.

The adapter freezes a runtime worker/core/configuration in each new output root, verifies its inputs and supports explicit resume. Such generated snapshots are audit artifacts, not duplicate source implementations. The full genus production coordinator remains an advanced entry; the default example does not launch it.

The archived PRJEB81106 complete-project record includes both ranks and valid zero candidates. The fixed-three-project historical score diagnostics were removed from the active tree, but the current nine-group scientific tests still exercise the unchanged core, workers, resume, thresholds and2-of-3. See [TESTING](../../TESTING.md) for commands and dated evidence.

Read [PIPELINE](../../PIPELINE.md), [the real matrix example](../../examples/real_bioproject/PRJEB81106/README.md) and [SCRIPT_REFERENCE](../../docs/SCRIPT_REFERENCE.md) for commands and formats.
