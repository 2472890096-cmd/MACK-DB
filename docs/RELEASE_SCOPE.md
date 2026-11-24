# Current scientific release scope

This working version follows the current manuscript: main Figures 1–4 and Table 1; Supplementary Methods S1–S8, Figures S1/S2 and Tables S1–S4. Its active inventory is [RELEASE_SCRIPT_MANIFEST.tsv](RELEASE_SCRIPT_MANIFEST.tsv); the historical 84/93/151 manifests remain immutable provenance.

- The real examples are the paired-FASTQ ERR14036375 upstream route and the complete archived 33-Run PRJEB81106 matrix workflow.
- EPI v3 has one shared R core, one shared Python worker and one species reference configuration; species/genus launchers select inputs and isolated output roots. Per-run frozen copies and configuration records are execution artifacts, not duplicate scientific implementations in the repository.
- Figure 4B subgroup calculations and their preparation/validation helpers are under `analysis/paper_network_comparison/original/matched_subgroups/`; the old `logs/` source location is retired. The different `stage4b_recurrence_recovery` task remains part of literature-support analysis.
- Six optional historical/extension scripts were removed from the active tree: Stage5A's directed complementarity builder, the two 04c logistic/recovery-extension scripts, the old Query forwarding helper and the two fixed-three-project EPI diagnostics. Their original bytes remain in Git history and the separate local preparation archive. This does not remove the null or support calculations used in Supplementary Figures S1/S2.
- EPI v3, unified-R DKI `MACKDB_DKI_R_20261012_v1`, accepted scientific parameters, synthetic expected values and real input/reference data are unchanged. The genus matrix follows the author-adopted species-label token rule.

This is a local structure revision on `paper-reproducibility`; no new commit, push, main merge or formal Release is implied. [TESTING](../TESTING.md) distinguishes recorded execution evidence from the checks performed for this revision. Large processed paper inputs, third-party source/data rights, licensing and persistent download links remain separate access decisions; see [paper reproduction](PAPER_REPRODUCTION.md).
