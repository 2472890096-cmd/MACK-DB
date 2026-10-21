# Figure 3 frozen inputs and reproduction

This guide describes the Figure 3 data boundary within the broader [MACK-DB research workflow](../../README.md). The final statistics use species/genus EPI v3 and unified-R DKI scoring **MACKDB_DKI_R_20261012_v1**. The main analysis covers 1,571 BioProjects; sensitivity analyses cover 1,570 and 1,568. Project eligibility is fixed by these lists, including legitimate zero-candidate projects.

## Included files and authorization

Six of the 70 reviewed v0.2 frozen input/reference files are included, copied byte-identically:

- [scopes/](scopes/) contains the three accepted 1,571/1,570/1,568 project lists.
- [expected/](expected/) contains the final unified-R main 22-value table, sensitivity results and old/new aggregate differences.

The [70-file availability manifest](../../docs/FIGURE3_INPUT_AVAILABILITY.tsv) records each file's SHA256 and include/hold decision. The other 64 files are held locally: 62 row-level research/provenance files await separate data-upload authorization, and two historical aggregate assertion tables are unnecessary without the full inputs. Held materials include candidate snapshots, the 61-row R-scoring replacement, taxon-support records and environment mappings/lists. Scientific acceptance or private repository visibility does not establish permission to upload them.

This subset is **not a self-contained real-data reproduction package**. Do not substitute a historical Lowmem candidate table for the accepted final science or fabricate missing records. Data availability and reuse limits are also documented in [DATA_AVAILABILITY_GUIDE.md](../../DATA_AVAILABILITY_GUIDE.md).

## Final statistics and input contract

The unchanged [figure3_release.py](../../analysis/keystone_concordance/concordance_release.py) checks the historical baseline and replaces the Species DKI candidates for PRJEB80533, PRJNA1345363, PRJNA1357645 and PRJEB72274 with all 61 accepted unified-R candidates. Other method inputs remain unchanged. The baseline must be the reviewed pre-replacement snapshot; passing an already-replaced snapshot as that baseline violates the program's assertions.

The full command and required directory layout are documented in the [Figure 3 section of PIPELINE.md](../../PIPELINE.md#figure3--final-code-retained-real-inputs-held). Real-data execution requires all of the following separately authorized inputs:

- `base_analysis/inputs/key_{species,genus}{,_epi,_dki}.tsv`: six frozen candidate tables with the reviewed name, BioProject and environment representation.
- `base_analysis/FIGURE3_ALL_FOUR_SCOPES_22_VALUES.tsv`: historical baseline assertions.
- `r_replacement/FOUR_PROJECT_R_CANDIDATES_FINAL.tsv`: the accepted 61 unified-R candidates.
- `impact_audit/FIGURE3_POTENTIAL_IMPACT.tsv`: accepted unified-R impact assertions.
- The three included project lists in `scopes/`.

After approval, supply the complete validated package under ignored `data/`, verify source hashes against the availability records and use a separate output directory. The statistics program reads saved candidates only; it does not rerun EPI, DKI, network inference or database queries. Rendering is a separate step.

## Reproduction evidence and limits

The preceding v0.2 full-input package reproduced the final 22 main values, 44 sensitivity values and complete Strict All3/human-gut support records. Main Strict All3 contains 25 species and 19 genera. That is prior full-data evidence; the row-level inputs needed to repeat it are not bundled in this checkout.

The synthetic Figure 3 set test passed here as one of the nine lightweight groups. It verifies mutually exclusive combinations, Broad/Strict support, inclusion of All3 in at-least-two support and taxon deduplication without real project data. See [TESTING.md](../../TESTING.md) and the overall [reproducibility scope](../../REPRODUCIBILITY.md). A synthetic pass is not a new full-data Figure 3 reproduction.
