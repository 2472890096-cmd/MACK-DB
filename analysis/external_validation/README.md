# External literature validation

Use the current [paper reproduction guide](../../docs/PAPER_REPRODUCTION.md) for manuscript numbering, actual inputs, execution commands and acceptance. MainFigures1–4/Table1, MethodsS1–S8, FiguresS1/S2 and TablesS1–S4 are the current scope. Older filenames retain their historical labels.

## Scientific entries

- `literature_networks/scripts/01_parse_external_networks.py`, `02_harmonize_external_taxonomy.py`, stage2b and stage2c parse and map external names to GTDB R226. Their original publication tables and some mapping helpers remain access/rights dependencies.
- stage3c and the shared Chen/Wallen protocol generate the two1,000-replicate nulls for11 association sets. [The stored-result replay](../paper_support/README.md) now reuses their exact summary function to reproduce Table1/TableS2/FigureS1 without generating new graphs. Null seeds are unique within each historical run family; the two families have different master seeds.
- stage4a–d implement publication membership, source exclusion, recurrence and support statistics. `04d_literature_recurrence_support_rank.py` and its independent recalculator underpin MethodS7/FigureS2/TableS4. Ranks use all653,071 original pre-exclusion pairs, while counts and support-fraction denominators exclude the relevant sources. The current workbook and MethodS7 now agree on this distinction.
- Stage5B provides actual subgroup preparation/comparison and degree-preserving rewiring functions. [Figure4 adapters](../paper_network_comparison/README.md) reuse the original B/C/D functions. B includes Guo Common runs pooled **before** FastSpar. C and D use complete-project and common-species scopes with explicit threshold/denominator differences.
- Figure4A now uses the stated six-study order and fixed nine-source exclusion. Its dedicated small set check compares every step with the historical order22 archive; it neither enumerates720 orders nor averages them.

## Running and access

Current isolated checks have recalculated the saved null summaries, empirical P values and support/rank statistics. They do not validate fresh random graph generation, new taxonomy resolution, network inference or model training. Final artwork is a separate operation with external scripts/fonts and documented dependencies.

Some original sources remain held after the author's12November2026 confirmation. See [current source access](../../docs/paper/SOURCE_ACCESS.tsv) and the historical [14-source record](../../docs/PENDING_FINAL_STATISTICS_CODE.tsv). No held file is redistributed or assigned a new license here. The processed data packages also have no published download identifier yet.

Many historical entries have no harmless `--help`; executing them can write outputs or launch thousands of simulations. Use the new bounded replay interfaces for current summary checking. To regenerate simulations, follow [the exact historical command/dependency boundary](../paper_support/README.md#rerunning-simulations), with separate output configuration, inputs and resource approval. Source exclusion requires Python3.13; the matrix example's Python version does not cover that whole historical path.

See [PIPELINE](../../PIPELINE.md), [data packages](../../docs/paper/DATA_PACKAGES.md) and [remaining release work](../../docs/FOLLOW_UP_ITEMS.md).
