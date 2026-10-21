> Historical planning document. For the current runnable complete-project matrix example, use [PRJEB81106](../real_bioproject/PRJEB81106/README.md); for the accepted single-Run upstream route, use [ERR14036375](../real_fastq/ERR14036375/README.md).

# Matrix-start PRJEB81106 plan

Not executed in phase1. Prerequisites: approved53x33species abundance/mapping, ordered33Run IDs, existingR environment, FastSpar/GNUparallel, taxonomy and full ZiPi dependencies. Existing38x33genus matrices may form a separate downstream branch; they do not validate canonical genus upstream generation.

1. Verify INPUT_REQUIREMENTS hashes, taxon/sample axes, finite/nonnegative values and identities. Reject mismatch; no guessed transposition or reordering.
2. Follow PIPELINE: count_sparcc.R -n1000000 → full run_fastspar.sh10000boot protocol → network_matrix_io → sparcc_corr_pvalue/corr_inter final|r|≥.3,P≤1e-4. A4–8CPU/8–16GiB budget and5–30min planning window are unbenchmarked; the recorded29-thread42.511s FastSpar run does not establish this budget.
3. Run the retained full signed ZiPi entry and candidate extraction, then tools/epi_project.py for EPIv3 with its complete success/summary/hash evidence. Historical reference: species25FastSpar edges(6positive/19negative), ZiPi6nodes/4edges/3modules and0candidates; EPI3eligible,0D1/D2/Q/2of3candidates. Zero keys without successful scientific evidence are not completion.
4. Default lightweight DKI option uses separately approved savedPtrain/Ptest/species_id/sample_id/qtrn/qtst → acceptedRscorer → name mapping → TopK/background selection. It is prediction replay and Rscoring validation, not training or prediction reproduction. Do not mix archived predictions with a numerically different newly produced matrix. Historical3DKI candidates must be checked from archived scores/candidate records, not chosen to satisfy an expected result.
5. A separately authorized full-model branch may use the existing093/089 entries with explicit1000updates/80patience,seed20260629,RK4.05. No such branch ran in phase1; matchingseed alone does not guarantee identical optimizer trajectories.

One project cannot reproduce an environmentf=2 network or the1571-project concordance. Full optimizer, external inference and real candidate agreement require a later actual example run; nine small tests do not establish them.
