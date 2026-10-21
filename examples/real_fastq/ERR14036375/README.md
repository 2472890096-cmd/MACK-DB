# ERR14036375: one real paired-end FASTQ to Species Profile

This is one Run from the **fixed archived set of 33 Runs from PRJEB81106** cohort, not today's 795-Run project. It contains 1,367,526 read pairs, 2,735,052 reads and 208,338,748 bases. The 2026-10-21 independent execution generated 31,236 Sketch keys with exactly the historical counts, 12 identical scientific Query records and a **byte-identical five-species Profile**. See [accepted evidence](accepted_evidence_summary.json). The original experiment and the [actual portable runner regression](portable_wrapper_regression.json) are distinct evidence records. The latter reused cached inputs in relocated code/output paths containing spaces; it did not download data.

Only small processed references, checksums and code are included. FASTQ/SRA, Sketch, MarkerDB and cohort matrices are obtained separately and must stay out of Git. No command below downloads the other 32 Runs, invokes AWS, generates genus matrices or runs a keystone model.

## 1. Obtain or reuse exactly two files

Reuse existing local files whenever available. Print the fixed URLs, bytes, MD5 and actual SHA256 without network access:

```bash
python pipeline/metakssd_profile/run_single_run.py data-info
```

If downloading this Run is authorized in your environment, these are manual retrieval commands, not an automatically invoked stage:

```bash
mkdir -p data/ERR14036375
curl --fail --location --proto '=https' --connect-timeout 20 --max-time 600 --retry 0 \
  --continue-at - --output data/ERR14036375/ERR14036375_1.fastq.gz.part \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR140/075/ERR14036375/ERR14036375_1.fastq.gz
curl --fail --location --proto '=https' --connect-timeout 20 --max-time 600 --retry 0 \
  --continue-at - --output data/ERR14036375/ERR14036375_2.fastq.gz.part \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR140/075/ERR14036375/ERR14036375_2.fastq.gz
printf '%s\n' \
  'f80a3e2e7acd28fb47c90bdd7e3783e1  data/ERR14036375/ERR14036375_1.fastq.gz.part' \
  '13cd3fe6f4058059e2f9ac7f81714b19  data/ERR14036375/ERR14036375_2.fastq.gz.part' | md5sum --check
```

**Stop if either curl or checksum command fails.** Inspect the retained partial file/log; do not rename an incomplete download. Only after both checks succeed, rename `.fastq.gz.part` to `.fastq.gz` without overwriting an existing file (`mv -n`), or point your private config directly at the already verified `.part` files. The runner independently checks exact bytes, MD5, SHA256, complete gzip streams, paired IDs and read/base totals. It never downloads or retries. The original retrieval experienced interrupted transfers; verified byte ranges were needed. No stable download speed is promised.

| Mate | Compressed bytes | SHA256 |
|---|---:|---|
| R1 |68,661,372|`a4dc2a29d01196c09337eeb2bb84a18b5d8a0a5d84716d4ae661deed75eb04f4`|
| R2 |71,325,061|`1e0fe61bb519461b9b1e1e5e6cea1a114631c26b599dba22f5bbb11fec6aa905`|

## 2. Configure existing software and reference paths

Use Linux, Python >=3.10, GNU gzip/zcat, `/usr/bin/time`, MetaKSSD **2.21**, matching L3K11 shuffle, the pinned GTDB R226 MarkerDB, and R with tidyverse/optparse. No installer is run. Official acquisition/build sources are in [EXTERNAL_DEPENDENCIES](../../../EXTERNAL_DEPENDENCIES.md).

```bash
cp -n config/single_run.example.json config/single_run.local.json
# Edit only the local JSON: FASTQ paths, executable paths, shuffle, MarkerDB and scratch root.
python tools/verify_release_candidate.py
python pipeline/metakssd_profile/run_single_run.py preflight --config config/single_run.local.json
```

Relative input/reference paths resolve from the JSON directory; executable names can use PATH. There are no Dell paths in the template. Choose a short writable scratch root such as `/tmp` with at least 3 GiB available; the runner creates its own safe directory there. Inputs and final output paths may contain spaces; internal native temporary paths are kept short and shell-safe. Reference paths must fit the native PATHLEN limit. Budget at least 4 GiB RAM for one worker; measured original Sketch peak was 579 MiB. No credentials belong in the JSON.

[reference_pins.json](reference_pins.json) records the exact accepted binary, shuffle and all three MarkerDB hashes. The default binary pin is the observed accepted Linux binary. If rebuilding 2.21 on another Linux system, review and record that build's actual SHA in the explicit local `metakssd_sha256` field. Version 2.21 and unchanged reference/core pins remain mandatory. Such a build is reported as different from the recorded binary; only an actual reference-matching run establishes output acceptance on that machine. A source commit alone does not prove compilation identity.

## 3. Run once in a new isolated directory

```bash
python pipeline/metakssd_profile/run_single_run.py run \
  --config config/single_run.local.json --output work/ERR14036375_first_run
python pipeline/metakssd_profile/run_single_run.py verify --output work/ERR14036375_first_run
```

The runner implements exactly this sequence:

1. Check compressed files/references/core/executable pins, MetaKSSD version, R dependencies and disk space; completely read both gzip streams and match every paired recordID. No filtering/trimming/reordering.
2. `gzip -cd R1 R2` creates **one collection of separate FASTQ records**, not concatenated nucleotide sequences. Both mates contribute counts to one Run.
3. Official `metakssd dist -L L3K11.shuf -A -p 1 -m 2 -o NEW_SKETCH RUN.fastq`;22-mer, shuffleID1481063132, nominal1/4096. One sample, count-preserving Sketch.
4. Official `metakssd composite -r PINNED_MARKERDB -q NEW_SKETCH` generates the headerless **seven-column** Query. The first column is the local input path.
5. Unchanged `pipeline/metakssd_profile/profile.R -l SINGLE_RUN.txt -i QUERY_DIR -s 12 -o PROFILE_DIR`: X3>=12, existing X5 Poisson lookup/subtract3, percentages rounded4decimals. **12 is shared-k-mer support, not k-mer length**; upstream generic Perl default 18 is not used.
6. Compare complete sorted key/count content, every scientific Query field by taxon, and the Profile SHA against the small accepted references. Any difference fails; no parameter tuning/retry.

The external stages have 600 s limits and one-thread environment controls. Failure/interruption retains logs and owned scratch, never claims success and never retries. Existing output roots are rejected. Success commits complete output directories, records final status and hashes, and removes only the ownership-checked scratch directory created by this run. `--keep-scratch` explicitly retains it. No original input/reference or historical result is removed.

## Outputs and acceptance

- `sketch/`: cofiles.stat, combco.0, combco.0.a and combco.index.0.
- `query/ERR14036375_all.txt`: seven headerless fields; Profile consumes X2taxon, X3shared count and X5depth proxy.
- `profile/ERR14036375_profile.txt`: headerless `ID_species_label<TAB>percentage`, identical to [expected_profile.tsv](expected_profile.tsv).
- `RUN_STATUS.json`, `OUTPUT_SHA256.tsv`, `logs/commands.jsonl`, individual stdout/stderr/time-v/status logs. Output files alone are not success evidence.

Sketch serialization order, embedded local path and C padding can differ. Sorted key/count content must match exactly; the original experiment found 15 different key positions but identical key/count mappings. Query input path is exempt; all scientific columns must match. The final Profile itself must be byte-identical. Its rounded percentage sum 100.0001 equals the accepted reference; do not renormalize it again.

Original measured Sketch/Query/Profile wall times were 4.62/3.42/1.82 s; pair checking+decompression+science total 27.1 s, excluding network and audit preparation. These are observations on one server, not performance guarantees or evidence for all 33 Runs.

## Separate scientific boundaries

The [33Run matrix example](../../real_bioproject/PRJEB81106/README.md) starts at archived 53×33 species / 38×33 genus matrices. It is separate from this fresh single-Run upstream acceptance. Genus follows the author-adopted original species-name token aggregation; this single-Run species example does not recompute genus. DKI saved-prediction scoring replay is not cNODE2 training. No complete new neural-training or full FASTQ-to-keystone acceptance is claimed.

The [historical SRA-stream review](../../../pipeline/metakssd_profile/HISTORICAL_SRA_STREAM.md) is an optional **NOT_VERIFIED** plan only. Its original AWS download/cleanup wrapper must not be inferred from log phrases or silently substituted for the tested paired-FASTQ route.
