# Historical SRA stream: source review, not execution acceptance

The recovered `20231227-20250105_run/sketch.sh` source fragment has SHA256 **4fd7a1406964bfbf47a71344b4879b5b5d5458d7e9f0ee4c7246207be4ff403f** (306bytes). Its relevant command is:

```bash
RUN=$1
OUT=${RUN}_sketch
LOG=${RUN}.log
ERR=${RUN}.err
$KSSD/kssd dist -A -L $SHUF -o $OUT --pipecmd "fastq-dump --skip-technical --split-spot -Z" $RUN 1>$LOG 2>$ERR
```

The command above preserves the actual unquoted historical argv and redirects; KSSD/SHUF were hard-coded historical paths, deliberately omitted from this portable documentation. It is source evidence, not a recommended safe executable wrapper. The source's `RUN` is passed as the positional input and appended to pipecmd by MetaKSSD. A local SRA path is therefore the intended streaming input. An accession can cause SRA Toolkit retrieval if no local file is present; the candidate never executes that behavior automatically. `--skip-technical --split-spot -Z` sends biological reads as separate FASTQ records to stdout, not one combined nucleotide sequence per pair.

The separate `script/sketch.sh` (264bytes, SHA256 **7ef87a502a253c86fbbcd5484f104f99aa3c14a48347c4567afba088ecabcf52**) takes `${SAMPLE}.fastq` directly and does not contain the streaming flag. They are not the same file/version or an interchangeable provenance claim.

ERR14036375 logs record AWS S3/SRA download messages, `mt_shortreads2koc`, shuffleID1481063132, decomposition of `/tmp/ERR14036375`,1,367,526spots and temporary-SRA cleanup. The complete source wrapper containing both AWS acquisition and cleanup was **not located in the bounded Dell source review**. Log phrases alone do not establish its exact S3URL, executable, all argv or cleanup implementation. No guessed AWS command is supplied. The recovered source fragment itself has no AWS download or removal command.

Inspect a nonexecuting plan:

```bash
python pipeline/metakssd_profile/run_single_run.py plan-sra \
  --sra-file data/ERR14036375/ERR14036375.sra --fastq-dump fastq-dump
```

This command prints **NOT_VERIFIED_SRA_STREAM_PLAN_ONLY** and performs no network access, conversion, Sketch creation or cleanup. It is not a tested SRA execution mode and has no execute switch. No SRA file was downloaded for this integration. The tested default remains [paired FASTQ](../../examples/real_fastq/ERR14036375/README.md).

Before a future executable SRA option is accepted, recover the complete wrapper or explicitly define a new independently validated route; establish SRA input identity/checksum and toolkit version; prevent implicit retrieval; check all pipe exits; isolate local SRA and counted Sketch paths; verify read/pair semantics and full scientific outputs; and restrict cleanup to owned temporary files after verified success. Historical destructive cleanup must not be copied without these guards. This review changes no scientific formula and grants no new code/data license.
