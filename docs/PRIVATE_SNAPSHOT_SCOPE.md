# Fixed source and historical entry boundaries

docs/ACCEPTED_84_SCRIPTS.tsv records all84 reviewed files,35 KEEP_MAIN,45 KEEP_SUPPLEMENTARY and4 KEEP_TEST, using exact v0.2 hashes. No scientific source or parameter was edited for this commit. The two EPI core paths and two worker compatibility paths remain present and byte-identical.

The local review found15 supplementary/historical accepted scripts containing source-layout references, old shebangs or historical batch defaults. None is a credential, IP address or connection configuration; none belongs to the35 portable main entries. They are retained under the author's explicit instruction to preserve all84 accepted scripts. Never launch historical whole-cohort scripts by default. Use the explicit single-project/scoring commands in PIPELINE.md. Retaining these strings is a disclosed provenance exception, not a claim that every supplementary script is path-independent.

The34 registered in-repository dependency edges and the small-test helpers are included. Some optional external taxonomy/source-exclusion helpers and real datasets were absent already in v0.2 and remain unavailable; DEPENDENCY_RESOLUTION.tsv identifies them. They are not invented or fetched from another source. Files named build_release.py or report generators inside logs/ are reviewed source code, not copied runtime logs; those exact selected source files are retained while server log files are excluded.

Three nested pipeline/epi_v3/*/audit/*.json files are retained because they are accepted nonsecret parameter/path templates used by actual code. The broad internal root audit/ tree is not copied. No root LICENSE is imposed, but docs/third_party/DKI_UPSTREAM_MIT_LICENSE.txt is preserved for upstream-derived code.
