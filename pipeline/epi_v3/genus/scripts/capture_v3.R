#!/usr/bin/env Rscript
# Audit wrapper only. Execute the pinned, unmodified core before exporting its objects.
core_path <- Sys.getenv("EPI_CAPTURE_CORE")
stopifnot(nzchar(core_path), file.exists(core_path))
source(core_path, local = .GlobalEnv)
stopifnot(identical(B1, select_outliers(D1)), identical(B2, select_outliers(D2)),
          identical(BQ, select_outliers(Q)))
saveRDS(list(name = species, prevalence = prevalence, D1 = D1, D2 = D2, Q = Q,
             D1_candidate = B1, D2_candidate = B2, Q_candidate = BQ),
        file.path(opts$outdir, "audit_scores.rds"), compress = FALSE, version = 2)
options(digits = 17)
write.table(data.frame(name = species, prevalence = prevalence, D1 = D1, D2 = D2,
                       Q = Q, D1_candidate = B1, D2_candidate = B2, Q_candidate = BQ),
            file.path(opts$outdir, "audit_scores.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat("R_version:", R.version.string, "\nvegan_version:", as.character(packageVersion("vegan")), "\n")
