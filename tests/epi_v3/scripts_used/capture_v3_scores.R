#!/usr/bin/env Rscript
# Run the byte-identical production core once; export its existing score objects afterward.
core_path <- Sys.getenv("EPI_CAPTURE_CORE")
stopifnot(nzchar(core_path), file.exists(core_path))
source(core_path, local = .GlobalEnv)
stopifnot(identical(B1, select_outliers(D1)),
          identical(B2, select_outliers(D2)),
          identical(BQ, select_outliers(Q)))
objects <- list(species=species, prevalence=prevalence, D1=D1, D2=D2, Q=Q,
                D1_candidate=B1, D2_candidate=B2, Q_candidate=BQ)
saveRDS(objects, file.path(opts$outdir, paste0(opts$output, "_scores.rds")),
        compress=FALSE, version=2)
options(digits=17)
scores <- data.frame(name=species, prevalence=prevalence,
                     eligible=prevalence >= opts$prevalence_min & prevalence <= opts$prevalence_max,
                     D1=D1, D2=D2, Q=Q,
                     D1_candidate=B1, D2_candidate=B2, Q_candidate=BQ)
write.table(scores, file.path(opts$outdir, paste0(opts$output, "_scores.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE, na="NA")
# These are post-run audit statistics, not inputs to the core or a second computation.
stats_rows <- lapply(c("D1", "D2", "Q"), function(metric) {
  values <- get(metric, envir=.GlobalEnv)
  ok <- is.finite(values)
  mu <- if (any(ok)) mean(values[ok]) else NA_real_
  sd_value <- if (sum(ok) >= 2) stats::sd(values[ok]) else NA_real_
  data.frame(metric=metric, finite_n=sum(ok), background_mean=mu,
             sample_SD=sd_value, candidate_threshold=mu+2*sd_value,
             comparison="strict >", minimum_finite_for_selection=2)
})
write.table(do.call(rbind, stats_rows),
            file.path(opts$outdir, paste0(opts$output, "_thresholds.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE, na="NA")
capture.output(sessionInfo(), file=file.path(opts$outdir, paste0(opts$output, "_sessionInfo.txt")))
cat("R_version:", R.version.string, "\nvegan_version:", as.character(packageVersion("vegan")), "\n")
