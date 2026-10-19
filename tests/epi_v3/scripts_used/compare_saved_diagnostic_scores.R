#!/usr/bin/env Rscript
# Read saved RDS only. This file never sources or executes an EPI algorithm.
args <- commandArgs(trailingOnly=TRUE)
root <- args[[1]]
audit <- args[[2]]
stopifnot(length(args) == 3L)
reference <- args[[3]]
rows <- list()
for (project in c("PRJEB81106", "PRJNA1233922", "PRJNA1231247")) {
  previous <- readRDS(file.path(reference, project, "v3_core", "full_precision_scores.rds"))
  reproduced <- readRDS(file.path(root, "scores", paste0(project, "_scores.rds")))
  stopifnot(identical(previous$species, reproduced$species),
            identical(previous$prevalence, reproduced$prevalence))
  for (metric in c("D1", "D2", "Q")) {
    a <- previous[[metric]]
    b <- reproduced[[metric]]
    finite <- is.finite(a) & is.finite(b)
    old_flag <- previous[[switch(metric, D1="B1", D2="B2", Q="BQ")]]
    new_flag <- reproduced[[paste0(metric, "_candidate")]]
    mean_difference <- abs(mean(a[finite])-mean(b[finite]))
    sd_difference <- abs(sd(a[finite])-sd(b[finite]))
    cutoff_difference <- abs((mean(a[finite])+2*sd(a[finite]))-(mean(b[finite])+2*sd(b[finite])))
    rows[[length(rows)+1L]] <- data.frame(
      project=project, metric=metric, reference_kind="PREPRODUCTION_DIAGNOSTIC_NOT_PRODUCTION_SCORES",
      reference_path=file.path(reference, project, "v3_core", "full_precision_scores.rds"),
      scored_taxa=sum(finite), names_and_prevalence_identical=TRUE,
      NA_pattern_identical=identical(is.na(a),is.na(b)),
      finite_pattern_identical=identical(is.finite(a),is.finite(b)),
      scores_identical=identical(a,b),max_abs_difference=max(abs(a[finite]-b[finite])),
      background_mean_abs_difference=mean_difference,sample_SD_abs_difference=sd_difference,
      candidate_threshold_abs_difference=cutoff_difference,
      flags_identical=identical(old_flag,new_flag))
  }
}
result <- do.call(rbind,rows)
stopifnot(all(result$scores_identical),all(result$flags_identical),
          all(result$NA_pattern_identical),all(result$max_abs_difference==0))
write.table(result,file.path(audit,"SUPPLEMENTARY_DIAGNOSTIC_SCORE_COMPARISON.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
print(result[,c("project","metric","scored_taxa","max_abs_difference","flags_identical")])
