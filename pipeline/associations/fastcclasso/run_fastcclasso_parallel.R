#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
})

DEFAULT_ABUNDANCE_DIR <- Sys.getenv("MACKDB_SPECIES_INPUT", "data/species") |> file.path("abundance")
DEFAULT_SPECIES_DIR <- Sys.getenv("MACKDB_SPECIES_INPUT", "data/species") |> file.path("species")
DEFAULT_FAST_CORE <- Sys.getenv("FASTCCLASSO_CORE", "external/fastCCLasso/RCode/fastCCLasso.R")

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "",
              help = "Project ID list. If empty, all abundance projects are used."),
  make_option(c("--project"), type = "character", default = "",
              help = "Run one project only."),
  make_option(c("-a", "--abundance_dir"), type = "character",
              default = DEFAULT_ABUNDANCE_DIR),
  make_option(c("-s", "--species_dir"), type = "character",
              default = DEFAULT_SPECIES_DIR),
  make_option(c("-o", "--output"), type = "character", default = "./"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.0001),
  make_option(c("-c", "--cor_cutoff"), type = "numeric", default = 0.1,
              help = "Keep pairs with abs(corr) >= cutoff."),
  make_option(c("-n", "--nboot"), type = "integer", default = 100),
  make_option(c("-k", "--kmax"), type = "integer", default = 20,
              help = "Maximum iterations for golden-section lambda search."),
  make_option(c("--workers"), type = "integer", default = 2,
              help = "Number of projects run in parallel."),
  make_option(c("--seed"), type = "integer", default = 20260804,
              help = "Base seed; each project receives a deterministic seed."),
  make_option(c("--fast_core"), type = "character",
              default = DEFAULT_FAST_CORE,
              help = "Path to official fastCCLasso.R."),
  make_option(c("--run_label"), type = "character", default = "all",
              help = "Label used in group-level summary filenames."),
  make_option(c("--skip_existing"), type = "character", default = "T"),
  make_option(c("--write_full_matrices"), type = "character", default = "T"),
  make_option(c("--write_pair_csv"), type = "character", default = "T")
)

opt <- parse_args(
  OptionParser(option_list = option_list),
  convert_hyphens_to_underscores = TRUE
)

truthy <- function(x) {
  toupper(as.character(x)) %in% c("T", "TRUE", "1", "Y", "YES")
}

if (opt$workers < 1L) stop("--workers must be >= 1")
if (opt$nboot < 1L) stop("--nboot must be >= 1")
if (opt$kmax < 1L) stop("--kmax must be >= 1")
if (!is.finite(opt$pvalue) || opt$pvalue < 0 || opt$pvalue > 1) {
  stop("--pvalue must be between 0 and 1")
}
if (!is.finite(opt$cor_cutoff) || opt$cor_cutoff < 0 || opt$cor_cutoff > 1) {
  stop("--cor-cutoff must be between 0 and 1")
}

opt$output <- normalizePath(opt$output, mustWork = FALSE)
opt$abundance_dir <- normalizePath(opt$abundance_dir, mustWork = TRUE)
opt$species_dir <- normalizePath(opt$species_dir, mustWork = TRUE)
opt$fast_core <- normalizePath(opt$fast_core, mustWork = TRUE)
opt$run_label <- gsub("[^A-Za-z0-9_.-]+", "_", opt$run_label)

source(opt$fast_core)
if (!exists("fastCCLasso", mode = "function")) {
  stop("Official core does not define fastCCLasso(): ", opt$fast_core)
}

matrix_path <- function(project) {
  file.path(opt$abundance_dir, paste0(project, "_abundance.matrix"))
}

species_path <- function(project) {
  file.path(opt$species_dir, paste0(project, "_species.matrix"))
}

list_projects <- function() {
  if (opt$project != "") return(opt$project)

  if (opt$file != "") {
    x <- trimws(readLines(opt$file, warn = FALSE))
    return(unique(x[x != "" & !grepl("^#", x)]))
  }

  files <- list.files(opt$abundance_dir, pattern = "_abundance\\.matrix$")
  sort(sub("_abundance\\.matrix$", "", files))
}

ensure_dirs <- function() {
  dirs <- c(
    file.path(opt$output, "fastCCLasso"),
    file.path(opt$output, "matrix", "cor"),
    file.path(opt$output, "matrix", "pv"),
    file.path(opt$output, "cor-pv_tsvformat"),
    file.path(opt$output, "cor-pv_csvformat"),
    file.path(opt$output, "proj_links"),
    file.path(opt$output, "project_logs"),
    file.path(opt$output, "project_status"),
    file.path(opt$output, "group_summaries")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

target_files <- function(project) {
  list(
    cor = file.path(opt$output, "matrix", "cor",
                    paste0(project, "_fastCCLasso_freq_cor.tsv")),
    pv = file.path(opt$output, "matrix", "pv",
                   paste0(project, "_fastCCLasso_freq_pv.tsv")),
    pairs_tsv = file.path(opt$output, "cor-pv_tsvformat",
                          paste0(project, "_fastCCLasso_freq_sp_pairs.tsv")),
    pairs_csv = file.path(opt$output, "cor-pv_csvformat",
                          paste0(project, "_fastCCLasso_freq_sp_pairs.csv")),
    links = file.path(opt$output, "proj_links",
                      paste0(project, "_links.txt")),
    status = file.path(opt$output, "project_status",
                       paste0(project, ".tsv")),
    log = file.path(opt$output, "project_logs",
                    paste0(project, ".log"))
  )
}

required_outputs <- function(files) {
  out <- c(files$pairs_tsv, files$links)
  if (truthy(opt$write_pair_csv)) out <- c(out, files$pairs_csv)
  if (truthy(opt$write_full_matrices)) out <- c(out, files$cor, files$pv)
  out
}

read_matrix <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(
      path, sep = "\t", header = TRUE, check.names = FALSE,
      data.table = TRUE, showProgress = FALSE
    )
    if (ncol(dt) < 2L) stop("Abundance file has fewer than 2 columns: ", path)
    ids <- dt[[1L]]
    dt[[1L]] <- NULL
    mat <- as.matrix(dt)
    rownames(mat) <- ids
  } else {
    mat <- as.matrix(read.table(
      path, header = TRUE, row.names = 1, sep = "\t",
      check.names = FALSE, quote = "", comment.char = ""
    ))
  }

  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) stop("Abundance matrix contains NA/NaN/Inf: ", path)
  if (any(mat < 0)) stop("Abundance matrix contains negative values: ", path)

  # Keep preprocessing identical to the previous CCLasso workflow:
  # taxa x samples -> samples x taxa, then add a tiny pseudocount.
  t(mat) + 1e-9
}

read_species <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(
      path, sep = "\t", header = FALSE,
      data.table = TRUE, showProgress = FALSE
    )
    if (ncol(dt) < 2L) stop("Species file has fewer than 2 columns: ", path)
    return(as.character(dt[[2L]]))
  }

  df <- read.delim(
    path, header = FALSE, sep = "\t",
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
  if (ncol(df) < 2L) stop("Species file has fewer than 2 columns: ", path)
  as.character(df[[2L]])
}

atomic_rename <- function(tmp, path) {
  if (file.rename(tmp, path)) return(invisible(TRUE))
  if (!file.copy(tmp, path, overwrite = TRUE)) {
    stop("Failed to move temporary output to: ", path)
  }
  unlink(tmp)
  invisible(TRUE)
}

atomic_write_table <- function(df, path, sep = "\t", col_names = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1L))
  on.exit(unlink(tmp), add = TRUE)

  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(
      df, tmp, sep = sep, col.names = col_names,
      quote = FALSE, na = "NA"
    )
  } else {
    write.table(
      df, tmp, sep = sep, col.names = col_names,
      row.names = FALSE, quote = FALSE, eol = "\n", na = "NA"
    )
  }
  atomic_rename(tmp, path)
}

atomic_write_empty <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1L))
  on.exit(unlink(tmp), add = TRUE)
  file.create(tmp)
  atomic_rename(tmp, path)
}

write_matrix_compatible <- function(mat, out_file) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  out <- as.data.frame(mat, check.names = FALSE)
  out <- data.frame(`#OUT_id` = 0:(nrow(mat) - 1L), out, check.names = FALSE)
  colnames(out) <- c("#OUT_id", as.character(0:(ncol(mat) - 1L)))
  atomic_write_table(out, out_file, sep = "\t", col_names = TRUE)
}

project_seed <- function(project) {
  bytes <- utf8ToInt(project)
  offset <- sum(bytes * seq_along(bytes))
  as.integer((opt$seed + offset) %% (.Machine$integer.max - 1L) + 1L)
}

write_status <- function(result, path) {
  atomic_write_table(
    as.data.frame(result, stringsAsFactors = FALSE),
    path, sep = "\t", col_names = TRUE
  )
}

fold_diagnostics <- function(x, k_cv = 3L, tol = 1e-15) {
  xx2 <- log(x) - rowMeans(log(x))
  n <- nrow(xx2)
  n_b <- floor(n / k_cv)

  full_var <- apply(xx2, 2L, stats::var)
  max_test_zero <- 0L
  max_train_zero <- 0L

  for (k in seq_len(k_cv)) {
    itest <- (n_b * (k - 1L) + 1L):(n_b * k)
    test_var <- apply(xx2[itest, , drop = FALSE], 2L, stats::var)
    train_var <- apply(xx2[-itest, , drop = FALSE], 2L, stats::var)
    max_test_zero <- max(max_test_zero, sum(!is.finite(test_var) | abs(test_var) <= tol))
    max_train_zero <- max(max_train_zero, sum(!is.finite(train_var) | abs(train_var) <= tol))
  }

  list(
    full_zero_variance_taxa = sum(!is.finite(full_var) | abs(full_var) <= tol),
    max_cv_test_zero_variance_taxa = max_test_zero,
    max_cv_train_zero_variance_taxa = max_train_zero
  )
}

run_project <- function(project) {
  files <- target_files(project)

  if (truthy(opt$skip_existing) && all(file.exists(required_outputs(files)))) {
    return(list(
      project = project, status = "skip", error = "",
      elapsed_seconds = 0, n_samples = NA_integer_, n_taxa = NA_integer_,
      selected_edges = NA_integer_, finite_pair_fraction = NA_real_,
      lambda_best = NA_real_, full_zero_variance_taxa = NA_integer_,
      max_cv_test_zero_variance_taxa = NA_integer_,
      max_cv_train_zero_variance_taxa = NA_integer_
    ))
  }

  log_con <- file(files$log, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
  }, add = TRUE)

  message("[start] ", project, " ", format(Sys.time(), "%F %T"))
  t0 <- proc.time()[["elapsed"]]

  mat_file <- matrix_path(project)
  spec_file <- species_path(project)

  if (!file.exists(mat_file)) stop("Missing matrix: ", mat_file)
  if (!file.exists(spec_file)) stop("Missing species: ", spec_file)

  x <- read_matrix(mat_file)
  species_names <- read_species(spec_file)

  if (nrow(x) < 6L) {
    stop("Fewer than 6 samples; 3-fold CV would contain folds with <2 samples: ", project)
  }
  if (ncol(x) < 2L) stop("Fewer than 2 taxa: ", project)
  if (ncol(x) != length(species_names)) {
    stop("Species count does not match abundance matrix columns: ", project)
  }
  if (anyDuplicated(species_names)) {
    message("[warning] duplicated species names: ", sum(duplicated(species_names)))
  }

  diag_info <- fold_diagnostics(x, k_cv = 3L)
  message(
    "[diagnostic] full_zero_variance_taxa=", diag_info$full_zero_variance_taxa,
    " max_cv_test_zero_variance_taxa=", diag_info$max_cv_test_zero_variance_taxa,
    " max_cv_train_zero_variance_taxa=", diag_info$max_cv_train_zero_variance_taxa
  )

  set.seed(project_seed(project))

  result <- fastCCLasso(
    xx = x,
    isCnt = FALSE,
    k_cv = 3,
    lam_min_ratio = 1e-4,
    k_max = opt$kmax,
    n_boot = opt$nboot
  )

  cor_mat <- result$rho
  pv_mat <- result$p_vals

  if (!is.matrix(cor_mat) || !is.matrix(pv_mat)) {
    stop("fastCCLasso did not return matrix outputs: ", project)
  }
  if (!identical(dim(cor_mat), c(ncol(x), ncol(x)))) {
    stop("Correlation matrix dimension mismatch: ", project)
  }
  if (!identical(dim(pv_mat), dim(cor_mat))) {
    stop("P-value matrix dimension mismatch: ", project)
  }

  if (truthy(opt$write_full_matrices)) {
    write_matrix_compatible(cor_mat, files$cor)
    write_matrix_compatible(pv_mat, files$pv)
  }

  idx <- which(upper.tri(cor_mat, diag = FALSE), arr.ind = TRUE)

  pair_df <- data.frame(
    from = species_names[idx[, 1L]],
    to = species_names[idx[, 2L]],
    corr = cor_mat[idx],
    pvalue = pv_mat[idx],
    stringsAsFactors = FALSE
  )

  finite_pairs <- is.finite(pair_df$corr) & is.finite(pair_df$pvalue)
  finite_pair_fraction <- mean(finite_pairs)

  pair_df <- pair_df[
    finite_pairs &
      pair_df$pvalue <= opt$pvalue &
      abs(pair_df$corr) >= opt$cor_cutoff,
    , drop = FALSE
  ]

  atomic_write_table(pair_df, files$pairs_tsv, sep = "\t", col_names = TRUE)
  if (truthy(opt$write_pair_csv)) {
    atomic_write_table(pair_df, files$pairs_csv, sep = ",", col_names = TRUE)
  }

  if (nrow(pair_df) == 0L) {
    atomic_write_empty(files$links)
  } else {
    links <- data.frame(
      from = pmin(pair_df$from, pair_df$to),
      to = pmax(pair_df$from, pair_df$to),
      stringsAsFactors = FALSE
    )
    atomic_write_table(links, files$links, sep = "\t", col_names = FALSE)
  }

  elapsed <- proc.time()[["elapsed"]] - t0

  out <- list(
    project = project,
    status = "pass",
    error = "",
    elapsed_seconds = elapsed,
    n_samples = nrow(x),
    n_taxa = ncol(x),
    selected_edges = nrow(pair_df),
    finite_pair_fraction = finite_pair_fraction,
    lambda_best = result$lambda_best,
    full_zero_variance_taxa = diag_info$full_zero_variance_taxa,
    max_cv_test_zero_variance_taxa = diag_info$max_cv_test_zero_variance_taxa,
    max_cv_train_zero_variance_taxa = diag_info$max_cv_train_zero_variance_taxa
  )

  write_status(out, files$status)
  message(
    "[finish] ", project,
    " elapsed=", round(elapsed, 2),
    " selected_edges=", nrow(pair_df),
    " finite_pair_fraction=", round(finite_pair_fraction, 6),
    " lambda_best=", signif(result$lambda_best, 6)
  )
  out
}

safe_run <- function(project) {
  tryCatch(
    run_project(project),
    error = function(e) {
      files <- target_files(project)
      msg <- conditionMessage(e)
      try(cat("[ERROR] ", project, " ", msg, "\n",
              file = files$log, append = TRUE), silent = TRUE)
      out <- list(
        project = project,
        status = "error",
        error = msg,
        elapsed_seconds = NA_real_,
        n_samples = NA_integer_,
        n_taxa = NA_integer_,
        selected_edges = NA_integer_,
        finite_pair_fraction = NA_real_,
        lambda_best = NA_real_,
        full_zero_variance_taxa = NA_integer_,
        max_cv_test_zero_variance_taxa = NA_integer_,
        max_cv_train_zero_variance_taxa = NA_integer_
      )
      write_status(out, files$status)
      out
    }
  )
}

ensure_dirs()
projects <- list_projects()
if (length(projects) == 0L) stop("No projects found.")

message("Using official fastCCLasso core: ", opt$fast_core)
message("projects=", length(projects), " workers=", opt$workers)
message("nboot=", opt$nboot, " kmax=", opt$kmax,
        " pvalue=", opt$pvalue, " cor_cutoff=", opt$cor_cutoff)
message("For workers > 1, keep BLAS/OpenMP threads at 1.")

if (opt$workers == 1L) {
  results <- lapply(projects, safe_run)
} else {
  results <- mclapply(
    projects,
    safe_run,
    mc.cores = opt$workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

result_df <- do.call(
  rbind,
  lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE))
)

summary_file <- file.path(
  opt$output, "group_summaries",
  paste0("project_timing_and_status_", opt$run_label, ".tsv")
)
atomic_write_table(result_df, summary_file, sep = "\t", col_names = TRUE)

writeLines(
  result_df$project[result_df$status == "pass"],
  file.path(opt$output, "group_summaries", paste0("pass_", opt$run_label, ".list"))
)
writeLines(
  result_df$project[result_df$status == "error"],
  file.path(opt$output, "fastCCLasso", paste0("list_err_", opt$run_label, ".txt"))
)

error_df <- result_df[
  result_df$status == "error",
  c("project", "error"),
  drop = FALSE
]
atomic_write_table(
  error_df,
  file.path(opt$output, "group_summaries", paste0("errors_", opt$run_label, ".tsv")),
  sep = "\t", col_names = TRUE
)

print(table(result_df$status, useNA = "ifany"))
message("Summary: ", summary_file)
