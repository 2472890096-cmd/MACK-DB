#!/usr/bin/env Rscript

# ============================================================
# DKI Structural Keystoneness calculation for new_DKI_rk4_fast
# - Reads original Ptrain/Ptest/species_id/sample_id from old project data/
# - Reads new qtst/qtrn from new_DKI_rk4_fast/<PROJECT>/results/
# - Writes Ks results to new_DKI_rk4_fast/<PROJECT>/results/
# ============================================================

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

# -------------------- argument parser --------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

get_bool_arg <- function(flag, default = FALSE) {
  val <- get_arg(flag, NA)
  if (is.na(val)) return(default)
  return(tolower(val) %in% c("true", "t", "1", "yes", "y"))
}

project_list_file <- get_arg(
  "--project-list",
  "data/dki/projects.txt"
)

data_base_dir <- get_arg(
  "--data-base",
  "data/dki"
)

pred_base_dir <- get_arg(
  "--pred-base",
  "work/dki"
)

threads <- as.integer(get_arg("--threads", "4"))
force <- get_bool_arg("--force", FALSE)

if (is.na(threads) || threads < 1) threads <- 1

status_file <- file.path(pred_base_dir, "reports", "ks_rk4_fast_status.tsv")
dir.create(file.path(pred_base_dir, "reports"), recursive = TRUE, showWarnings = FALSE)

# -------------------- helpers --------------------
read_matrix_csv <- function(file) {
  if (!file.exists(file)) stop(sprintf("Missing file: %s", file))

  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- as.matrix(data.table::fread(file, header = FALSE, data.table = FALSE))
  } else {
    x <- as.matrix(read.table(file, header = FALSE, sep = ",", check.names = FALSE))
  }

  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  return(x)
}

safe_norm_cols <- function(x, eps = 1e-12) {
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  cs <- colSums(x)
  valid <- cs > eps
  out <- x
  if (any(valid)) {
    out[, valid] <- sweep(out[, valid, drop = FALSE], 2, cs[valid], "/")
  }
  if (any(!valid)) {
    out[, !valid] <- 0
  }
  return(list(mat = out, valid = valid, sums = cs))
}

safe_norm_rows <- function(x, eps = 1e-12) {
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  rs <- rowSums(x)
  valid <- rs > eps
  out <- x
  if (any(valid)) {
    out[valid, ] <- sweep(out[valid, , drop = FALSE], 1, rs[valid], "/")
  }
  if (any(!valid)) {
    out[!valid, ] <- 0
  }
  return(list(mat = out, valid = valid, sums = rs))
}

append_reason <- function(reason_vec, idx, reason) {
  if (!any(idx)) return(reason_vec)
  reason_vec[idx] <- ifelse(
    reason_vec[idx] == "",
    reason,
    paste(reason_vec[idx], reason, sep = ";")
  )
  return(reason_vec)
}

detect_index_base_and_adjust <- function(idx, max_allowed, name) {
  idx <- as.integer(idx)

  if (length(idx) == 0) {
    stop(sprintf("%s is empty", name))
  }

  if (any(is.na(idx))) {
    stop(sprintf("%s contains NA", name))
  }

  if (min(idx) == 0) {
    idx_adj <- idx + 1L
    index_base <- "0-based converted to 1-based"
  } else {
    idx_adj <- idx
    index_base <- "1-based"
  }

  if (min(idx_adj) < 1 || max(idx_adj) > max_allowed) {
    stop(sprintf(
      "%s index out of range after adjustment: min=%d max=%d allowed=1..%d",
      name, min(idx_adj), max(idx_adj), max_allowed
    ))
  }

  return(list(idx = idx_adj, index_base = index_base))
}

write_empty_outputs <- function(out_file1, out_file2, out_file3) {
  empty_detail <- data.frame(
    species = integer(),
    sample = integer(),
    original_abundance = numeric(),
    BC_pred = numeric(),
    Ks_pred = numeric(),
    BC_reference = numeric(),
    Ks_reference = numeric(),
    Ks_true = numeric(),
    valid = logical(),
    invalid_reason = character()
  )
  write.csv(empty_detail, out_file1, row.names = FALSE)

  empty_summary <- data.frame(
    species = integer(),
    n_valid = integer(),
    median_Ks_pred = numeric(),
    mean_Ks_pred = numeric(),
    q90_Ks_pred = numeric(),
    max_Ks_pred = numeric(),
    median_BC_pred = numeric(),
    mean_BC_pred = numeric(),
    rank_median_Ks_pred = integer(),
    rank_q90_Ks_pred = integer(),
    rank_max_Ks_pred = integer()
  )
  write.csv(empty_summary, out_file2, row.names = FALSE)

  empty_legacy <- data.frame(
    species = integer(),
    Ks_true = numeric(),
    Ks_pred = numeric()
  )
  write.csv(empty_legacy, out_file3, row.names = FALSE)
}

# -------------------- core calculation --------------------
process_project <- function(project) {
  start_time <- Sys.time()

  data_dir <- file.path(data_base_dir, project, "data")
  result_dir <- file.path(pred_base_dir, project, "results")

  f_sp   <- file.path(data_dir, "species_id.csv")
  f_sa   <- file.path(data_dir, "sample_id.csv")
  f_ptrn <- file.path(data_dir, "Ptrain.csv")
  f_ptst <- file.path(data_dir, "Ptest.csv")

  f_qtst <- file.path(result_dir, "qtst.csv")
  f_qtrn <- file.path(result_dir, "qtrn.csv")

  out_file_detail <- file.path(result_dir, "Keystoneness_structural.csv")
  out_file_summary <- file.path(result_dir, "Keystoneness_structural_summary.csv")
  out_file_legacy_median <- file.path(result_dir, "Keystoneness_structural_median.csv")

  req_files <- c(f_sp, f_sa, f_ptrn, f_ptst, f_qtst, f_qtrn)
  missing <- req_files[!file.exists(req_files)]

  if (length(missing) > 0) {
    return(data.frame(
      project = project,
      status = "MISSING_INPUT",
      n_rows = NA_integer_,
      n_valid = NA_integer_,
      n_invalid = NA_integer_,
      n_species = NA_integer_,
      n_samples = NA_integer_,
      median_Ks_pred_global = NA_real_,
      max_Ks_pred_global = NA_real_,
      output_detail = out_file_detail,
      output_summary = out_file_summary,
      message = paste("Missing:", paste(missing, collapse = ";")),
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      stringsAsFactors = FALSE
    ))
  }

  if (!force && file.exists(out_file_detail) && file.exists(out_file_summary) && file.exists(out_file_legacy_median)) {
    return(data.frame(
      project = project,
      status = "SKIP",
      n_rows = NA_integer_,
      n_valid = NA_integer_,
      n_invalid = NA_integer_,
      n_species = NA_integer_,
      n_samples = NA_integer_,
      median_Ks_pred_global = NA_real_,
      max_Ks_pred_global = NA_real_,
      output_detail = out_file_detail,
      output_summary = out_file_summary,
      message = "existing outputs",
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

    species_id <- read.table(f_sp, header = FALSE, sep = ",")
    sample_id  <- read.table(f_sa, header = FALSE, sep = ",")

    Ptrain_raw <- read_matrix_csv(f_ptrn)
    Ptest_raw  <- read_matrix_csv(f_ptst)
    qtst_raw   <- read_matrix_csv(f_qtst)
    qtrn_raw   <- read_matrix_csv(f_qtrn)

    # Normalize matrices safely
    Ptrain_norm <- safe_norm_cols(Ptrain_raw)
    Ptest_norm  <- safe_norm_cols(Ptest_raw)
    qtst_norm   <- safe_norm_rows(qtst_raw)
    qtrn_norm   <- safe_norm_rows(qtrn_raw)

    Ptrain <- Ptrain_norm$mat
    ptst   <- Ptest_norm$mat
    qtst   <- qtst_norm$mat
    qtrn   <- qtrn_norm$mat

    n_species <- nrow(Ptrain)
    n_samples <- ncol(Ptrain)

    if (nrow(ptst) != n_species) {
      stop(sprintf("Ptest species rows mismatch: Ptrain=%d Ptest=%d", n_species, nrow(ptst)))
    }
    if (ncol(qtrn) != n_species) {
      stop(sprintf("qtrn species columns mismatch: Ptrain species=%d qtrn columns=%d", n_species, ncol(qtrn)))
    }
    if (ncol(qtst) != n_species) {
      stop(sprintf("qtst species columns mismatch: Ptrain species=%d qtst columns=%d", n_species, ncol(qtst)))
    }

    N_exp <- nrow(qtst)

    if (ncol(ptst) != N_exp) {
      stop(sprintf("Ptest perturbation count mismatch: ncol(Ptest)=%d nrow(qtst)=%d", ncol(ptst), N_exp))
    }

    if (nrow(species_id) != N_exp || nrow(sample_id) != N_exp) {
      stop(sprintf(
        "species_id/sample_id length mismatch: species_id=%d sample_id=%d qtst_rows=%d",
        nrow(species_id), nrow(sample_id), N_exp
      ))
    }

    sp_check <- detect_index_base_and_adjust(species_id$V1, n_species, "species_id")
    sa_check <- detect_index_base_and_adjust(sample_id$V1, nrow(qtrn), "sample_id")

    sp_idx <- sp_check$idx
    s_idx  <- sa_check$idx

    idx_mat <- cbind(seq_len(N_exp), sp_idx)

    # Predicted null composition
    q_null_mat <- qtrn[s_idx, , drop = FALSE]
    q_null_mat[idx_mat] <- 0
    q_null_rs <- rowSums(q_null_mat)
    q_null_valid <- q_null_rs > 1e-12
    if (any(q_null_valid)) {
      q_null_mat[q_null_valid, ] <- sweep(
        q_null_mat[q_null_valid, , drop = FALSE],
        1,
        q_null_rs[q_null_valid],
        "/"
      )
    }
    if (any(!q_null_valid)) {
      q_null_mat[!q_null_valid, ] <- 0
    }

    # Reference null composition
    p_null_mat <- t(Ptrain)[s_idx, , drop = FALSE]
    p_null_mat[idx_mat] <- 0
    p_null_rs <- rowSums(p_null_mat)
    p_null_valid <- p_null_rs > 1e-12
    if (any(p_null_valid)) {
      p_null_mat[p_null_valid, ] <- sweep(
        p_null_mat[p_null_valid, , drop = FALSE],
        1,
        p_null_rs[p_null_valid],
        "/"
      )
    }
    if (any(!p_null_valid)) {
      p_null_mat[!p_null_valid, ] <- 0
    }

    ptst_t <- t(ptst)

    qtst_valid <- rowSums(qtst) > 1e-12
    ptst_valid <- rowSums(ptst_t) > 1e-12

    # Bray-Curtis predicted
    BC_pred <- rep(NA_real_, N_exp)
    bc_pred_den <- rowSums(abs(q_null_mat + qtst))
    bc_pred_num <- rowSums(abs(q_null_mat - qtst))
    ok_pred <- q_null_valid & qtst_valid & bc_pred_den > 1e-12
    BC_pred[ok_pred] <- bc_pred_num[ok_pred] / bc_pred_den[ok_pred]

    # Bray-Curtis reference
    BC_reference <- rep(NA_real_, N_exp)
    bc_ref_den <- rowSums(abs(p_null_mat + ptst_t))
    bc_ref_num <- rowSums(abs(p_null_mat - ptst_t))
    ok_ref <- p_null_valid & ptst_valid & bc_ref_den > 1e-12
    BC_reference[ok_ref] <- bc_ref_num[ok_ref] / bc_ref_den[ok_ref]

    original_abundance <- Ptrain[cbind(sp_idx, s_idx)]

    Ks_pred <- BC_pred * (1 - original_abundance)
    Ks_reference <- BC_reference * (1 - original_abundance)

    valid <- is.finite(Ks_pred)

    invalid_reason <- rep("", N_exp)
    invalid_reason <- append_reason(invalid_reason, !q_null_valid, "q_null_zero_after_removal")
    invalid_reason <- append_reason(invalid_reason, !qtst_valid, "qtst_zero_row")
    invalid_reason <- append_reason(invalid_reason, bc_pred_den <= 1e-12, "bc_pred_den_zero")
    invalid_reason <- append_reason(invalid_reason, !is.finite(Ks_pred), "Ks_pred_nonfinite")
    invalid_reason[valid & invalid_reason == ""] <- "OK"

    detail <- data.frame(
      species = sp_idx,
      sample = s_idx,
      original_abundance = original_abundance,
      BC_pred = BC_pred,
      Ks_pred = Ks_pred,
      BC_reference = BC_reference,
      Ks_reference = Ks_reference,
      Ks_true = Ks_reference,   # legacy-compatible name; reference only
      valid = valid,
      invalid_reason = invalid_reason,
      stringsAsFactors = FALSE
    )

    write.csv(detail, out_file_detail, row.names = FALSE)

    valid_df <- detail %>%
      filter(valid, is.finite(Ks_pred))

    if (nrow(valid_df) == 0) {
      write_empty_outputs(out_file_detail, out_file_summary, out_file_legacy_median)
      status <- "NO_VALID_ROWS"
      med_global <- NA_real_
      max_global <- NA_real_
      n_valid <- 0L
      n_invalid <- N_exp
    } else {
      summary_df <- valid_df %>%
        group_by(species) %>%
        summarise(
          n_valid = n(),
          median_Ks_pred = median(Ks_pred, na.rm = TRUE),
          mean_Ks_pred = mean(Ks_pred, na.rm = TRUE),
          q90_Ks_pred = as.numeric(quantile(Ks_pred, 0.9, na.rm = TRUE)),
          max_Ks_pred = max(Ks_pred, na.rm = TRUE),
          median_BC_pred = median(BC_pred, na.rm = TRUE),
          mean_BC_pred = mean(BC_pred, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          rank_median_Ks_pred = rank(-median_Ks_pred, ties.method = "min"),
          rank_q90_Ks_pred = rank(-q90_Ks_pred, ties.method = "min"),
          rank_max_Ks_pred = rank(-max_Ks_pred, ties.method = "min")
        ) %>%
        arrange(rank_median_Ks_pred, rank_q90_Ks_pred)

      write.csv(summary_df, out_file_summary, row.names = FALSE)

      legacy_median <- valid_df %>%
        group_by(species) %>%
        summarise(
          Ks_true = median(Ks_reference, na.rm = TRUE),
          Ks_pred = median(Ks_pred, na.rm = TRUE),
          .groups = "drop"
        )

      write.csv(legacy_median, out_file_legacy_median, row.names = FALSE)

      status <- "DONE"
      med_global <- median(valid_df$Ks_pred, na.rm = TRUE)
      max_global <- max(valid_df$Ks_pred, na.rm = TRUE)
      n_valid <- nrow(valid_df)
      n_invalid <- N_exp - n_valid
    }

    return(data.frame(
      project = project,
      status = status,
      n_rows = N_exp,
      n_valid = n_valid,
      n_invalid = n_invalid,
      n_species = n_species,
      n_samples = n_samples,
      median_Ks_pred_global = med_global,
      max_Ks_pred_global = max_global,
      output_detail = out_file_detail,
      output_summary = out_file_summary,
      message = sprintf(
        "sample_index=%s; species_index=%s",
        sa_check$index_base,
        sp_check$index_base
      ),
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    return(data.frame(
      project = project,
      status = "FAILED",
      n_rows = NA_integer_,
      n_valid = NA_integer_,
      n_invalid = NA_integer_,
      n_species = NA_integer_,
      n_samples = NA_integer_,
      median_Ks_pred_global = NA_real_,
      max_Ks_pred_global = NA_real_,
      output_detail = out_file_detail,
      output_summary = out_file_summary,
      message = e$message,
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      stringsAsFactors = FALSE
    ))
  })
}

# -------------------- main --------------------
if (!file.exists(project_list_file)) {
  stop(sprintf("Project list not found: %s", project_list_file))
}

projects <- readLines(project_list_file)
projects <- trimws(projects)
projects <- projects[projects != ""]

if (length(projects) == 0) {
  stop("No projects found in project list.")
}

max_cores <- max(1, parallel::detectCores() - 2)
if (threads > max_cores) threads <- max_cores
if (threads < 1) threads <- 1

cat(sprintf("Project list: %s\n", project_list_file))
cat(sprintf("Data base: %s\n", data_base_dir))
cat(sprintf("Prediction base: %s\n", pred_base_dir))
cat(sprintf("Projects: %d\n", length(projects)))
cat(sprintf("Threads: %d\n", threads))
cat(sprintf("Force: %s\n", force))

if (threads == 1) {
  results <- lapply(projects, process_project)
} else {
  results <- parallel::mclapply(projects, process_project, mc.cores = threads)
}

status_df <- bind_rows(results)

write.table(
  status_df,
  status_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("\nStatus written to: %s\n", status_file))
cat("\nStatus summary:\n")
print(table(status_df$status, useNA = "ifany"))

cat("\n🎉 Structural Keystoneness calculation finished.\n")
