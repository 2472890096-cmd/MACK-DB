#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    abundance = NULL,
    species = NULL,
    output = NULL,
    outdir = ".",
    mode = "all",
    prevalence_min = 0.25,
    prevalence_max = 0.75,
    network_quantile = 0.1
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value_keys <- c(
      "-a", "--abundance", "-s", "--species", "-o", "--output",
      "-d", "--outdir", "--mode", "--prevalence-min",
      "--prevalence-max", "--network-quantile"
    )

    if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript get_keystone_EPI_single_safe_qecdf_v3_fast.R ",
        "-a abundance.matrix -s species.matrix -o project_id -d outdir ",
        "--mode all|prevalence [--prevalence-min 0.25 --prevalence-max 0.75]\n",
        sep = ""
      )
      quit(status = 0)
    }

    if (!(key %in% value_keys) || i == length(args)) {
      stop("Invalid or incomplete argument: ", key)
    }

    value <- args[[i + 1]]

    if (key %in% c("-a", "--abundance")) opts$abundance <- value
    if (key %in% c("-s", "--species")) opts$species <- value
    if (key %in% c("-o", "--output")) opts$output <- value
    if (key %in% c("-d", "--outdir")) opts$outdir <- value
    if (key == "--mode") opts$mode <- value
    if (key == "--prevalence-min") opts$prevalence_min <- as.numeric(value)
    if (key == "--prevalence-max") opts$prevalence_max <- as.numeric(value)
    if (key == "--network-quantile") opts$network_quantile <- as.numeric(value)

    i <- i + 2
  }

  required <- c("abundance", "species", "output")
  missing <- required[vapply(opts[required], is.null, logical(1))]

  if (length(missing) > 0) {
    stop("Missing required option(s): ", paste(missing, collapse = ", "))
  }

  if (!(opts$mode %in% c("all", "prevalence"))) {
    stop("--mode must be all or prevalence")
  }

  if (!is.finite(opts$prevalence_min) ||
      !is.finite(opts$prevalence_max)) {
    stop("Prevalence bounds must be numeric")
  }

  if (opts$prevalence_min < 0 ||
      opts$prevalence_max > 1 ||
      opts$prevalence_min > opts$prevalence_max) {
    stop("Invalid prevalence bounds")
  }

  if (!is.finite(opts$network_quantile) ||
      opts$network_quantile < 0 ||
      opts$network_quantile > 1) {
    stop("--network-quantile must be between 0 and 1")
  }

  opts
}


read_abundance <- function(path) {
  otu <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  otu <- as.matrix(otu)
  storage.mode(otu) <- "numeric"
  otu[!is.finite(otu)] <- 0

  otu
}


read_species <- function(path, expected_n) {
  sp <- read.table(
    path,
    sep = "\t",
    header = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (ncol(sp) >= 2) {
    species <- sp[[2]]
  } else {
    species <- sp[[1]]
  }

  species <- as.character(species)

  if (length(species) != expected_n) {
    stop(
      "Species count does not match abundance rows: species=",
      length(species),
      ", abundance_rows=",
      expected_n
    )
  }

  species
}


normalize_columns <- function(mat) {
  mat <- as.matrix(mat)

  col_total <- colSums(mat)

  out <- mat

  positive <- is.finite(col_total) & col_total > 0

  if (any(positive)) {
    out[, positive] <- sweep(
      mat[, positive, drop = FALSE],
      2,
      col_total[positive],
      "/"
    )
  }

  if (any(!positive)) {
    out[, !positive] <- 0
  }

  out[!is.finite(out)] <- 0

  out
}


bray_two <- function(a, b) {
  denom <- sum(a + b)

  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }

  value <- sum(abs(a - b)) / denom

  if (is.finite(value)) {
    value
  } else {
    0
  }
}


safe_bray_dist <- function(species_by_sample) {
  sample_by_species <- t(species_by_sample)

  bc_dist <- vegan::vegdist(
    sample_by_species,
    method = "bray"
  )

  bad <- !is.finite(bc_dist)

  if (any(bad)) {
    bc_dist[bad] <- 0
  }

  pair_dist <- as.numeric(bc_dist)

  bc <- as.matrix(bc_dist)

  bc[!is.finite(bc)] <- 0

  bc <- (bc + t(bc)) / 2

  diag(bc) <- 0

  list(
    matrix = bc,
    pair_dist = pair_dist
  )
}


ecdf_threshold_fast <- function(pair_dist, p) {
  x <- pair_dist

  if (length(x) == 0) {
    return(NA_real_)
  }

  finite <- is.finite(x)

  if (!all(finite)) {
    x <- x[finite]
  }

  n <- length(x)

  if (n == 0) {
    return(NA_real_)
  }

  if (p <= 0) {
    return(min(x))
  }

  limit <- p * n

  tol <- 64 * .Machine$double.eps * max(1, abs(limit))

  k <- floor(limit + tol)

  if (k < 1) {
    return(min(x))
  }

  if (k >= n) {
    return(max(x))
  }

  partial_x <- sort(
    x,
    partial = k
  )

  candidate <- partial_x[k]

  cumulative_count <- sum(x <= candidate)

  if (cumulative_count <= limit + tol) {
    return(candidate)
  }

  lower_mask <- x < candidate

  if (!any(lower_mask)) {
    return(min(x))
  }

  max(x[lower_mask])
}


select_outliers <- function(values) {
  ok <- is.finite(values)

  if (sum(ok) < 2) {
    return(rep(FALSE, length(values)))
  }

  threshold <- mean(values[ok]) +
    2 * stats::sd(values[ok])

  if (!is.finite(threshold)) {
    return(rep(FALSE, length(values)))
  }

  ok & values > threshold
}


write_key <- function(values, path) {
  values <- as.character(
    values[!is.na(values)]
  )

  con <- file(
    path,
    open = "w"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  if (length(values) > 0) {
    writeLines(
      values,
      con = con
    )
  }
}


opts <- parse_args(
  commandArgs(
    trailingOnly = TRUE
  )
)


if (!requireNamespace(
  "vegan",
  quietly = TRUE
)) {
  stop(
    "Required R package 'vegan' is not installed. ",
    "Install it first, for example: install.packages('vegan')"
  )
}


dir.create(
  opts$outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


otu <- read_abundance(
  opts$abundance
)

species <- read_species(
  opts$species,
  nrow(otu)
)


N <- nrow(otu)

M <- ncol(otu)


one_zero_otu <- ifelse(
  otu > 0,
  1,
  0
)


prevalence <- rowSums(
  one_zero_otu,
  na.rm = TRUE
) / M


D1 <- rep(
  NA_real_,
  N
)

D2 <- rep(
  NA_real_,
  N
)

Q <- rep(
  NA_real_,
  N
)


cat(
  sprintf(
    "Project: %s\n",
    opts$output
  )
)

cat(
  sprintf(
    "Mode: %s\n",
    opts$mode
  )
)

cat(
  sprintf(
    "Species: %d; Samples: %d\n",
    N,
    M
  )
)


for (i in seq_len(N)) {

  present_count <- sum(
    one_zero_otu[i, ],
    na.rm = TRUE
  )

  if (opts$mode == "all") {

    should_compute <-
      present_count > 0 &&
      present_count < M

  } else {

    should_compute <-
      prevalence[i] >= opts$prevalence_min &&
      prevalence[i] <= opts$prevalence_max
  }


  if (!should_compute) {
    next
  }


  ind_pres <-
    one_zero_otu[i, ] != 0


  out2 <- normalize_columns(
    otu[-i, , drop = FALSE]
  )


  bray_result <- safe_bray_dist(
    out2
  )


  bc <- bray_result$matrix

  pair_dist <-
    bray_result$pair_dist


  bc2 <- bc[
    ind_pres,
    !ind_pres,
    drop = FALSE
  ]


  finite_bc2 <-
    bc2[is.finite(bc2)]


  if (length(finite_bc2) > 0) {
    D1[i] <-
      mean(finite_bc2)
  }


  threshold <- ecdf_threshold_fast(
    pair_dist,
    opts$network_quantile
  )


  if (is.finite(threshold)) {

    B_i <- bc <= threshold

    B_i[is.na(B_i)] <- FALSE

    B_i <- B_i | t(B_i)

    diag(B_i) <- FALSE


    s_i <- as.numeric(
      one_zero_otu[i, ]
    )

    s_i[s_i == 0] <- -1


    d <- rowSums(B_i)

    q_sum <- sum(B_i) / 2


    if (is.finite(q_sum) &&
        q_sum > 0) {

      expected <-
        (d %*% t(d)) /
        (2 * q_sum)


      Q[i] <- as.numeric(
        (
          t(s_i) %*%
            (B_i - expected) %*%
            s_i
        ) /
          (4 * q_sum)
      )

    } else {

      Q[i] <- 0
    }

  } else {

    Q[i] <- 0
  }


  otu_pres <- normalize_columns(
    otu[
      -i,
      ind_pres,
      drop = FALSE
    ]
  )


  otu_abs <- normalize_columns(
    otu[
      -i,
      !ind_pres,
      drop = FALSE
    ]
  )


  D2[i] <- bray_two(
    rowMeans(otu_pres),
    rowMeans(otu_abs)
  )
}


B1 <- select_outliers(D1)

B2 <- select_outliers(D2)

BQ <- select_outliers(Q)


write_key(
  species[which(B1)],
  file.path(
    opts$outdir,
    paste0(
      opts$output,
      "_EPI-D1.key"
    )
  )
)


write_key(
  species[which(B2)],
  file.path(
    opts$outdir,
    paste0(
      opts$output,
      "_EPI-D2.key"
    )
  )
)


write_key(
  species[which(BQ)],
  file.path(
    opts$outdir,
    paste0(
      opts$output,
      "_EPI-Q.key"
    )
  )
)


summary_path <- file.path(
  opts$outdir,
  paste0(
    opts$output,
    "_EPI_summary.tsv"
  )
)


summary_df <- data.frame(
  project = opts$output,
  mode = opts$mode,
  species_total = N,
  sample_total = M,
  species_scored =
    sum(
      is.finite(D1) |
      is.finite(D2) |
      is.finite(Q)
    ),
  D1_key_count =
    sum(B1),
  D2_key_count =
    sum(B2),
  Q_key_count =
    sum(BQ),
  stringsAsFactors = FALSE
)


write.table(
  summary_df,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  sprintf(
    "Done: D1=%d, D2=%d, Q=%d\n",
    summary_df$D1_key_count,
    summary_df$D2_key_count,
    summary_df$Q_key_count
  )
)
