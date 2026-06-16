#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Required R package not installed: data.table")
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Required R package not installed: digest")
  }
})

library(data.table)

log_msg <- function(...) {
  cat(sprintf(...), "\n", sep = "")
  flush.console()
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name)
  key_eq <- paste0(key, "=")

  hit_eq <- startsWith(args, key_eq)
  if (any(hit_eq)) {
    return(sub(key_eq, "", args[which(hit_eq)[1]], fixed = TRUE))
  }

  hit <- which(args == key)
  if (length(hit) > 0 && hit[1] < length(args)) {
    return(args[hit[1] + 1])
  }

  default
}

has_flag <- function(name) {
  args <- commandArgs(trailingOnly = TRUE)
  paste0("--", name) %in% args
}

INPUT_FILE <- get_arg(
  "input",
  "/home/minhao/knn/wgmlst_B_scheme/04_allelecall/target_foodenv_human_allelecall_contigs50_seed/results_alleles.tsv"
)

OUT_DIR <- get_arg(
  "outdir",
  "/home/minhao/knn/wgmlst_B_scheme/05_locus_filter/redundancy99"
)

METADATA_COLS_ARG <- get_arg("metadata-cols", "auto")
MISSING_THRESHOLD <- as.numeric(get_arg("missing-threshold", "0.99"))
MAJOR_THRESHOLD <- as.numeric(get_arg("major-threshold", "0.995"))
REDUNDANCY_THRESHOLD <- as.numeric(get_arg("redundancy-threshold", "0.99"))
N_CHUNKS <- as.integer(get_arg("n-chunks", "250"))
CANDIDATE_BATCH_SIZE <- as.integer(get_arg("candidate-batch-size", "512"))
EXPECTED_FINAL_LOCI <- as.integer(get_arg("expected-final-loci", "3586"))
WRITE_INTERMEDIATE_MATRICES <- has_flag("write-intermediate-matrices")

SPECIAL_MISSING <- c(
  "", "-", "NA", "N/A", "NAN", "NULL", "NONE",
  "LNF", "NIPH", "NIPHEM", "ALM", "ASM", "PLOT3", "PLOT5", "LOTSC", "PAMA"
)

detect_metadata_columns <- function(dt, metadata_cols_arg) {
  if (tolower(metadata_cols_arg) == "auto") {
    return(names(dt)[1])
  }

  cols <- trimws(unlist(strsplit(metadata_cols_arg, ",", fixed = TRUE)))
  cols <- cols[cols != ""]
  if (length(cols) == 0) {
    stop("No metadata columns were provided.")
  }

  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Metadata columns not found in input: ", paste(missing_cols, collapse = ", "))
  }

  cols
}

clean_allele_vector <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z)] <- ""
  zu <- toupper(z)

  missing_mask <- zu %in% SPECIAL_MISSING
  inf_mask <- grepl("^INF[-_]", zu)

  z2 <- z
  z2[inf_mask] <- sub("^INF[-_]", "", z2[inf_mask], ignore.case = TRUE)

  num <- suppressWarnings(as.integer(as.numeric(z2)))
  num[missing_mask | is.na(num) | num < 0L] <- -1L
  as.integer(num)
}

calculate_locus_stats <- function(X, loci) {
  n_samples <- nrow(X)
  n_loci <- ncol(X)

  stats <- vector("list", n_loci)

  for (j in seq_len(n_loci)) {
    col <- X[, j]
    valid <- col[col != -1L]
    n_valid <- length(valid)
    n_missing <- n_samples - n_valid
    missing_rate <- n_missing / n_samples
    call_rate <- n_valid / n_samples

    if (n_valid == 0L) {
      n_alleles <- 0L
      major_allele <- -1L
      major_count <- 0L
      major_allele_rate <- NA_real_
    } else {
      tab <- table(valid)
      max_pos <- which.max(tab)
      n_alleles <- length(tab)
      major_allele <- as.integer(names(tab)[max_pos])
      major_count <- as.integer(tab[max_pos])
      major_allele_rate <- major_count / n_valid
    }

    stats[[j]] <- data.table(
      locus = loci[j],
      original_order = j,
      n_samples = n_samples,
      n_valid = n_valid,
      n_missing = n_missing,
      missing_rate = missing_rate,
      call_rate = call_rate,
      n_alleles_nonmissing = n_alleles,
      major_allele = major_allele,
      major_allele_count = major_count,
      major_allele_rate = major_allele_rate
    )

    if (j %% 500 == 0 || j == n_loci) {
      log_msg("[stats] loci processed: %d/%d", j, n_loci)
    }
  }

  rbindlist(stats)
}

apply_basic_filters <- function(stats, missing_threshold, major_threshold) {
  stats[, pass_missing_rate := missing_rate <= missing_threshold]
  stats[, pass_polymorphic := n_alleles_nonmissing >= 2L]
  stats[, pass_major_allele := !is.na(major_allele_rate) & major_allele_rate <= major_threshold]
  stats[, pass_basic_locus_filters := pass_missing_rate & pass_polymorphic & pass_major_allele]
  stats
}

make_chunk_indices <- function(n_rows, n_chunks) {
  n_chunks <- max(1L, min(as.integer(n_chunks), as.integer(n_rows)))
  edges <- unique(as.integer(round(seq(1, n_rows + 1, length.out = n_chunks + 1))))
  chunks <- vector("list", length(edges) - 1L)

  for (i in seq_len(length(edges) - 1L)) {
    a <- edges[i]
    b <- edges[i + 1L] - 1L
    if (a <= b) {
      chunks[[i]] <- a:b
    }
  }

  chunks[!vapply(chunks, is.null, logical(1))]
}

hash_integer_vector <- function(x) {
  digest::digest(x, algo = "xxhash64", serialize = TRUE)
}

build_column_block_hashes <- function(X, chunks) {
  n_loci <- ncol(X)
  n_blocks <- length(chunks)

  hashes <- vector("list", n_loci)
  for (j in seq_len(n_loci)) {
    hashes[[j]] <- character(n_blocks)
  }

  for (b in seq_along(chunks)) {
    idx <- chunks[[b]]
    for (j in seq_len(n_loci)) {
      hashes[[j]][b] <- hash_integer_vector(X[idx, j])
    }
    log_msg("[redundancy] block hashes computed: %d/%d", b, n_blocks)
  }

  hashes
}

compute_similarity_to_candidates <- function(X, query_idx, candidate_indices, batch_size) {
  if (length(candidate_indices) == 0L) {
    return(list(best_idx = NA_integer_, best_sim = -1))
  }

  q <- X[, query_idx]
  best_idx <- NA_integer_
  best_sim <- -1

  starts <- seq(1L, length(candidate_indices), by = batch_size)
  for (s in starts) {
    e <- min(s + batch_size - 1L, length(candidate_indices))
    batch <- candidate_indices[s:e]

    sims <- colMeans(X[, batch, drop = FALSE] == q)
    local_pos <- which.max(sims)
    local_sim <- sims[local_pos]

    if (local_sim > best_sim) {
      best_sim <- as.numeric(local_sim)
      best_idx <- as.integer(batch[local_pos])
    }
  }

  list(best_idx = best_idx, best_sim = best_sim)
}

greedy_redundancy_filter <- function(
  X,
  loci,
  stats,
  similarity_threshold,
  n_chunks,
  batch_size
) {
  n_samples <- nrow(X)
  n_loci <- ncol(X)

  chunks <- make_chunk_indices(n_samples, n_chunks)
  log_msg(
    "[redundancy] samples=%d, loci_before_redundancy=%d, chunks=%d",
    n_samples, n_loci, length(chunks)
  )

  block_hashes <- build_column_block_hashes(X, chunks)

  order_dt <- copy(stats)
  order_dt[, sort_major_allele_rate := fifelse(is.na(major_allele_rate), 1.0, major_allele_rate)]
  setorder(order_dt, missing_rate, sort_major_allele_rate, -n_alleles_nonmissing, original_order)
  processing_order <- order_dt$matrix_index

  retained <- integer(0)
  removed <- list()
  lookup <- new.env(hash = TRUE, parent = emptyenv())

  for (step in seq_along(processing_order)) {
    j <- processing_order[step]
    candidates <- integer(0)

    hvec <- block_hashes[[j]]
    for (b in seq_along(hvec)) {
      key <- paste0(b, "|", hvec[b])
      if (exists(key, envir = lookup, inherits = FALSE)) {
        candidates <- c(candidates, get(key, envir = lookup, inherits = FALSE))
      }
    }

    if (length(candidates) > 0L) {
      candidates <- sort(unique(candidates))
    }

    sim_res <- compute_similarity_to_candidates(X, j, candidates, batch_size)
    best_idx <- sim_res$best_idx
    best_sim <- sim_res$best_sim

    if (!is.na(best_idx) && best_sim >= similarity_threshold) {
      removed[[length(removed) + 1L]] <- data.table(
        removed_locus = loci[j],
        removed_matrix_index = j,
        representative_locus = loci[best_idx],
        representative_matrix_index = best_idx,
        allele_pattern_similarity = best_sim,
        reason = paste0("similarity_ge_", similarity_threshold)
      )
    } else {
      retained <- c(retained, j)
      for (b in seq_along(hvec)) {
        key <- paste0(b, "|", hvec[b])
        if (exists(key, envir = lookup, inherits = FALSE)) {
          assign(key, c(get(key, envir = lookup, inherits = FALSE), j), envir = lookup)
        } else {
          assign(key, j, envir = lookup)
        }
      }
    }

    if (step == 1L || step %% 100L == 0L || step == length(processing_order)) {
      log_msg(
        "[redundancy] processed=%d/%d | retained=%d | removed=%d | candidates_last=%d",
        step, length(processing_order), length(retained), length(removed), length(candidates)
      )
    }
  }

  retained <- sort(retained)

  if (length(removed) == 0L) {
    removed_dt <- data.table(
      removed_locus = character(),
      removed_matrix_index = integer(),
      representative_locus = character(),
      representative_matrix_index = integer(),
      allele_pattern_similarity = numeric(),
      reason = character()
    )
  } else {
    removed_dt <- rbindlist(removed)
  }

  list(retained = retained, removed = removed_dt)
}

write_matrix <- function(path, metadata, X, loci) {
  out_dt <- as.data.table(X)
  setnames(out_dt, loci)
  out_dt <- cbind(metadata, out_dt)
  fwrite(out_dt, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(TRUE)
}

start_time <- Sys.time()

input_file <- normalizePath(INPUT_FILE, mustWork = FALSE)
out_dir <- normalizePath(OUT_DIR, mustWork = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("===== wgMLST locus filtering and redundancy reduction =====")
log_msg("[input] %s", input_file)
log_msg("[outdir] %s", out_dir)
log_msg("[params] missing_threshold=%s", MISSING_THRESHOLD)
log_msg("[params] major_threshold=%s", MAJOR_THRESHOLD)
log_msg("[params] redundancy_threshold=%s", REDUNDANCY_THRESHOLD)
log_msg("[params] n_chunks=%d", N_CHUNKS)
log_msg("[params] candidate_batch_size=%d", CANDIDATE_BATCH_SIZE)

log_msg("[step 1] reading allele table")
dt <- fread(input_file, sep = "\t", header = TRUE, colClasses = "character", data.table = TRUE, showProgress = TRUE)

metadata_cols <- detect_metadata_columns(dt, METADATA_COLS_ARG)
locus_cols <- setdiff(names(dt), metadata_cols)

if (length(locus_cols) == 0L) {
  stop("No locus columns detected.")
}

metadata <- copy(dt[, ..metadata_cols])
log_msg(
  "[input] rows=%d, columns=%d, metadata_cols=%s, locus_cols=%d",
  nrow(dt), ncol(dt), paste(metadata_cols, collapse = ","), length(locus_cols)
)

log_msg("[step 2] cleaning allele calls")
X <- matrix(-1L, nrow = nrow(dt), ncol = length(locus_cols))
colnames(X) <- locus_cols

for (j in seq_along(locus_cols)) {
  locus <- locus_cols[j]
  X[, j] <- clean_allele_vector(dt[[locus]])

  if (j %% 500L == 0L || j == length(locus_cols)) {
    log_msg("[clean] loci cleaned: %d/%d", j, length(locus_cols))
  }
}

rm(dt)
gc()

log_msg("[step 3] calculating locus statistics")
stats <- calculate_locus_stats(X, locus_cols)
stats <- apply_basic_filters(stats, MISSING_THRESHOLD, MAJOR_THRESHOLD)

stats_file <- file.path(out_dir, "wgmlst_locus_filter_stats.tsv")
fwrite(stats, stats_file, sep = "\t", quote = FALSE, na = "NA")

if (WRITE_INTERMEDIATE_MATRICES) {
  cleaned_file <- file.path(out_dir, "wgmlst_cleaned_all_loci.tsv")
  log_msg("[write] cleaned matrix: %s", cleaned_file)
  write_matrix(cleaned_file, metadata, X, locus_cols)
}

keep_basic_mask <- stats$pass_basic_locus_filters
basic_indices <- which(keep_basic_mask)

removed_basic <- copy(stats[!keep_basic_mask])
removed_basic_file <- file.path(out_dir, "wgmlst_removed_by_basic_filters.tsv")
fwrite(removed_basic, removed_basic_file, sep = "\t", quote = FALSE, na = "NA")

X_basic <- X[, basic_indices, drop = FALSE]
loci_basic <- locus_cols[basic_indices]
stats_basic <- copy(stats[keep_basic_mask])
stats_basic[, matrix_index := seq_len(.N)]

log_msg(
  "[basic filters] loci_before=%d, loci_after=%d, removed=%d",
  length(locus_cols), length(loci_basic), length(locus_cols) - length(loci_basic)
)
log_msg("[write] stats: %s", stats_file)
log_msg("[write] removed_by_basic_filters: %s", removed_basic_file)

if (WRITE_INTERMEDIATE_MATRICES) {
  basic_file <- file.path(out_dir, "wgmlst_filtered_missing99_major995.tsv")
  log_msg("[write] pre-redundancy matrix: %s", basic_file)
  write_matrix(basic_file, metadata, X_basic, loci_basic)
}

rm(X)
gc()

log_msg("[step 4] removing near-redundant loci")
redundancy_res <- greedy_redundancy_filter(
  X = X_basic,
  loci = loci_basic,
  stats = stats_basic,
  similarity_threshold = REDUNDANCY_THRESHOLD,
  n_chunks = N_CHUNKS,
  batch_size = CANDIDATE_BATCH_SIZE
)

retained_basic_indices <- redundancy_res$retained
redundancy_map <- redundancy_res$removed

X_final <- X_basic[, retained_basic_indices, drop = FALSE]
loci_final <- loci_basic[retained_basic_indices]

retained_loci_file <- file.path(out_dir, "wgmlst_retained_loci_redundancy99.txt")
removed_redundant_file <- file.path(out_dir, "wgmlst_removed_redundant_loci_redundancy99.tsv")
final_file <- file.path(out_dir, "wgmlst_filtered_missing99_major995_redundancy99.tsv")
summary_file <- file.path(out_dir, "wgmlst_filter_redundancy99_summary.tsv")

writeLines(loci_final, retained_loci_file)
fwrite(redundancy_map, removed_redundant_file, sep = "\t", quote = FALSE, na = "NA")

log_msg("[write] final matrix: %s", final_file)
write_matrix(final_file, metadata, X_final, loci_final)

elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

summary_dt <- data.table(
  metric = c(
    "input_file",
    "output_directory",
    "n_samples",
    "n_metadata_columns",
    "metadata_columns",
    "n_loci_raw",
    "n_loci_after_basic_filters",
    "n_loci_removed_by_basic_filters",
    "n_loci_removed_by_redundancy99",
    "n_loci_final",
    "missing_threshold",
    "major_allele_threshold",
    "redundancy_similarity_threshold",
    "expected_final_loci",
    "final_locus_count_matches_expected",
    "final_matrix",
    "retained_loci_file",
    "removed_redundant_loci_file",
    "locus_stats_file",
    "removed_by_basic_filters_file",
    "elapsed_seconds"
  ),
  value = as.character(c(
    input_file,
    out_dir,
    nrow(metadata),
    length(metadata_cols),
    paste(metadata_cols, collapse = ","),
    length(locus_cols),
    length(loci_basic),
    length(locus_cols) - length(loci_basic),
    length(loci_basic) - length(loci_final),
    length(loci_final),
    MISSING_THRESHOLD,
    MAJOR_THRESHOLD,
    REDUNDANCY_THRESHOLD,
    EXPECTED_FINAL_LOCI,
    length(loci_final) == EXPECTED_FINAL_LOCI,
    final_file,
    retained_loci_file,
    removed_redundant_file,
    stats_file,
    removed_basic_file,
    elapsed_seconds
  ))
)

fwrite(summary_dt, summary_file, sep = "\t", quote = FALSE, na = "NA")

log_msg("===== summary =====")
for (i in seq_len(nrow(summary_dt))) {
  log_msg("%s: %s", summary_dt$metric[i], summary_dt$value[i])
}

if (EXPECTED_FINAL_LOCI > 0L && length(loci_final) != EXPECTED_FINAL_LOCI) {
  log_msg(
    "[WARNING] final locus count is %d, but expected %d. Check input data and redundancy settings.",
    length(loci_final), EXPECTED_FINAL_LOCI
  )
}

log_msg("[DONE]")
