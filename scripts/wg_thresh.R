#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

setDTthreads(30)

# ============================================================
# Configuration
# ============================================================
tag <- "wg"
dataset_label <- "wgMLST-Top800"

train_file <- "/home/minhao/knn/xinxinxinknn/data/wgmlst/xunlianji/top/top_matrices/wgmlst_top800_knn_hamming_importance_R.with_label.tsv"
model_file <- "/home/minhao/knn/xinxinxinknn/moxing/knn/wg/moxing/knn_hamming_top800_k11_inv1_model.rds"

out_root <- "/home/minhao/knn/xinxinxinknn/data/yuzhi"
raw_out_dir <- file.path(out_root, "50wg")
normalized_out_dir <- file.path(out_root, "50wgnormalized")
dir.create(raw_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normalized_out_dir, recursive = TRUE, showWarnings = FALSE)

# Optional out-of-fold prediction files. If one exists and use_correct_oof_only=TRUE,
# thresholds are estimated only from correctly classified training isolates.
use_correct_oof_only <- FALSE
oof_file_candidates <- c(
  "/home/minhao/knn/xinxinxinknn/moxing/knn/wg/moxing/oof_probs_top800.tsv",
  "/home/minhao/knn/xinxinxinknn/moxing/knn/wgmlst_ranked_topN/oof_predictions_wgmlst_top800.tsv",
  "/home/minhao/knn/xinxinxinknn/moxing/knn/wg/paixu/oof_predictions_wgmlst_top800.tsv"
)

# Optional CC filtering. Set require_cc_ge50=TRUE and provide cc_file to calculate
# thresholds only from samples belonging to clonal complexes with at least 50 isolates.
require_cc_ge50 <- FALSE
cc_file <- ""
cc_min_count <- 50L
cc_sample_col <- "sample"
cc_col <- "CC"

# If TRUE, write an extra meanTop11 normalized threshold file using the legacy
# CCge50-style name when CC filtering is not active. Keep FALSE unless you need
# to feed an older downstream filtering script without editing its threshold path.
write_legacy_ccge50_filename <- FALSE

k_top <- 11L
missing_code <- -1L
n_cores <- 30L
chunk_size <- 25L
expected_features <- 800L

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

logf <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf(...))
  cat("\n")
  flush.console()
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("-", "_", x)
}

clean_allele_vec <- function(x, missing_code = -1L) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^INF[-_]", "", x)

  non_call <- c(
    "", "-", "NA", "N/A", "NaN", "nan", "NULL", "null",
    "LNF", "NIPH", "NIPHEM", "ALM", "ASM",
    "PLOT3", "PLOT5", "LOTSC", "PAMA",
    "INF", "Inf", "inf"
  )

  x[x %in% non_call] <- NA_character_
  v <- suppressWarnings(as.integer(as.numeric(x)))
  v[is.na(v)] <- missing_code
  as.integer(v)
}

calc_stats <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(data.table(
      n_valid = 0L,
      min_value = NA_real_,
      q1_value = NA_real_,
      median_value = NA_real_,
      mean_value = NA_real_,
      sd_value = NA_real_,
      q90_value = NA_real_,
      q95_value = NA_real_,
      q975_value = NA_real_,
      max_value = NA_real_
    ))
  }

  data.table(
    n_valid = length(x),
    min_value = min(x),
    q1_value = as.numeric(quantile(x, 0.25, names = FALSE, type = 7)),
    median_value = median(x),
    mean_value = mean(x),
    sd_value = ifelse(length(x) > 1, sd(x), NA_real_),
    q90_value = as.numeric(quantile(x, 0.90, names = FALSE, type = 7)),
    q95_value = as.numeric(quantile(x, 0.95, names = FALSE, type = 7)),
    q975_value = as.numeric(quantile(x, 0.975, names = FALSE, type = 7)),
    max_value = max(x)
  )
}

find_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

infer_predicted_from_probs <- function(oof_dt, label_levels) {
  prob_cols <- intersect(label_levels, names(oof_dt))
  if (length(prob_cols) == 0) {
    return(rep(NA_character_, nrow(oof_dt)))
  }
  prob_mat <- as.matrix(oof_dt[, ..prob_cols])
  prob_mat <- apply(prob_mat, 2, as.numeric)
  label_levels[max.col(prob_mat, ties.method = "first")]
}

read_oof_correct_samples <- function(oof_file, train_dt, label_levels) {
  if (is.na(oof_file) || !file.exists(oof_file)) {
    return(NULL)
  }

  logf("Reading OOF predictions: %s", oof_file)
  oof <- fread(oof_file)
  if (!"sample" %in% names(oof)) {
    stop("OOF file must contain a sample column: ", oof_file)
  }
  oof[, sample := as.character(sample)]
  oof[, sample_clean := clean_id(sample)]

  if (!"true_label" %in% names(oof)) {
    if ("label" %in% names(oof)) {
      setnames(oof, "label", "true_label")
    } else {
      oof <- merge(
        oof,
        unique(train_dt[, .(sample_clean, true_label = as.character(label))]),
        by = "sample_clean",
        all.x = TRUE
      )
    }
  }

  if (!"predicted_label" %in% names(oof)) {
    oof[, predicted_label := infer_predicted_from_probs(.SD, label_levels), .SDcols = names(oof)]
  }

  if (anyDuplicated(oof$sample_clean)) {
    dup <- oof[duplicated(sample_clean) | duplicated(sample_clean, fromLast = TRUE)]
    dup_file <- file.path(raw_out_dir, paste0(tag, "_duplicated_oof_sample_clean_ids.tsv"))
    fwrite(dup, dup_file, sep = "\t")
    stop("Duplicated sample_clean IDs in OOF file. See: ", dup_file)
  }

  oof[, true_label := as.character(true_label)]
  oof[, predicted_label := as.character(predicted_label)]
  oof[!is.na(true_label) & !is.na(predicted_label) & true_label == predicted_label, sample_clean]
}

apply_cc_filter <- function(train_dt) {
  if (!require_cc_ge50) {
    return(list(active = FALSE, keep_samples = train_dt$sample_clean, tag_suffix = "sameSource_allTraining"))
  }

  if (!nzchar(cc_file) || !file.exists(cc_file)) {
    stop("require_cc_ge50=TRUE but cc_file is missing or does not exist: ", cc_file)
  }

  cc_dt <- fread(cc_file)
  if (!all(c(cc_sample_col, cc_col) %in% names(cc_dt))) {
    stop("CC file must contain columns: ", cc_sample_col, ", ", cc_col)
  }

  setnames(cc_dt, c(cc_sample_col, cc_col), c("sample", "CC"))
  cc_dt[, sample := as.character(sample)]
  cc_dt[, sample_clean := clean_id(sample)]
  cc_dt[, CC := as.character(CC)]
  cc_dt <- cc_dt[!is.na(CC) & CC != "" & CC != "ND"]

  cc_counts <- cc_dt[, .N, by = CC]
  cc_keep <- cc_counts[N >= cc_min_count, CC]
  keep <- unique(cc_dt[CC %in% cc_keep, sample_clean])

  logf("CC filter active: retained %d samples in CC >= %d", length(keep), cc_min_count)
  list(active = TRUE, keep_samples = keep, tag_suffix = paste0("sameSource_CCge", cc_min_count))
}

calc_one_query <- function(q_idx, X, sample_ids, candidate_idx, k_top, missing_code) {
  q <- X[q_idx, ]
  cand <- candidate_idx[candidate_idx != q_idx]

  if (length(cand) < k_top) {
    return(data.table(
      sample = sample_ids[q_idx],
      status = "excluded_less_than_11_same_source_candidates",
      n_candidates = length(cand),
      n_valid_distance_candidates = 0L
    ))
  }

  cand_x <- X[cand, , drop = FALSE]
  q_valid <- q != missing_code
  valid <- sweep(cand_x != missing_code, 2, q_valid, `&`)
  valid_counts <- rowSums(valid)

  mismatch <- rowSums(valid & sweep(cand_x, 2, q, `!=`))
  ok <- valid_counts > 0

  if (sum(ok) < k_top) {
    return(data.table(
      sample = sample_ids[q_idx],
      status = "excluded_less_than_11_valid_distance_candidates",
      n_candidates = length(cand),
      n_valid_distance_candidates = sum(ok)
    ))
  }

  cand_ok <- cand[ok]
  mismatch_ok <- mismatch[ok]
  valid_ok <- valid_counts[ok]
  normalized_ok <- mismatch_ok / valid_ok

  raw_order <- order(mismatch_ok, normalized_ok, sample_ids[cand_ok], na.last = NA)
  raw_top <- raw_order[seq_len(k_top)]

  norm_order <- order(normalized_ok, mismatch_ok, sample_ids[cand_ok], na.last = NA)
  norm_top <- norm_order[seq_len(k_top)]

  data.table(
    sample = sample_ids[q_idx],
    status = "OK",
    n_candidates = length(cand),
    n_valid_distance_candidates = sum(ok),

    raw_d1_dismatch_count = as.numeric(mismatch_ok[raw_top][1]),
    raw_d11_dismatch_count = as.numeric(mismatch_ok[raw_top][k_top]),
    raw_mean_top11_dismatch_count = mean(mismatch_ok[raw_top]),
    raw_median_top11_dismatch_count = median(mismatch_ok[raw_top]),
    raw_max_top11_dismatch_count = max(mismatch_ok[raw_top]),
    raw_mean_top11_valid_loci = mean(valid_ok[raw_top]),
    raw_nearest_train_samples = paste(sample_ids[cand_ok[raw_top]], collapse = ";"),
    raw_nearest_dismatch_counts = paste(mismatch_ok[raw_top], collapse = ";"),
    raw_nearest_normalized_distances = paste(signif(normalized_ok[raw_top], 10), collapse = ";"),
    raw_nearest_valid_loci = paste(valid_ok[raw_top], collapse = ";"),

    norm_d1_normalized = as.numeric(normalized_ok[norm_top][1]),
    norm_d11_normalized = as.numeric(normalized_ok[norm_top][k_top]),
    norm_mean_top11_normalized = mean(normalized_ok[norm_top]),
    norm_median_top11_normalized = median(normalized_ok[norm_top]),
    norm_max_top11_normalized = max(normalized_ok[norm_top]),
    norm_d1_dismatch_count = as.numeric(mismatch_ok[norm_top][1]),
    norm_d11_dismatch_count = as.numeric(mismatch_ok[norm_top][k_top]),
    norm_mean_top11_dismatch_count = mean(mismatch_ok[norm_top]),
    norm_mean_top11_valid_loci = mean(valid_ok[norm_top]),
    norm_nearest_train_samples = paste(sample_ids[cand_ok[norm_top]], collapse = ";"),
    norm_nearest_normalized_distances = paste(signif(normalized_ok[norm_top], 10), collapse = ";"),
    norm_nearest_dismatch_counts = paste(mismatch_ok[norm_top], collapse = ";"),
    norm_nearest_valid_loci = paste(valid_ok[norm_top], collapse = ";")
  )
}

make_chunks <- function(x, chunk_size) {
  split(x, ceiling(seq_along(x) / chunk_size))
}

compute_same_source_distances <- function(train_dt, X, eligible_samples) {
  sample_ids <- train_dt$sample_clean
  labels <- as.character(train_dt$label)
  eligible_set <- unique(eligible_samples)
  query_indices_all <- which(sample_ids %in% eligible_set)

  out <- list()
  idx <- 1L

  for (src in sort(unique(labels))) {
    source_candidate_idx <- which(labels == src)
    source_query_idx <- query_indices_all[labels[query_indices_all] == src]

    logf("Source %s | candidates=%d | query_samples=%d", src, length(source_candidate_idx), length(source_query_idx))

    if (length(source_query_idx) == 0) next

    chunks <- make_chunks(source_query_idx, chunk_size)
    res_list <- mclapply(
      chunks,
      function(ch) {
        rbindlist(lapply(
          ch,
          calc_one_query,
          X = X,
          sample_ids = sample_ids,
          candidate_idx = source_candidate_idx,
          k_top = k_top,
          missing_code = missing_code
        ), use.names = TRUE, fill = TRUE)
      },
      mc.cores = n_cores
    )

    src_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
    src_dt[, source := src]
    out[[idx]] <- src_dt
    idx <- idx + 1L
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

make_thresholds_by_source <- function(dt, value_col, distance_unit, metric_name) {
  ok_dt <- dt[status == "OK" & !is.na(get(value_col))]
  ans <- ok_dt[, calc_stats(get(value_col)), by = source]
  ans[, `:=`(
    dataset = dataset_label,
    distance_unit = distance_unit,
    metric = metric_name,
    value_column = value_col,
    k_top = k_top
  )]
  setcolorder(ans, c("dataset", "source", "distance_unit", "metric", "value_column", "k_top"))
  setorder(ans, source)
  ans
}

write_threshold_table <- function(thr, out_dir, filename) {
  fwrite(thr, file.path(out_dir, filename), sep = "\t", quote = FALSE)
}

# ============================================================
# Main
# ============================================================
start_time <- Sys.time()

for (ff in c(train_file, model_file)) {
  if (!file.exists(ff)) stop("Missing input file: ", ff)
}

logf("Dataset: %s", dataset_label)
logf("Training matrix: %s", train_file)
logf("Model file: %s", model_file)
logf("Raw output dir: %s", raw_out_dir)
logf("Normalized output dir: %s", normalized_out_dir)

model <- readRDS(model_file)
if (!all(c("feature_order", "label_levels") %in% names(model))) {
  stop("Model RDS must contain feature_order and label_levels.")
}
feature_order <- as.character(model$feature_order)
label_levels <- as.character(model$label_levels)
if (length(feature_order) != expected_features) {
  logf("WARNING: expected %d features, but model contains %d features", expected_features, length(feature_order))
}

train_dt <- fread(train_file, sep = "\t", data.table = TRUE, showProgress = TRUE)
if (ncol(train_dt) < 3) stop("Training matrix must contain sample, label, and feature columns.")
setnames(train_dt, names(train_dt)[1:2], c("sample", "label"))
train_dt[, sample := as.character(sample)]
train_dt[, sample_clean := clean_id(sample)]
train_dt[, label := as.character(label)]

if (anyDuplicated(train_dt$sample_clean)) {
  dup <- train_dt[duplicated(sample_clean) | duplicated(sample_clean, fromLast = TRUE), .(sample, sample_clean, label)]
  dup_file <- file.path(raw_out_dir, paste0(tag, "_duplicated_training_sample_clean_ids.tsv"))
  fwrite(dup, dup_file, sep = "\t")
  stop("Duplicated sample_clean IDs in training data. See: ", dup_file)
}

missing_features <- setdiff(feature_order, names(train_dt))
if (length(missing_features) > 0) {
  miss_file <- file.path(raw_out_dir, paste0(tag, "_missing_model_features_in_training.log"))
  writeLines(missing_features, miss_file)
  stop("Training matrix is missing model features. See: ", miss_file)
}

train_dt <- train_dt[, c("sample", "sample_clean", "label", feature_order), with = FALSE]
for (cc in feature_order) {
  train_dt[[cc]] <- clean_allele_vec(train_dt[[cc]], missing_code)
}
X <- as.matrix(train_dt[, ..feature_order])
storage.mode(X) <- "integer"

label_count <- train_dt[, .N, by = label][order(label)]
fwrite(label_count, file.path(raw_out_dir, paste0(tag, "_training_label_counts.tsv")), sep = "\t")
fwrite(label_count, file.path(normalized_out_dir, paste0(tag, "_training_label_counts.tsv")), sep = "\t")

eligible_samples <- train_dt$sample_clean
filter_notes <- c("all_training_samples")

cc_res <- apply_cc_filter(train_dt)
eligible_samples <- intersect(eligible_samples, cc_res$keep_samples)
analysis_tag <- cc_res$tag_suffix
if (cc_res$active) filter_notes <- c(filter_notes, paste0("CC_ge_", cc_min_count))

oof_file <- find_existing(oof_file_candidates)
if (use_correct_oof_only && !is.na(oof_file)) {
  correct_samples <- read_oof_correct_samples(oof_file, train_dt, label_levels)
  eligible_samples <- intersect(eligible_samples, correct_samples)
  filter_notes <- c(filter_notes, "OOF_correct_only")
  logf("OOF correct filter active: retained %d samples", length(correct_samples))
} else if (use_correct_oof_only && is.na(oof_file)) {
  logf("WARNING: use_correct_oof_only=TRUE but no OOF file was found. Thresholds will use all eligible training samples.")
  filter_notes <- c(filter_notes, "OOF_file_not_found")
} else {
  filter_notes <- c(filter_notes, "OOF_filter_not_used")
}

if (length(eligible_samples) == 0) {
  stop("No eligible samples remain for threshold calculation.")
}

eligible_dt <- train_dt[sample_clean %in% eligible_samples, .N, by = label][order(label)]
fwrite(eligible_dt, file.path(raw_out_dir, paste0(tag, "_eligible_query_samples_by_source.tsv")), sep = "\t")
fwrite(eligible_dt, file.path(normalized_out_dir, paste0(tag, "_eligible_query_samples_by_source.tsv")), sep = "\t")

logf("Eligible query samples: %d", length(eligible_samples))

same_source_dt <- compute_same_source_distances(train_dt, X, eligible_samples)
same_source_dt <- merge(
  same_source_dt,
  train_dt[, .(sample = sample_clean, original_sample = sample, true_label = label)],
  by = "sample",
  all.x = TRUE
)
setcolorder(same_source_dt, c("sample", "original_sample", "true_label", "source"))

raw_dist_file <- file.path(raw_out_dir, paste0(tag, "_", analysis_tag, "_training_sameSource_raw_and_normalized_distances.tsv"))
norm_dist_file <- file.path(normalized_out_dir, paste0(tag, "_", analysis_tag, "_training_sameSource_raw_and_normalized_distances.tsv"))
fwrite(same_source_dt, raw_dist_file, sep = "\t", quote = FALSE)
fwrite(same_source_dt, norm_dist_file, sep = "\t", quote = FALSE)

status_summary <- same_source_dt[, .N, by = status][order(status)]
fwrite(status_summary, file.path(raw_out_dir, paste0(tag, "_", analysis_tag, "_status_summary.tsv")), sep = "\t")
fwrite(status_summary, file.path(normalized_out_dir, paste0(tag, "_", analysis_tag, "_status_summary.tsv")), sep = "\t")

raw_d1_thr <- make_thresholds_by_source(same_source_dt, "raw_d1_dismatch_count", "raw_dismatch_count", "d1")
raw_mean_thr <- make_thresholds_by_source(same_source_dt, "raw_mean_top11_dismatch_count", "raw_dismatch_count", "meanTop11")
raw_d11_thr <- make_thresholds_by_source(same_source_dt, "raw_d11_dismatch_count", "raw_dismatch_count", "d11")

norm_d1_thr <- make_thresholds_by_source(same_source_dt, "norm_d1_normalized", "normalized_hamming", "d1")
norm_mean_thr <- make_thresholds_by_source(same_source_dt, "norm_mean_top11_normalized", "normalized_hamming", "meanTop11")
norm_d11_thr <- make_thresholds_by_source(same_source_dt, "norm_d11_normalized", "normalized_hamming", "d11")

write_threshold_table(raw_d1_thr, raw_out_dir, paste0(tag, "_", analysis_tag, "_d1_dismatch_thresholds_by_source.tsv"))
write_threshold_table(raw_mean_thr, raw_out_dir, paste0(tag, "_", analysis_tag, "_meanTop11_dismatch_thresholds_by_source.tsv"))
write_threshold_table(raw_d11_thr, raw_out_dir, paste0(tag, "_", analysis_tag, "_d11_dismatch_thresholds_by_source.tsv"))

write_threshold_table(norm_d1_thr, normalized_out_dir, paste0(tag, "_", analysis_tag, "_d1_normalized_thresholds_by_source.tsv"))
write_threshold_table(norm_mean_thr, normalized_out_dir, paste0(tag, "_", analysis_tag, "_meanTop11_normalized_thresholds_by_source.tsv"))
write_threshold_table(norm_d11_thr, normalized_out_dir, paste0(tag, "_", analysis_tag, "_d11_normalized_thresholds_by_source.tsv"))

all_thr <- rbindlist(
  list(raw_d1_thr, raw_mean_thr, raw_d11_thr, norm_d1_thr, norm_mean_thr, norm_d11_thr),
  use.names = TRUE,
  fill = TRUE
)
fwrite(all_thr, file.path(raw_out_dir, paste0(tag, "_", analysis_tag, "_all_thresholds_by_source_long.tsv")), sep = "\t", quote = FALSE)
fwrite(all_thr, file.path(normalized_out_dir, paste0(tag, "_", analysis_tag, "_all_thresholds_by_source_long.tsv")), sep = "\t", quote = FALSE)

if (write_legacy_ccge50_filename && !cc_res$active) {
  legacy_file <- file.path(normalized_out_dir, paste0(tag, "_sameSource_CCge50_meanTop11_normalized_thresholds_by_source.tsv"))
  fwrite(norm_mean_thr, legacy_file, sep = "\t", quote = FALSE)
  logf("Legacy compatibility file written: %s", legacy_file)
}

run_summary <- data.table(
  item = c(
    "dataset",
    "tag",
    "train_file",
    "model_file",
    "oof_file_used",
    "use_correct_oof_only",
    "require_cc_ge50",
    "cc_file",
    "cc_min_count",
    "analysis_tag",
    "filter_notes",
    "n_training_samples",
    "n_eligible_query_samples",
    "n_features",
    "k_top",
    "distance_rule_raw",
    "distance_rule_normalized",
    "raw_out_dir",
    "normalized_out_dir",
    "elapsed_seconds"
  ),
  value = c(
    dataset_label,
    tag,
    train_file,
    model_file,
    ifelse(is.na(oof_file), "", oof_file),
    as.character(use_correct_oof_only),
    as.character(require_cc_ge50),
    cc_file,
    as.character(cc_min_count),
    analysis_tag,
    paste(filter_notes, collapse = ";"),
    as.character(nrow(train_dt)),
    as.character(length(eligible_samples)),
    as.character(length(feature_order)),
    as.character(k_top),
    "dismatch_count among shared non-missing loci; nearest neighbors sorted by raw dismatch_count",
    "dismatch_count / valid_shared_loci; nearest neighbors sorted by normalized_hamming",
    raw_out_dir,
    normalized_out_dir,
    as.character(round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
  )
)

fwrite(run_summary, file.path(raw_out_dir, paste0(tag, "_", analysis_tag, "_threshold_run_summary.tsv")), sep = "\t", quote = FALSE)
fwrite(run_summary, file.path(normalized_out_dir, paste0(tag, "_", analysis_tag, "_threshold_run_summary.tsv")), sep = "\t", quote = FALSE)

cat("\n===== Raw meanTop11 dismatch thresholds =====\n")
print(raw_mean_thr)
cat("\n===== Normalized meanTop11 thresholds =====\n")
print(norm_mean_thr)
cat("\n[DONE] Outputs saved to:\n")
cat(raw_out_dir, "\n")
cat(normalized_out_dir, "\n")
