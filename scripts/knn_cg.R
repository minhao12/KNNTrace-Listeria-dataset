#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_pkgs <- c("data.table", "caret", "ggplot2", "foreach", "doParallel", "progress", "parallel")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  library(data.table)
  library(caret)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(progress)
  library(parallel)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

set.seed(42)
setDTthreads(as.integer(Sys.getenv("DATA_TABLE_THREADS", "30")))

input_file <- "/home/minhao/xiezuo/data/cgmlst/ml_ready_dataset.tsv"
output_dir <- "/home/minhao/knn/xinxinxinknn/moxing/knn/cgmlst1748"
dataset_name <- "cgmlst1748"
expected_features <- 1748L
k_neighbors <- 11L
weight_mode <- "inv1"
n_cores <- as.integer(Sys.getenv("N_CORES", "30"))
missing_code <- -1L
confidence_eps <- 1e-6

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n", sep = "")
  flush.console()
}

clean_allele_vector <- function(x) {
  v <- trimws(as.character(x))
  v <- sub("^INF[-_]", "", v)
  invalid <- is.na(v) | v == "" | v == "-" |
    v %in% c("NA", "N/A", "NaN", "nan", "NULL", "null",
             "LNF", "NIPH", "NIPHEM", "ALM", "ASM",
             "PLOT3", "PLOT5", "LOTSC", "PAMA") |
    grepl("^[A-Za-z]", v)
  v[invalid] <- NA_character_
  z <- suppressWarnings(as.integer(as.numeric(v)))
  z[is.na(z)] <- missing_code
  as.integer(z)
}

clean_allele_matrix <- function(dt, feature_cols) {
  for (col in feature_cols) {
    dt[[col]] <- clean_allele_vector(dt[[col]])
  }
  invisible(dt)
}

hamming_ignore_missing <- function(a, b) {
  ok <- (a != -1L) & (b != -1L)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(1)
  1 - sum(a[ok] == b[ok]) / n_ok
}

get_weights <- function(nn_dists, mode = "inv1", eps = 1e-6) {
  if (mode == "inv1") {
    return(1 / (nn_dists + eps))
  }
  if (mode == "inv2") {
    return(1 / (nn_dists + eps)^2)
  }
  stop("Unsupported weight mode: ", mode)
}

calc_macro_f1 <- function(pred_label, true_label, label_levels) {
  tab <- table(
    factor(pred_label, levels = label_levels),
    factor(true_label, levels = label_levels)
  )
  f1_values <- numeric(length(label_levels))
  for (i in seq_along(label_levels)) {
    cls <- label_levels[i]
    tp <- tab[cls, cls]
    fp <- sum(tab[cls, ]) - tp
    fn <- sum(tab[, cls]) - tp
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    f1_values[i] <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  }
  mean(f1_values)
}

make_or_read_folds <- function(dt, fold_file, k = 5L) {
  if (file.exists(fold_file)) {
    fold_dt <- fread(fold_file)
    required_cols <- c("sample", "label", "fold")
    if (!all(required_cols %in% names(fold_dt))) {
      stop("Existing fold file does not contain required columns: ", fold_file)
    }
    idx <- match(dt$sample, fold_dt$sample)
    if (any(is.na(idx))) {
      stop("Existing fold file is missing samples from the input matrix: ", fold_file)
    }
    fold_dt <- fold_dt[idx]
  } else {
    folds <- caret::createFolds(dt$label, k = k, list = TRUE, returnTrain = FALSE)
    fold_dt <- rbindlist(lapply(seq_along(folds), function(i) {
      data.table(sample = dt$sample[folds[[i]]], label = as.character(dt$label[folds[[i]]]), fold = i)
    }))
    fwrite(fold_dt, fold_file, sep = "\t", quote = FALSE)
  }
  split(seq_len(nrow(dt)), fold_dt$fold)
}

run_knn_cv <- function(dt, fold_file) {
  feature_cols <- setdiff(names(dt), c("sample", "label"))
  label_levels <- levels(dt$label)
  folds <- make_or_read_folds(dt, fold_file, k = 5L)

  clean_allele_matrix(dt, feature_cols)
  x_all <- as.matrix(dt[, ..feature_cols])
  storage.mode(x_all) <- "integer"
  y_all <- dt$label

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)

  parallel::clusterExport(
    cl,
    c("hamming_ignore_missing", "get_weights", "k_neighbors", "weight_mode", "confidence_eps"),
    envir = environment()
  )

  metric_rows <- vector("list", length(folds))
  class_rows <- vector("list", length(folds))
  oof_rows <- vector("list", length(folds))

  pb <- progress::progress_bar$new(
    format = "[:bar] :current/:total :percent fold=:fold",
    total = length(folds), width = 60, clear = FALSE
  )

  for (fold_id in seq_along(folds)) {
    pb$tick(tokens = list(fold = fold_id))

    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(dt)), test_idx)

    x_train <- x_all[train_idx, , drop = FALSE]
    y_train <- y_all[train_idx]
    x_test <- x_all[test_idx, , drop = FALSE]
    y_test <- y_all[test_idx]

    fold_res <- foreach::foreach(
      j = seq_len(nrow(x_test)),
      .combine = rbind,
      .export = c("label_levels", "k_neighbors", "weight_mode", "confidence_eps",
                  "hamming_ignore_missing", "get_weights")
    ) %dopar% {
      x <- x_test[j, ]
      dist_vec <- apply(x_train, 1L, hamming_ignore_missing, b = x)
      k_eff <- min(k_neighbors, length(dist_vec))
      nn_idx <- order(dist_vec, decreasing = FALSE)[seq_len(k_eff)]
      nn_labels <- y_train[nn_idx]
      nn_dists <- dist_vec[nn_idx]

      weights <- get_weights(nn_dists, mode = weight_mode, eps = confidence_eps)
      weights <- weights / sum(weights)
      prob_vec <- setNames(rep(0, length(label_levels)), label_levels)
      for (cls in unique(nn_labels)) {
        prob_vec[as.character(cls)] <- sum(weights[nn_labels == cls])
      }
      predicted <- names(which.max(prob_vec))
      c(predicted_label = predicted, prob_vec)
    }

    fold_res <- as.data.table(fold_res)
    prob_cols <- setdiff(names(fold_res), "predicted_label")
    for (col in prob_cols) fold_res[[col]] <- as.numeric(fold_res[[col]])

    pred_factor <- factor(fold_res$predicted_label, levels = label_levels)
    cm <- caret::confusionMatrix(pred_factor, y_test)

    metric_rows[[fold_id]] <- data.table(
      dataset = dataset_name,
      feature_set = "cgMLST-1748",
      n_features = length(feature_cols),
      model = "KNN_Hamming",
      k = k_neighbors,
      weight_mode = weight_mode,
      fold = fold_id,
      Accuracy = as.numeric(cm$overall["Accuracy"]),
      Kappa = as.numeric(cm$overall["Kappa"]),
      MacroF1 = calc_macro_f1(pred_factor, y_test, label_levels)
    )

    by_class <- as.data.table(cm$byClass, keep.rownames = TRUE)
    setnames(by_class, "rn", "Class")
    by_class[, Class := sub("^Class: ", "", Class)]
    if (!"F1" %in% names(by_class) && all(c("Sensitivity", "Pos Pred Value") %in% names(by_class))) {
      by_class[, F1 := 2 * (`Pos Pred Value` * Sensitivity) / pmax(`Pos Pred Value` + Sensitivity, 1e-9)]
    }

    class_rows[[fold_id]] <- by_class[, .(
      dataset = dataset_name,
      feature_set = "cgMLST-1748",
      n_features = length(feature_cols),
      model = "KNN_Hamming",
      k = k_neighbors,
      weight_mode = weight_mode,
      fold = fold_id,
      Class,
      Precision = `Pos Pred Value`,
      Recall = Sensitivity,
      F1
    )]

    oof_dt <- data.table(
      sample = dt$sample[test_idx],
      true_label = as.character(y_test),
      predicted_label = as.character(pred_factor),
      fold = fold_id,
      dataset = dataset_name,
      feature_set = "cgMLST-1748",
      n_features = length(feature_cols),
      model = "KNN_Hamming",
      k = k_neighbors,
      weight_mode = weight_mode
    )
    oof_dt <- cbind(oof_dt, fold_res[, ..label_levels])
    oof_rows[[fold_id]] <- oof_dt

    log_msg("Fold %d | Accuracy=%.4f | Kappa=%.4f | MacroF1=%.4f",
            fold_id, metric_rows[[fold_id]]$Accuracy, metric_rows[[fold_id]]$Kappa, metric_rows[[fold_id]]$MacroF1)
  }

  list(
    metrics_by_fold = rbindlist(metric_rows, fill = TRUE),
    classwise_by_fold = rbindlist(class_rows, fill = TRUE),
    oof_predictions = rbindlist(oof_rows, fill = TRUE),
    X_all = x_all,
    y_all = y_all,
    sample_ids = dt$sample,
    feature_order = feature_cols,
    label_levels = label_levels
  )
}

summarize_results <- function(res) {
  metrics_summary <- res$metrics_by_fold[, .(
    Mean_Accuracy = mean(Accuracy),
    SD_Accuracy = sd(Accuracy),
    Mean_Kappa = mean(Kappa),
    SD_Kappa = sd(Kappa),
    Mean_MacroF1 = mean(MacroF1),
    SD_MacroF1 = sd(MacroF1),
    n_folds = .N
  ), by = .(dataset, feature_set, n_features, model, k, weight_mode)]

  classwise_summary <- res$classwise_by_fold[, .(
    Mean_Precision = mean(Precision, na.rm = TRUE),
    SD_Precision = sd(Precision, na.rm = TRUE),
    Mean_Recall = mean(Recall, na.rm = TRUE),
    SD_Recall = sd(Recall, na.rm = TRUE),
    Mean_F1 = mean(F1, na.rm = TRUE),
    SD_F1 = sd(F1, na.rm = TRUE),
    n_folds = .N
  ), by = .(dataset, feature_set, n_features, model, k, weight_mode, Class)][order(Class)]

  list(metrics_summary = metrics_summary, classwise_summary = classwise_summary)
}

save_model_bundle <- function(res, model_file) {
  label_distribution <- data.table(label = as.character(res$y_all))[, .N, by = label][order(label)]
  model_bundle <- list(
    model_type = "KNN_Hamming",
    dataset = dataset_name,
    feature_set = "cgMLST-1748",
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    input_file = input_file,
    n_samples = nrow(res$X_all),
    n_features = ncol(res$X_all),
    k_neighbors = k_neighbors,
    weight_mode = weight_mode,
    missing_code = missing_code,
    distance_rule = "Normalized Hamming distance calculated over loci with valid allele calls in both isolates; missing value -1 is ignored.",
    vote_rule = "Inverse-distance weighted voting with weight = 1 / (distance + 1e-6).",
    sample_ids = res$sample_ids,
    feature_order = res$feature_order,
    label_levels = res$label_levels,
    label_distribution = label_distribution,
    X_train = res$X_all,
    y_train = res$y_all
  )
  saveRDS(model_bundle, model_file, compress = "gzip")
  invisible(model_bundle)
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

log_msg("Reading input: %s", input_file)
dt <- fread(input_file, sep = "\t", data.table = TRUE, showProgress = TRUE)
if (ncol(dt) < 3) stop("Input matrix must contain at least three columns.")
setnames(dt, names(dt)[1:2], c("sample", "label"))
dt[, sample := as.character(sample)]
dt[, label := factor(label)]
feature_cols <- setdiff(names(dt), c("sample", "label"))
if (length(feature_cols) != expected_features) {
  stop(sprintf("Expected %d features, but detected %d features.", expected_features, length(feature_cols)))
}

label_counts <- dt[, .N, by = label][order(label)]
fwrite(label_counts, file.path(output_dir, "label_counts_cgmlst1748.tsv"), sep = "\t", quote = FALSE)

fold_file <- file.path(output_dir, "fold_assignment_cgmlst1748.tsv")
res <- run_knn_cv(dt, fold_file)
sum_res <- summarize_results(res)

fwrite(res$metrics_by_fold, file.path(output_dir, "metrics_by_fold_cgmlst1748.tsv"), sep = "\t", quote = FALSE)
fwrite(sum_res$metrics_summary, file.path(output_dir, "metrics_summary_cgmlst1748.tsv"), sep = "\t", quote = FALSE)
fwrite(res$classwise_by_fold, file.path(output_dir, "classwise_by_fold_cgmlst1748.tsv"), sep = "\t", quote = FALSE)
fwrite(sum_res$classwise_summary, file.path(output_dir, "classwise_summary_cgmlst1748.tsv"), sep = "\t", quote = FALSE)
fwrite(res$oof_predictions, file.path(output_dir, "oof_predictions_cgmlst1748.tsv"), sep = "\t", quote = FALSE)
writeLines(res$feature_order, file.path(output_dir, "feature_order_cgmlst1748.txt"))
writeLines(res$label_levels, file.path(output_dir, "label_levels_cgmlst1748.txt"))

model_file <- file.path(output_dir, "knn_cgmlst1748_k11_inv1_model.rds")
model_bundle <- save_model_bundle(res, model_file)

model_summary <- data.table(
  item = c("model_type", "dataset", "feature_set", "created_at", "input_file", "n_samples", "n_features", "k_neighbors", "weight_mode", "missing_code", "model_file"),
  value = c(model_bundle$model_type, model_bundle$dataset, model_bundle$feature_set, model_bundle$created_at,
            model_bundle$input_file, as.character(model_bundle$n_samples), as.character(model_bundle$n_features),
            as.character(model_bundle$k_neighbors), model_bundle$weight_mode, as.character(model_bundle$missing_code), model_file)
)
fwrite(model_summary, file.path(output_dir, "model_summary_cgmlst1748.tsv"), sep = "\t", quote = FALSE)

p_metrics <- melt(
  sum_res$metrics_summary,
  id.vars = c("dataset", "feature_set", "n_features", "model", "k", "weight_mode", "n_folds"),
  measure.vars = c("Mean_Accuracy", "Mean_Kappa", "Mean_MacroF1"),
  variable.name = "Metric",
  value.name = "Value"
)
p <- ggplot(p_metrics, aes(x = Metric, y = Value)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.4f", Value)), vjust = -0.35, size = 4) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "KNN-Hamming cgMLST-1748 5-fold performance", x = NULL, y = "Mean value") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5))
ggsave(file.path(output_dir, "knn_cgmlst1748_5fold_metrics.png"), p, width = 7, height = 5, dpi = 300)

writeLines(c(
  "KNN-Hamming cgMLST-1748 5-fold CV",
  sprintf("Created at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Input file: %s", input_file),
  sprintf("Output directory: %s", output_dir),
  sprintf("Samples: %d", nrow(res$X_all)),
  sprintf("Features: %d", ncol(res$X_all)),
  sprintf("K: %d", k_neighbors),
  sprintf("Weight mode: %s", weight_mode),
  sprintf("Model file: %s", model_file)
), file.path(output_dir, "run_summary_cgmlst1748.txt"))

log_msg("Done. Results written to: %s", output_dir)
log_msg("Model saved to: %s", model_file)
