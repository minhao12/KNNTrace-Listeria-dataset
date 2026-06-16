#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
})

setDTthreads(30)

DATASET_NAME <- "xgb_cgmlst1748"
FEATURE_REPRESENTATION <- "cgMLST-1748"
INPUT_FILE <- "/home/minhao/xiezuo/data/cgmlst/ml_ready_dataset.tsv"
OUTPUT_DIR <- "/home/minhao/knn/xinxinxinknn/moxing/xgb/cgmlst1748"
SHARED_FOLD_FILE <- "/home/minhao/knn/xinxinxinknn/daima/shared_5fold_assignment.tsv"
N_FOLDS <- 5L
SEED <- 20250609L
MISSING_CODE <- -1

XGB_PARAMS <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  max_depth = 7,
  eta = 0.08,
  min_child_weight = 2,
  subsample = 0.8,
  colsample_bytree = 0.8,
  nthread = 10
)
NROUNDS <- 2000L

OUTPUT_PREFIX <- file.path(OUTPUT_DIR, DATASET_NAME)
FOLD_OUT_FILE <- paste0(OUTPUT_PREFIX, "_fold_assignment.tsv")
LABEL_COUNTS_FILE <- paste0(OUTPUT_PREFIX, "_label_counts.tsv")
OOF_FILE <- paste0(OUTPUT_PREFIX, "_oof_predictions.tsv")
METRICS_BY_FOLD_FILE <- paste0(OUTPUT_PREFIX, "_metrics_by_fold.tsv")
METRICS_SUMMARY_FILE <- paste0(OUTPUT_PREFIX, "_metrics_summary.tsv")
CLASSWISE_BY_FOLD_FILE <- paste0(OUTPUT_PREFIX, "_classwise_by_fold.tsv")
CLASSWISE_SUMMARY_FILE <- paste0(OUTPUT_PREFIX, "_classwise_summary.tsv")
PARAMS_FILE <- paste0(OUTPUT_PREFIX, "_parameters.tsv")
MODEL_RDS_FILE <- paste0(OUTPUT_PREFIX, "_model.rds")
MODEL_BINARY_FILE <- paste0(OUTPUT_PREFIX, "_model.xgb")
SUMMARY_FILE <- paste0(OUTPUT_PREFIX, "_run_summary.txt")
METRICS_PLOT_FILE <- paste0(OUTPUT_PREFIX, "_5fold_metrics.png")
RECALL_PLOT_FILE <- paste0(OUTPUT_PREFIX, "_recall_by_class.png")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
}

check_required_packages <- function() {
  required <- c("data.table", "xgboost")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
}

clean_feature_vector <- function(x) {
  z <- as.character(x)
  z <- trimws(z)
  z <- gsub("^INF[-_]", "", z, ignore.case = TRUE)
  invalid <- is.na(z) | z == "" | z %in% c("-", ".", "NA", "N/A", "NaN", "NULL")
  invalid <- invalid | grepl("^(LNF|NIPH|NIPHEM|ALM|ASM|PLOT3|PLOT5|LOTSC|PAMA)$", z, ignore.case = TRUE)
  suppressWarnings(num <- as.numeric(z))
  num[invalid | is.na(num)] <- MISSING_CODE
  num
}

read_input_matrix <- function(input_file) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  dt <- fread(input_file)
  if (ncol(dt) < 3L) {
    stop("Input file must contain at least sample, label, and one feature column.")
  }
  setnames(dt, names(dt)[1:2], c("sample", "label"))
  dt[, sample := as.character(sample)]
  dt[, label := as.character(label)]
  feature_cols <- setdiff(names(dt), c("sample", "label"))
  if (length(feature_cols) < 1L) {
    stop("No feature columns were found after sample and label columns.")
  }
  log_msg("Cleaning %d feature columns.", length(feature_cols))
  dt[, (feature_cols) := lapply(.SD, clean_feature_vector), .SDcols = feature_cols]
  list(dt = dt, feature_cols = feature_cols)
}

create_stratified_folds <- function(samples, labels, n_folds, seed) {
  set.seed(seed)
  fold_dt <- data.table(sample = as.character(samples), label = as.character(labels), fold = NA_integer_)
  for (lv in sort(unique(fold_dt$label))) {
    idx <- which(fold_dt$label == lv)
    idx <- sample(idx, length(idx))
    fold_dt$fold[idx] <- rep(seq_len(n_folds), length.out = length(idx))
  }
  fold_dt[order(sample)]
}

get_or_create_folds <- function(samples, labels, shared_fold_file, fold_out_file) {
  dir.create(dirname(shared_fold_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(fold_out_file), recursive = TRUE, showWarnings = FALSE)

  samples <- as.character(samples)
  labels <- as.character(labels)

  if (file.exists(shared_fold_file)) {
    fold_dt <- fread(shared_fold_file)
    if (!all(c("sample", "fold") %in% names(fold_dt))) {
      stop("Shared fold file must contain columns: sample and fold.")
    }
    fold_dt[, sample := as.character(sample)]
    missing_samples <- setdiff(samples, fold_dt$sample)
    if (length(missing_samples) > 0L) {
      stop("Shared fold file does not contain all input samples. Missing examples: ", paste(head(missing_samples, 10), collapse = ", "))
    }
    fold_dt <- fold_dt[match(samples, sample)]
    if (any(is.na(fold_dt$fold))) {
      stop("Invalid fold assignment detected for matched samples.")
    }
    fold_dt[, label := labels]
    fold_dt <- fold_dt[, .(sample, label, fold)]
    fwrite(fold_dt, fold_out_file, sep = "\t")
    return(fold_dt)
  }

  fold_dt <- create_stratified_folds(samples, labels, N_FOLDS, SEED)
  fold_dt <- fold_dt[match(samples, sample)]
  fwrite(fold_dt, shared_fold_file, sep = "\t")
  fwrite(fold_dt, fold_out_file, sep = "\t")
  fold_dt
}

safe_div <- function(num, den) {
  out <- ifelse(den == 0, 0, num / den)
  out[is.na(out) | is.nan(out) | is.infinite(out)] <- 0
  out
}

calculate_metrics <- function(true_labels, pred_labels, label_levels) {
  true_factor <- factor(true_labels, levels = label_levels)
  pred_factor <- factor(pred_labels, levels = label_levels)
  cm <- table(true_factor, pred_factor)
  n <- sum(cm)
  accuracy <- ifelse(n == 0, NA_real_, sum(diag(cm)) / n)
  row_marg <- rowSums(cm)
  col_marg <- colSums(cm)
  expected <- ifelse(n == 0, NA_real_, sum(row_marg * col_marg) / (n * n))
  kappa <- ifelse(is.na(expected) || abs(1 - expected) < 1e-12, NA_real_, (accuracy - expected) / (1 - expected))

  tp <- diag(cm)
  precision <- safe_div(tp, colSums(cm))
  recall <- safe_div(tp, rowSums(cm))
  f1 <- safe_div(2 * precision * recall, precision + recall)

  classwise <- data.table(
    Class = label_levels,
    Precision = as.numeric(precision),
    Recall = as.numeric(recall),
    F1 = as.numeric(f1),
    Support = as.numeric(rowSums(cm))
  )

  list(
    overall = data.table(
      Accuracy = accuracy,
      Kappa = kappa,
      MacroF1 = mean(classwise$F1, na.rm = TRUE)
    ),
    classwise = classwise
  )
}

make_dmatrix <- function(x, y = NULL) {
  if (is.null(y)) {
    xgb.DMatrix(data = x, missing = NA)
  } else {
    xgb.DMatrix(data = x, label = y, missing = NA)
  }
}

parse_softprob <- function(raw_pred, n_rows, label_levels) {
  n_class <- length(label_levels)
  if (is.matrix(raw_pred)) {
    prob <- raw_pred
  } else {
    prob <- matrix(raw_pred, ncol = n_class, byrow = TRUE)
  }
  if (nrow(prob) != n_rows || ncol(prob) != n_class) {
    stop("Unexpected prediction probability shape from xgboost.")
  }
  colnames(prob) <- label_levels
  prob
}

save_optional_plots <- function(metrics_by_fold, classwise_summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  suppressPackageStartupMessages(library(ggplot2))

  mlong <- melt(
    copy(metrics_by_fold),
    id.vars = c("dataset", "feature_representation", "model", "fold"),
    measure.vars = c("Accuracy", "Kappa", "MacroF1"),
    variable.name = "Metric",
    value.name = "Value"
  )
  p1 <- ggplot(mlong, aes(x = Metric, y = Value)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 2) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Score", title = paste0("XGBoost ", FEATURE_REPRESENTATION, " 5-fold CV"))
  ggsave(METRICS_PLOT_FILE, p1, width = 6.5, height = 4.2, dpi = 300)

  p2 <- ggplot(classwise_summary, aes(x = reorder(Class, Mean_Recall), y = Mean_Recall)) +
    geom_col(width = 0.7) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Mean recall", title = paste0("XGBoost ", FEATURE_REPRESENTATION, " recall by class"))
  ggsave(RECALL_PLOT_FILE, p2, width = 6.5, height = 4.8, dpi = 300)

  invisible(TRUE)
}

run_xgb_cv <- function(dt, feature_cols, fold_dt) {
  label_levels <- sort(unique(dt$label))
  n_class <- length(label_levels)
  if (n_class < 2L) {
    stop("At least two classes are required for multiclass XGBoost.")
  }

  params <- XGB_PARAMS
  params$num_class <- n_class

  x_all <- as.matrix(dt[, ..feature_cols])
  storage.mode(x_all) <- "numeric"
  y_all <- match(dt$label, label_levels) - 1L

  oof_list <- list()
  metrics_list <- list()
  classwise_list <- list()

  for (fold_id in sort(unique(fold_dt$fold))) {
    log_msg("Running fold %s.", fold_id)
    test_idx <- which(fold_dt$fold == fold_id)
    train_idx <- setdiff(seq_len(nrow(dt)), test_idx)

    dtrain <- make_dmatrix(x_all[train_idx, , drop = FALSE], y_all[train_idx])
    dtest <- make_dmatrix(x_all[test_idx, , drop = FALSE])

    model <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = NROUNDS,
      verbose = 0
    )

    prob <- parse_softprob(
      raw_pred = predict(model, dtest),
      n_rows = length(test_idx),
      label_levels = label_levels
    )
    pred_idx <- max.col(prob, ties.method = "first")
    pred_labels <- label_levels[pred_idx]
    sorted_prob <- t(apply(prob, 1, sort, decreasing = TRUE))
    max_prob <- sorted_prob[, 1]
    second_prob <- if (ncol(sorted_prob) >= 2L) sorted_prob[, 2] else rep(NA_real_, nrow(sorted_prob))

    oof_dt <- data.table(
      sample = dt$sample[test_idx],
      true_label = dt$label[test_idx],
      predicted_label = pred_labels,
      fold = fold_id,
      max_prob = max_prob,
      second_prob = second_prob,
      margin = max_prob - second_prob
    )
    prob_dt <- as.data.table(prob)
    oof_list[[as.character(fold_id)]] <- cbind(oof_dt, prob_dt)

    met <- calculate_metrics(dt$label[test_idx], pred_labels, label_levels)
    metrics_list[[as.character(fold_id)]] <- cbind(
      data.table(
        dataset = DATASET_NAME,
        feature_representation = FEATURE_REPRESENTATION,
        model = "XGBoost",
        fold = fold_id,
        n_features = length(feature_cols),
        nrounds = NROUNDS
      ),
      met$overall
    )
    classwise_list[[as.character(fold_id)]] <- cbind(
      data.table(
        dataset = DATASET_NAME,
        feature_representation = FEATURE_REPRESENTATION,
        model = "XGBoost",
        fold = fold_id,
        n_features = length(feature_cols)
      ),
      met$classwise
    )
  }

  oof <- rbindlist(oof_list, use.names = TRUE, fill = TRUE)
  metrics_by_fold <- rbindlist(metrics_list, use.names = TRUE, fill = TRUE)
  classwise_by_fold <- rbindlist(classwise_list, use.names = TRUE, fill = TRUE)

  metrics_summary <- metrics_by_fold[, .(
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    SD_Accuracy = sd(Accuracy, na.rm = TRUE),
    Mean_Kappa = mean(Kappa, na.rm = TRUE),
    SD_Kappa = sd(Kappa, na.rm = TRUE),
    Mean_MacroF1 = mean(MacroF1, na.rm = TRUE),
    SD_MacroF1 = sd(MacroF1, na.rm = TRUE),
    n_folds = .N
  ), by = .(dataset, feature_representation, model, n_features, nrounds)]

  classwise_summary <- classwise_by_fold[, .(
    Mean_Precision = mean(Precision, na.rm = TRUE),
    SD_Precision = sd(Precision, na.rm = TRUE),
    Mean_Recall = mean(Recall, na.rm = TRUE),
    SD_Recall = sd(Recall, na.rm = TRUE),
    Mean_F1 = mean(F1, na.rm = TRUE),
    SD_F1 = sd(F1, na.rm = TRUE),
    Mean_Support = mean(Support, na.rm = TRUE),
    n_folds = .N
  ), by = .(dataset, feature_representation, model, n_features, Class)]

  list(
    label_levels = label_levels,
    x_all = x_all,
    y_all = y_all,
    params = params,
    oof = oof,
    metrics_by_fold = metrics_by_fold,
    metrics_summary = metrics_summary,
    classwise_by_fold = classwise_by_fold,
    classwise_summary = classwise_summary
  )
}

train_and_save_final_model <- function(x_all, y_all, label_levels, feature_cols, params, dt) {
  log_msg("Training final model on the full dataset.")
  dtrain <- make_dmatrix(x_all, y_all)
  final_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = NROUNDS,
    verbose = 0
  )

  bundle <- list(
    model = final_model,
    model_type = "XGBoost",
    dataset = DATASET_NAME,
    feature_representation = FEATURE_REPRESENTATION,
    input_file = INPUT_FILE,
    sample_ids = dt$sample,
    y_train = dt$label,
    label_levels = label_levels,
    feature_order = feature_cols,
    missing_code = MISSING_CODE,
    params = params,
    nrounds = NROUNDS,
    created_at = as.character(Sys.time())
  )

  saveRDS(bundle, MODEL_RDS_FILE)
  tryCatch(
    xgb.save(final_model, MODEL_BINARY_FILE),
    error = function(e) warning("Could not save XGBoost binary model: ", conditionMessage(e))
  )
}

write_parameters <- function(params, feature_cols, label_levels) {
  param_dt <- rbindlist(lapply(names(params), function(nm) {
    data.table(parameter = nm, value = as.character(params[[nm]]))
  }))
  extra <- data.table(
    parameter = c("dataset", "feature_representation", "input_file", "n_features", "n_classes", "label_levels", "nrounds", "n_folds", "seed", "missing_code", "shared_fold_file"),
    value = c(DATASET_NAME, FEATURE_REPRESENTATION, INPUT_FILE, length(feature_cols), length(label_levels), paste(label_levels, collapse = ";"), NROUNDS, N_FOLDS, SEED, MISSING_CODE, SHARED_FOLD_FILE)
  )
  fwrite(rbind(param_dt, extra, fill = TRUE), PARAMS_FILE, sep = "\t")
}

write_summary <- function(dt, feature_cols, fold_dt) {
  lines <- c(
    sprintf("XGBoost %s 5-fold CV", FEATURE_REPRESENTATION),
    sprintf("Created at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Input file: %s", INPUT_FILE),
    sprintf("Output directory: %s", OUTPUT_DIR),
    sprintf("Dataset name: %s", DATASET_NAME),
    sprintf("Samples: %d", nrow(dt)),
    sprintf("Features: %d", length(feature_cols)),
    sprintf("Classes: %s", paste(sort(unique(dt$label)), collapse = ", ")),
    sprintf("Shared fold file: %s", SHARED_FOLD_FILE),
    sprintf("Fold assignment file: %s", FOLD_OUT_FILE),
    sprintf("Folds used: %s", paste(sort(unique(fold_dt$fold)), collapse = ", ")),
    sprintf("Model RDS file: %s", MODEL_RDS_FILE),
    sprintf("Model binary file: %s", MODEL_BINARY_FILE),
    sprintf("OOF prediction file: %s", OOF_FILE),
    sprintf("Metrics summary file: %s", METRICS_SUMMARY_FILE),
    "Parameters: max_depth=7, eta=0.08, nrounds=2000, min_child_weight=2, subsample=0.8, colsample_bytree=0.8, nthread=10",
    "Validation: stratified 5-fold cross-validation; no inner validation; no early stopping."
  )
  writeLines(lines, SUMMARY_FILE)
}

main <- function() {
  check_required_packages()
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  log_msg("Starting XGBoost %s analysis.", FEATURE_REPRESENTATION)
  input <- read_input_matrix(INPUT_FILE)
  dt <- input$dt
  feature_cols <- input$feature_cols

  label_counts <- dt[, .N, by = label][order(-N)]
  fwrite(label_counts, LABEL_COUNTS_FILE, sep = "\t")

  fold_dt <- get_or_create_folds(dt$sample, dt$label, SHARED_FOLD_FILE, FOLD_OUT_FILE)
  result <- run_xgb_cv(dt, feature_cols, fold_dt)

  fwrite(result$oof, OOF_FILE, sep = "\t")
  fwrite(result$metrics_by_fold, METRICS_BY_FOLD_FILE, sep = "\t")
  fwrite(result$metrics_summary, METRICS_SUMMARY_FILE, sep = "\t")
  fwrite(result$classwise_by_fold, CLASSWISE_BY_FOLD_FILE, sep = "\t")
  fwrite(result$classwise_summary, CLASSWISE_SUMMARY_FILE, sep = "\t")

  write_parameters(result$params, feature_cols, result$label_levels)
  train_and_save_final_model(result$x_all, result$y_all, result$label_levels, feature_cols, result$params, dt)
  save_optional_plots(result$metrics_by_fold, result$classwise_summary)
  write_summary(dt, feature_cols, fold_dt)

  log_msg("Finished XGBoost %s analysis.", FEATURE_REPRESENTATION)
}

main()
