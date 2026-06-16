#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c(
    "data.table", "caret", "glmnet", "Matrix",
    "doParallel", "foreach", "parallel", "ggplot2"
  )
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Please install them before running this script."
    )
  }
  library(data.table)
  library(caret)
  library(glmnet)
  library(Matrix)
  library(doParallel)
  library(foreach)
  library(parallel)
  library(ggplot2)
})

options(glmnet.parallel = FALSE)

config <- list(
  dataset_name = "wgmlst3586",
  input_file = "/home/minhao/knn/xinxinxinknn/data/wgmlst/xunlianji/top/top_matrices/wgmlst_top3586_knn_hamming_importance_R.with_label.tsv",
  output_dir = "/home/minhao/knn/xinxinxinknn/moxing/glm/wgmlst3586",
  shared_fold_file = "/home/minhao/knn/xinxinxinknn/daima/glmnet_shared_5fold_assignment.tsv",
  alpha = 1,
  top_k_alleles = 10L,
  lambda_rule = "min",
  nlambda = 60L,
  outer_folds = 5L,
  inner_folds = 5L,
  n_cores = 30L,
  random_seed = 42L,
  missing_codes = c("", "-", "NA", "N/A", "NaN", "NULL", "-1")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) config$input_file <- args[1]
if (length(args) >= 2 && nzchar(args[2])) config$output_dir <- args[2]
if (length(args) >= 3 && nzchar(args[3])) config$shared_fold_file <- args[3]

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

set.seed(config$random_seed)

log_message <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...), "\n", sep = ""
  )
  flush.console()
}

normalize_allele <- function(x, missing_codes) {
  y <- trimws(as.character(x))
  y[is.na(y)] <- ""
  y[toupper(y) %in% toupper(missing_codes)] <- "MISSING"
  y[grepl("^[A-Za-z]", y) & y != "MISSING"] <- "MISSING"
  y[y == ""] <- "MISSING"
  y
}

make_topk_mapper <- function(x_train, k, missing_codes) {
  x <- normalize_allele(x_train, missing_codes)
  tab <- sort(table(x[x != "MISSING"]), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(k, length(tab)))]
  levels_out <- c(keep, "OTHER", "MISSING")

  function(x_new) {
    y <- normalize_allele(x_new, missing_codes)
    y <- ifelse(y %in% keep, y, ifelse(y == "MISSING", "MISSING", "OTHER"))
    factor(y, levels = levels_out)
  }
}

class_weights <- function(y_factor) {
  tab <- table(y_factor)
  inv <- 1 / as.numeric(tab)
  names(inv) <- names(tab)
  w <- inv[as.character(y_factor)]
  as.numeric(w / mean(w))
}

align_sparse_columns <- function(mat, target_cols) {
  missing_cols <- setdiff(target_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    zero_block <- Matrix::Matrix(0, nrow(mat), length(missing_cols), sparse = TRUE)
    colnames(zero_block) <- missing_cols
    mat <- cbind(mat, zero_block)
  }
  mat[, target_cols, drop = FALSE]
}

calculate_metrics <- function(predicted, observed) {
  predicted <- factor(predicted, levels = levels(observed))
  observed <- factor(observed, levels = levels(observed))
  cm <- caret::confusionMatrix(predicted, observed)

  by_class <- as.data.table(cm$byClass, keep.rownames = TRUE)
  setnames(by_class, "rn", "Class")
  by_class[, Class := sub("^Class: ", "", Class)]

  if (!"F1" %in% names(by_class) &&
      all(c("Sensitivity", "Pos Pred Value") %in% names(by_class))) {
    by_class[, F1 := 2 * (`Pos Pred Value` * Sensitivity) /
               pmax(`Pos Pred Value` + Sensitivity, 1e-9)]
  }

  out_class <- by_class[, .(
    Class,
    Precision = fifelse(is.na(`Pos Pred Value`), 0, `Pos Pred Value`),
    Recall = fifelse(is.na(Sensitivity), 0, Sensitivity),
    F1 = fifelse(is.na(F1), 0, F1)
  )]

  list(
    accuracy = as.numeric(cm$overall["Accuracy"]),
    kappa = as.numeric(cm$overall["Kappa"]),
    macro_f1 = mean(out_class$F1),
    by_class = out_class
  )
}

load_or_create_folds <- function(dt, config) {
  if (file.exists(config$shared_fold_file)) {
    fold_dt <- fread(config$shared_fold_file)
    required_cols <- c("sample", "label", "fold")
    if (!all(required_cols %in% names(fold_dt))) {
      stop("Invalid fold file. Required columns: sample, label, fold")
    }
    missing_samples <- setdiff(dt$sample, fold_dt$sample)
    if (length(missing_samples) > 0) {
      stop("Fold file does not contain all input samples. Missing examples: ",
           paste(head(missing_samples, 10), collapse = ", "))
    }
    fold_dt <- fold_dt[match(dt$sample, fold_dt$sample)]
    folds <- split(seq_len(nrow(dt)), fold_dt$fold)
    folds <- folds[order(as.integer(names(folds)))]
    return(list(folds = folds, fold_assignment = fold_dt))
  }

  folds <- caret::createFolds(
    dt$label,
    k = config$outer_folds,
    list = TRUE,
    returnTrain = FALSE
  )
  fold_assignment <- rbindlist(lapply(seq_along(folds), function(i) {
    data.table(
      sample = dt$sample[folds[[i]]],
      label = as.character(dt$label[folds[[i]]]),
      fold = i
    )
  }))
  dir.create(dirname(config$shared_fold_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(fold_assignment, config$shared_fold_file, sep = "\t")
  list(folds = folds, fold_assignment = fold_assignment)
}

run_glmnet_cv <- function(dt, feature_cols, config) {
  label_levels <- levels(dt$label)
  fold_info <- load_or_create_folds(dt, config)
  folds <- fold_info$folds

  cl <- parallel::makeCluster(min(config$n_cores, length(folds)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  parallel::clusterExport(
    cl,
    c(
      "dt", "feature_cols", "label_levels", "folds", "config",
      "make_topk_mapper", "class_weights", "align_sparse_columns",
      "calculate_metrics", "normalize_allele"
    ),
    envir = environment()
  )

  foreach(
    fold_id = seq_along(folds),
    .packages = c("data.table", "caret", "glmnet", "Matrix")
  ) %dopar% {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(dt)), test_idx)

    mappers <- lapply(feature_cols, function(col) {
      make_topk_mapper(
        dt[[col]][train_idx],
        k = config$top_k_alleles,
        missing_codes = config$missing_codes
      )
    })
    names(mappers) <- feature_cols

    encode_matrix <- function(idx) {
      encoded_df <- setNames(
        lapply(feature_cols, function(col) mappers[[col]](dt[[col]][idx])),
        feature_cols
      )
      encoded_df <- as.data.frame(encoded_df, stringsAsFactors = TRUE)
      Matrix::sparse.model.matrix(~ . - 1, data = encoded_df)
    }

    x_train <- encode_matrix(train_idx)
    x_test <- encode_matrix(test_idx)
    all_columns <- union(colnames(x_train), colnames(x_test))
    x_train <- align_sparse_columns(x_train, all_columns)
    x_test <- align_sparse_columns(x_test, all_columns)

    y_train <- factor(dt$label[train_idx], levels = label_levels)
    y_test <- factor(dt$label[test_idx], levels = label_levels)
    weights_train <- class_weights(y_train)

    fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "multinomial",
      alpha = config$alpha,
      nlambda = config$nlambda,
      type.multinomial = "ungrouped",
      standardize = FALSE,
      weights = weights_train,
      parallel = FALSE,
      nfolds = config$inner_folds
    )

    lambda_value <- if (config$lambda_rule == "1se") fit$lambda.1se else fit$lambda.min

    pred_class <- factor(
      as.vector(predict(fit, newx = x_test, s = lambda_value, type = "class")),
      levels = label_levels
    )

    prob_arr <- predict(fit, newx = x_test, s = lambda_value, type = "response")
    prob_mat <- as.matrix(prob_arr[, , 1])
    colnames(prob_mat) <- label_levels

    metrics <- calculate_metrics(pred_class, y_test)

    oof <- data.table(
      sample = dt$sample[test_idx],
      true_label = as.character(y_test),
      predicted_label = as.character(pred_class),
      fold = fold_id,
      lambda_used = as.numeric(lambda_value)
    )
    oof <- cbind(oof, as.data.table(prob_mat))

    list(
      metrics = data.table(
        dataset = config$dataset_name,
        n_features = length(feature_cols),
        model = "GLMNET",
        fold = fold_id,
        alpha = config$alpha,
        top_k_alleles = config$top_k_alleles,
        lambda_rule = config$lambda_rule,
        lambda_used = as.numeric(lambda_value),
        Accuracy = metrics$accuracy,
        Kappa = metrics$kappa,
        MacroF1 = metrics$macro_f1
      ),
      classwise = metrics$by_class[, `:=`(
        dataset = config$dataset_name,
        n_features = length(feature_cols),
        model = "GLMNET",
        fold = fold_id,
        alpha = config$alpha,
        top_k_alleles = config$top_k_alleles,
        lambda_rule = config$lambda_rule
      )],
      oof = oof
    )
  }
}

write_outputs <- function(results, dt, feature_cols, config) {
  prefix <- paste0("glmnet_", config$dataset_name)

  metrics_by_fold <- rbindlist(lapply(results, `[[`, "metrics"))
  classwise_by_fold <- rbindlist(lapply(results, `[[`, "classwise"))
  oof_predictions <- rbindlist(lapply(results, `[[`, "oof"))

  metrics_summary <- metrics_by_fold[, .(
    Mean_Accuracy = mean(Accuracy),
    SD_Accuracy = sd(Accuracy),
    Mean_Kappa = mean(Kappa),
    SD_Kappa = sd(Kappa),
    Mean_MacroF1 = mean(MacroF1),
    SD_MacroF1 = sd(MacroF1),
    n_folds = .N
  ), by = .(dataset, n_features, model, alpha, top_k_alleles, lambda_rule)]

  classwise_summary <- classwise_by_fold[, .(
    Mean_Precision = mean(Precision, na.rm = TRUE),
    SD_Precision = sd(Precision, na.rm = TRUE),
    Mean_Recall = mean(Recall, na.rm = TRUE),
    SD_Recall = sd(Recall, na.rm = TRUE),
    Mean_F1 = mean(F1, na.rm = TRUE),
    SD_F1 = sd(F1, na.rm = TRUE),
    n_folds = .N
  ), by = .(dataset, n_features, model, alpha, top_k_alleles, lambda_rule, Class)][order(Class)]

  label_counts <- dt[, .N, by = label][order(-N)]
  parameter_table <- data.table(
    dataset = config$dataset_name,
    model = "GLMNET",
    input_file = config$input_file,
    output_dir = config$output_dir,
    shared_fold_file = config$shared_fold_file,
    samples = nrow(dt),
    features = length(feature_cols),
    classes = length(levels(dt$label)),
    alpha = config$alpha,
    top_k_alleles = config$top_k_alleles,
    lambda_rule = config$lambda_rule,
    nlambda = config$nlambda,
    outer_folds = config$outer_folds,
    inner_folds = config$inner_folds,
    n_cores = config$n_cores,
    random_seed = config$random_seed
  )

  fwrite(label_counts, file.path(config$output_dir, paste0(prefix, "_label_counts.tsv")), sep = "\t")
  fwrite(metrics_by_fold, file.path(config$output_dir, paste0(prefix, "_metrics_by_fold.tsv")), sep = "\t")
  fwrite(metrics_summary, file.path(config$output_dir, paste0(prefix, "_metrics_summary.tsv")), sep = "\t")
  fwrite(classwise_by_fold, file.path(config$output_dir, paste0(prefix, "_classwise_by_fold.tsv")), sep = "\t")
  fwrite(classwise_summary, file.path(config$output_dir, paste0(prefix, "_classwise_summary.tsv")), sep = "\t")
  fwrite(oof_predictions, file.path(config$output_dir, paste0(prefix, "_oof_predictions.tsv")), sep = "\t")
  fwrite(parameter_table, file.path(config$output_dir, paste0(prefix, "_parameters.tsv")), sep = "\t")

  metric_long <- melt(
    metrics_by_fold,
    id.vars = c("dataset", "n_features", "model", "fold"),
    measure.vars = c("Accuracy", "Kappa", "MacroF1"),
    variable.name = "Metric",
    value.name = "Value"
  )

  p_metrics <- ggplot(metric_long, aes(x = Metric, y = Value)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 2) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste0("GLMNET ", config$dataset_name, " 5-fold CV"), x = NULL, y = "Score") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(config$output_dir, paste0(prefix, "_metrics_5fold.png")),
    p_metrics,
    width = 8,
    height = 5,
    dpi = 300
  )

  p_recall <- ggplot(classwise_summary, aes(x = Class, y = Mean_Recall)) +
    geom_col(width = 0.7) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste0("GLMNET ", config$dataset_name, " class-wise recall"), x = "Class", y = "Mean recall") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(config$output_dir, paste0(prefix, "_classwise_recall.png")),
    p_recall,
    width = 10,
    height = 6,
    dpi = 300
  )

  writeLines(c(
    paste0("GLMNET ", config$dataset_name, " 5-fold cross-validation summary"),
    sprintf("Created at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Input file: %s", config$input_file),
    sprintf("Output directory: %s", config$output_dir),
    sprintf("Shared fold file: %s", config$shared_fold_file),
    sprintf("Samples: %d", nrow(dt)),
    sprintf("Features: %d", length(feature_cols)),
    sprintf("Classes: %s", paste(levels(dt$label), collapse = ", ")),
    sprintf("alpha: %s", config$alpha),
    sprintf("top_k_alleles: %d", config$top_k_alleles),
    sprintf("lambda_rule: %s", config$lambda_rule),
    sprintf("nlambda: %d", config$nlambda),
    sprintf("outer_folds: %d", config$outer_folds),
    sprintf("inner_folds: %d", config$inner_folds),
    sprintf("n_cores: %d", config$n_cores),
    "Output files:",
    paste0(prefix, "_label_counts.tsv"),
    paste0(prefix, "_metrics_by_fold.tsv"),
    paste0(prefix, "_metrics_summary.tsv"),
    paste0(prefix, "_classwise_by_fold.tsv"),
    paste0(prefix, "_classwise_summary.tsv"),
    paste0(prefix, "_oof_predictions.tsv"),
    paste0(prefix, "_parameters.tsv"),
    paste0(prefix, "_metrics_5fold.png"),
    paste0(prefix, "_classwise_recall.png")
  ), file.path(config$output_dir, paste0(prefix, "_run_summary.txt")))

  list(metrics_summary = metrics_summary, classwise_summary = classwise_summary)
}

log_message("Starting GLMNET %s 5-fold CV", config$dataset_name)
log_message("Input file: %s", config$input_file)
log_message("Output directory: %s", config$output_dir)

if (!file.exists(config$input_file)) {
  stop("Input file not found: ", config$input_file)
}

dt <- fread(config$input_file)
if (ncol(dt) < 3) stop("Input file must contain at least three columns.")
setnames(dt, names(dt)[1:2], c("sample", "label"))
dt[, sample := as.character(sample)]
dt[, label := factor(label)]
feature_cols <- setdiff(names(dt), c("sample", "label"))

log_message("Samples: %d", nrow(dt))
log_message("Features: %d", length(feature_cols))
log_message("Classes: %s", paste(levels(dt$label), collapse = ", "))

if (length(feature_cols) != 3586L) {
  warning(sprintf("Expected 3586 wgMLST features, but found %d features.", length(feature_cols)))
}

results <- run_glmnet_cv(dt, feature_cols, config)
summary_outputs <- write_outputs(results, dt, feature_cols, config)

log_message("Completed GLMNET %s 5-fold CV", config$dataset_name)
log_message("Metrics summary: %s", file.path(config$output_dir, paste0("glmnet_", config$dataset_name, "_metrics_summary.tsv")))
print(summary_outputs$metrics_summary)
