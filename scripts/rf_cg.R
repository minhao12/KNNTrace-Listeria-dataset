#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_pkgs <- c("data.table", "caret", "ranger", "foreach", "doParallel", "parallel", "ggplot2")
  missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(missing_pkgs) > 0) {
    message("[ERROR] Missing packages: ", paste(missing_pkgs, collapse = ", "))
    quit(status = 1)
  }

  library(data.table)
  library(caret)
  library(ranger)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(ggplot2)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

set.seed(42)

input_file <- "/home/minhao/xiezuo/data/cgmlst/ml_ready_dataset.tsv"
output_dir <- "/home/minhao/knn/xinxinxinknn/moxing/rf/cgmlst1748"
shared_fold_file <- "/home/minhao/knn/xinxinxinknn/daima/rf_shared_5fold_assignment.tsv"
model_file <- "/home/minhao/knn/xinxinxinknn/moxing/rf/cgmlst1748/rf_cgmlst1748_model.rds"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(shared_fold_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_file), recursive = TRUE, showWarnings = FALSE)

rf_num_trees <- 500L
rf_mtry_base <- 582L
rf_min_node_size <- 5L
rf_sample_fraction <- 1
rf_splitrule <- "gini"
n_cores <- 30L

log_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "
", sep = "")
  flush.console()
}

clean_allele_matrix <- function(dt, feature_cols) {
  for (feature in feature_cols) {
    value <- as.character(dt[[feature]])
    value <- trimws(value)
    value <- sub("^INF[-_]", "", value, ignore.case = TRUE)
    invalid <- is.na(value) | value == "" | value == "-" |
      toupper(value) %in% c("NA", "N/A", "NAN", "NULL", "LNF", "NIPH", "NIPHEM", "ALM", "ASM", "PLOT3", "PLOT5", "LOTSC", "PAMA") |
      !grepl("^-?[0-9]+$", value)
    value[invalid] <- NA_character_
    numeric_value <- suppressWarnings(as.integer(value))
    numeric_value[is.na(numeric_value)] <- -1L
    dt[[feature]] <- as.numeric(numeric_value)
  }
  invisible(dt)
}

macro_f1_from_confusion <- function(confusion_obj) {
  by_class <- as.data.table(confusion_obj$byClass, keep.rownames = TRUE)
  setnames(by_class, "rn", "Class")
  if (!"F1" %in% names(by_class) && all(c("Sensitivity", "Pos Pred Value") %in% names(by_class))) {
    by_class[, F1 := 2 * (`Pos Pred Value` * Sensitivity) / pmax(`Pos Pred Value` + Sensitivity, 1e-9)]
  }
  by_class[is.na(F1), F1 := 0]
  mean(by_class$F1)
}

classwise_from_confusion <- function(confusion_obj) {
  by_class <- as.data.table(confusion_obj$byClass, keep.rownames = TRUE)
  setnames(by_class, "rn", "Class")
  if (!"F1" %in% names(by_class) && all(c("Sensitivity", "Pos Pred Value") %in% names(by_class))) {
    by_class[, F1 := 2 * (`Pos Pred Value` * Sensitivity) / pmax(`Pos Pred Value` + Sensitivity, 1e-9)]
  }
  out <- by_class[, .(
    Class,
    Precision = `Pos Pred Value`,
    Recall = Sensitivity,
    F1
  )]
  out[is.na(Precision), Precision := 0]
  out[is.na(Recall), Recall := 0]
  out[is.na(F1), F1 := 0]
  out
}

load_or_create_folds <- function(dt, shared_file, k = 5L) {
  if (file.exists(shared_file)) {
    fold_dt <- fread(shared_file)
    if ("accession" %in% names(fold_dt) && !"sample" %in% names(fold_dt)) {
      setnames(fold_dt, "accession", "sample")
    }
    if (all(c("sample", "fold") %in% names(fold_dt))) {
      matched <- fold_dt[match(dt$sample, fold_dt$sample)]
      if (!any(is.na(matched$fold))) {
        fold_ids <- sort(unique(matched$fold))
        folds <- lapply(fold_ids, function(fold_id) which(matched$fold == fold_id))
        names(folds) <- paste0("Fold", fold_ids)
        assignment <- data.table(
          sample = dt$sample,
          label = as.character(dt$label),
          fold = matched$fold
        )
        return(list(folds = folds, assignment = assignment, source = shared_file))
      }
    }
  }

  folds <- caret::createFolds(dt$label, k = k, list = TRUE, returnTrain = FALSE)
  assignment <- rbindlist(lapply(seq_along(folds), function(i) {
    data.table(
      sample = dt$sample[folds[[i]]],
      label = as.character(dt$label[folds[[i]]]),
      fold = i
    )
  }))
  fwrite(assignment, shared_file, sep = "	")
  list(folds = folds, assignment = assignment, source = "created_by_caret_createFolds")
}

log_message("===== Random Forest cgmlst1748 5-fold CV started =====")
log_message("Input file: %s", input_file)
log_message("Output directory: %s", output_dir)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dt <- fread(input_file)
if (ncol(dt) < 3L) {
  stop("Input file must contain at least three columns: sample, label, and feature columns")
}
setnames(dt, names(dt)[1:2], c("sample", "label"))
dt[, label := factor(label)]

label_levels <- levels(dt$label)
feature_cols <- setdiff(names(dt), c("sample", "label"))
nsamp <- nrow(dt)
nfeat <- length(feature_cols)

if (nfeat != 1748L) {
  stop(sprintf("Expected 1748 features, but found %d", nfeat))
}

rf_mtry <- min(rf_mtry_base, nfeat)

log_message("Samples: %d", nsamp)
log_message("Features: %d", nfeat)
log_message("Classes: %d", length(label_levels))
log_message("Class levels: %s", paste(label_levels, collapse = ", "))
log_message("RF parameters: num.trees=%d, mtry=%d, min.node.size=%d, sample.fraction=%s, splitrule=%s", rf_num_trees, rf_mtry, rf_min_node_size, rf_sample_fraction, rf_splitrule)

label_counts <- dt[, .N, by = label][order(-N)]
fwrite(label_counts, file.path(output_dir, "label_counts_cgmlst1748.tsv"), sep = "	")

dt <- clean_allele_matrix(dt, feature_cols)

fold_obj <- load_or_create_folds(dt, shared_fold_file, k = 5L)
folds <- fold_obj$folds
fold_assignment <- fold_obj$assignment
fold_source <- fold_obj$source

saveRDS(folds, file.path(output_dir, "folds_5cv_cgmlst1748.rds"))
fwrite(fold_assignment, file.path(output_dir, "fold_assignment_cgmlst1748.tsv"), sep = "	")
log_message("Fold source: %s", fold_source)

fold_cores <- min(length(folds), n_cores)
rf_threads <- max(1L, floor(n_cores / fold_cores))
log_message("Parallel setting: fold_cores=%d, ranger_threads_per_fold=%d", fold_cores, rf_threads)

cluster_obj <- makeCluster(fold_cores)
registerDoParallel(cluster_obj)

clusterExport(
  cluster_obj,
  c("dt", "feature_cols", "label_levels", "nsamp", "folds", "rf_num_trees", "rf_mtry", "rf_min_node_size", "rf_sample_fraction", "rf_splitrule", "rf_threads", "macro_f1_from_confusion", "classwise_from_confusion"),
  envir = environment()
)

cv_results <- foreach(
  i = seq_along(folds),
  .packages = c("data.table", "caret", "ranger")
) %dopar% {
  test_index <- folds[[i]]
  train_index <- setdiff(seq_len(nsamp), test_index)

  train_df <- as.data.frame(dt[train_index, c("label", feature_cols), with = FALSE])
  test_x <- as.data.frame(dt[test_index, feature_cols, with = FALSE])
  y_test <- factor(dt$label[test_index], levels = label_levels)

  rf_fit <- ranger(
    dependent.variable.name = "label",
    data = train_df,
    num.trees = rf_num_trees,
    mtry = rf_mtry,
    min.node.size = rf_min_node_size,
    sample.fraction = rf_sample_fraction,
    splitrule = rf_splitrule,
    probability = TRUE,
    classification = TRUE,
    write.forest = TRUE,
    num.threads = rf_threads,
    seed = 1000 + i
  )

  probability_matrix <- predict(rf_fit, test_x)$predictions
  missing_classes <- setdiff(label_levels, colnames(probability_matrix))
  if (length(missing_classes) > 0) {
    add_matrix <- matrix(0, nrow = nrow(probability_matrix), ncol = length(missing_classes))
    colnames(add_matrix) <- missing_classes
    probability_matrix <- cbind(probability_matrix, add_matrix)
  }
  probability_matrix <- probability_matrix[, label_levels, drop = FALSE]

  predicted_label <- factor(colnames(probability_matrix)[max.col(probability_matrix, ties.method = "first")], levels = label_levels)
  confusion_obj <- caret::confusionMatrix(predicted_label, y_test)
  classwise_dt <- classwise_from_confusion(confusion_obj)

  oof_dt <- data.table(
    sample = dt$sample[test_index],
    true_label = as.character(y_test),
    predicted_label = as.character(predicted_label),
    fold = i
  )
  oof_dt <- cbind(oof_dt, as.data.table(probability_matrix))

  list(
    metrics = data.table(
      dataset = "cgmlst1748",
      n_features = length(feature_cols),
      model = "RandomForest",
      fold = i,
      num.trees = rf_num_trees,
      mtry = rf_mtry,
      min.node.size = rf_min_node_size,
      sample.fraction = rf_sample_fraction,
      splitrule = rf_splitrule,
      Accuracy = as.numeric(confusion_obj$overall["Accuracy"]),
      Kappa = as.numeric(confusion_obj$overall["Kappa"]),
      MacroF1 = macro_f1_from_confusion(confusion_obj)
    ),
    classwise = classwise_dt[, `:=`(
      dataset = "cgmlst1748",
      n_features = length(feature_cols),
      model = "RandomForest",
      fold = i,
      num.trees = rf_num_trees,
      mtry = rf_mtry,
      min.node.size = rf_min_node_size,
      sample.fraction = rf_sample_fraction,
      splitrule = rf_splitrule
    )],
    oof = oof_dt
  )
}

stopCluster(cluster_obj)

metrics_by_fold <- rbindlist(lapply(cv_results, `[[`, "metrics"))
classwise_by_fold <- rbindlist(lapply(cv_results, `[[`, "classwise"))
oof_probabilities <- rbindlist(lapply(cv_results, `[[`, "oof"))

fwrite(metrics_by_fold, file.path(output_dir, "metrics_by_fold_cgmlst1748.tsv"), sep = "	")
fwrite(classwise_by_fold, file.path(output_dir, "classwise_by_fold_cgmlst1748.tsv"), sep = "	")
fwrite(oof_probabilities, file.path(output_dir, "oof_probabilities_cgmlst1748.tsv"), sep = "	")

metrics_summary <- metrics_by_fold[, .(
  Mean_Accuracy = mean(Accuracy),
  SD_Accuracy = sd(Accuracy),
  Mean_Kappa = mean(Kappa),
  SD_Kappa = sd(Kappa),
  Mean_MacroF1 = mean(MacroF1),
  SD_MacroF1 = sd(MacroF1),
  n_folds = .N
), by = .(dataset, n_features, model, num.trees, mtry, min.node.size, sample.fraction, splitrule)]

fwrite(metrics_summary, file.path(output_dir, "metrics_summary_cgmlst1748.tsv"), sep = "	")

classwise_summary <- classwise_by_fold[, .(
  Mean_Precision = mean(Precision, na.rm = TRUE),
  SD_Precision = sd(Precision, na.rm = TRUE),
  Mean_Recall = mean(Recall, na.rm = TRUE),
  SD_Recall = sd(Recall, na.rm = TRUE),
  Mean_F1 = mean(F1, na.rm = TRUE),
  SD_F1 = sd(F1, na.rm = TRUE),
  n_folds = .N
), by = .(dataset, n_features, model, num.trees, mtry, min.node.size, sample.fraction, splitrule, Class)][order(Class)]

fwrite(classwise_summary, file.path(output_dir, "classwise_summary_cgmlst1748.tsv"), sep = "	")

parameter_table <- data.table(
  dataset = "cgmlst1748",
  model = "RandomForest",
  num.trees = rf_num_trees,
  mtry = rf_mtry,
  mtry_base = rf_mtry_base,
  min.node.size = rf_min_node_size,
  sample.fraction = rf_sample_fraction,
  splitrule = rf_splitrule,
  outer_cv = "5-fold stratified",
  n_cores = n_cores,
  fold_source = fold_source
)
fwrite(parameter_table, file.path(output_dir, "fixed_parameters_cgmlst1748.tsv"), sep = "	")

metric_long <- melt(
  metrics_by_fold,
  id.vars = c("dataset", "n_features", "model", "fold"),
  measure.vars = c("Accuracy", "Kappa", "MacroF1"),
  variable.name = "Metric",
  value.name = "Value"
)

metric_plot <- ggplot(metric_long, aes(x = Metric, y = Value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Random Forest cgmlst1748 5-fold CV", x = NULL, y = "Score") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "rf_cgmlst1748_5fold_metrics.png"), metric_plot, width = 8, height = 5, dpi = 300)

recall_plot <- ggplot(classwise_summary, aes(x = Class, y = Mean_Recall)) +
  geom_col(width = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Random Forest cgmlst1748 class-wise recall", x = "Class", y = "Mean Recall") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "rf_cgmlst1748_recall_by_class.png"), recall_plot, width = 10, height = 6, dpi = 300)

log_message("Training final deployment model on the full dataset")
full_train_df <- as.data.frame(dt[, c("label", feature_cols), with = FALSE])
full_rf_model <- ranger(
  dependent.variable.name = "label",
  data = full_train_df,
  num.trees = rf_num_trees,
  mtry = rf_mtry,
  min.node.size = rf_min_node_size,
  sample.fraction = rf_sample_fraction,
  splitrule = rf_splitrule,
  probability = TRUE,
  classification = TRUE,
  write.forest = TRUE,
  num.threads = n_cores,
  seed = 4242
)

model_bundle <- list(
  model = full_rf_model,
  dataset = "cgmlst1748",
  model_type = "RandomForest_ranger_probability",
  feature_order = feature_cols,
  label_levels = label_levels,
  missing_code = -1L,
  sample_ids = dt$sample,
  rf_parameters = list(
    num.trees = rf_num_trees,
    mtry = rf_mtry,
    mtry_base = rf_mtry_base,
    min.node.size = rf_min_node_size,
    sample.fraction = rf_sample_fraction,
    splitrule = rf_splitrule
  ),
  training_input_file = input_file,
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
saveRDS(model_bundle, model_file)

summary_file <- file.path(output_dir, "run_summary_cgmlst1748.txt")
writeLines(c(
  "===== Random Forest 5-fold CV summary =====",
  sprintf("Created at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Dataset: %s", "cgmlst1748"),
  sprintf("Input file: %s", input_file),
  sprintf("Output directory: %s", output_dir),
  sprintf("Shared fold file: %s", shared_fold_file),
  sprintf("Fold source: %s", fold_source),
  sprintf("Samples: %d", nsamp),
  sprintf("Features: %d", nfeat),
  sprintf("Classes: %d", length(label_levels)),
  sprintf("Class levels: %s", paste(label_levels, collapse = ", ")),
  "",
  "Fixed RF parameters:",
  sprintf("num.trees = %d", rf_num_trees),
  sprintf("mtry = %d", rf_mtry),
  sprintf("mtry_base = %d", rf_mtry_base),
  sprintf("min.node.size = %d", rf_min_node_size),
  sprintf("sample.fraction = %s", rf_sample_fraction),
  sprintf("splitrule = %s", rf_splitrule),
  "",
  "Output files:",
  sprintf("metrics_by_fold_%s.tsv", "cgmlst1748"),
  sprintf("metrics_summary_%s.tsv", "cgmlst1748"),
  sprintf("classwise_by_fold_%s.tsv", "cgmlst1748"),
  sprintf("classwise_summary_%s.tsv", "cgmlst1748"),
  sprintf("oof_probabilities_%s.tsv", "cgmlst1748"),
  sprintf("fold_assignment_%s.tsv", "cgmlst1748"),
  sprintf("fixed_parameters_%s.tsv", "cgmlst1748"),
  sprintf("rf_%s_5fold_metrics.png", "cgmlst1748"),
  sprintf("rf_%s_recall_by_class.png", "cgmlst1748"),
  sprintf("Model file: %s", model_file)
), summary_file)

log_message("===== Random Forest cgmlst1748 DONE =====")
log_message("Metrics summary: %s", file.path(output_dir, "metrics_summary_cgmlst1748.tsv"))
log_message("Model file: %s", model_file)
print(metrics_summary)
