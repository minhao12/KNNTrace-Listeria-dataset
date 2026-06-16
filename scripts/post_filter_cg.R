#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(30)

# ============================================================
# Config
# ============================================================
tag <- "cg"
dataset_label <- "cgMLST-1748"

out_dir <- "/home/minhao/knn/xinxinxinknn/result/renyuancg/normalizedshaixuanhou"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

human_distance_file <- "/home/minhao/knn/xinxinxinknn/result/juli/cgnormalized/cg_human_to_predicted_source_d1_meanTop11_normalized_all.tsv"
threshold_file <- "/home/minhao/knn/xinxinxinknn/data/yuzhi/50cgnormalized/cg_sameSource_CCge50_meanTop11_normalized_thresholds_by_source.tsv"

threshold_cols <- c("q1_value", "median_value", "mean_value", "q90_value", "q95_value", "q975_value")

logf <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf(...))
  cat("\n")
  flush.console()
}

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  as.character(x) %in% c("TRUE", "True", "true", "T", "1", "YES", "Yes", "yes")
}

safe_rate <- function(a, b) {
  ifelse(b > 0, a / b, NA_real_)
}

# ============================================================
# Read files
# ============================================================
for (ff in c(human_distance_file, threshold_file)) {
  if (!file.exists(ff)) stop("[ERROR] Missing file: ", ff)
}

logf("Dataset: %s", dataset_label)
logf("Human normalized distance file: %s", human_distance_file)
logf("Normalized threshold file: %s", threshold_file)
logf("Output dir: %s", out_dir)

dt <- fread(human_distance_file)
thr <- fread(threshold_file)

required_dt_cols <- c(
  "sample",
  "predicted_source",
  "confidence_pass_0.5",
  "max_prob",
  "d1_normalized",
  "mean_top11_normalized"
)

miss_dt <- setdiff(required_dt_cols, names(dt))
if (length(miss_dt) > 0) {
  stop("[ERROR] Human distance file missing columns: ", paste(miss_dt, collapse = ", "))
}

miss_thr <- setdiff(c("source", threshold_cols), names(thr))
if (length(miss_thr) > 0) {
  stop("[ERROR] Threshold file missing columns: ", paste(miss_thr, collapse = ", "))
}

dt[, confidence_pass_0.5 := parse_bool(`confidence_pass_0.5`)]
dt[, predicted_source := as.character(predicted_source)]

thr_keep <- thr[, c("source", threshold_cols), with = FALSE]
setnames(thr_keep, threshold_cols, paste0("threshold_", threshold_cols))

m <- merge(
  dt,
  thr_keep,
  by.x = "predicted_source",
  by.y = "source",
  all.x = TRUE
)

if (any(is.na(m$threshold_q95_value))) {
  missing_sources <- unique(m[is.na(threshold_q95_value), predicted_source])
  writeLines(missing_sources, file.path(out_dir, paste0(tag, "_missing_threshold_sources.log")))
  stop("[ERROR] Some predicted sources do not have thresholds.")
}

distance_metrics <- c(
  d1 = "d1_normalized",
  meanTop11 = "mean_top11_normalized"
)

threshold_names <- sub("_value$", "", threshold_cols)

# ============================================================
# Long comparison table
# ============================================================
long_list <- list()
idx <- 1L

for (metric_name in names(distance_metrics)) {
  dcol <- distance_metrics[[metric_name]]

  for (i in seq_along(threshold_cols)) {
    th_col <- paste0("threshold_", threshold_cols[i])
    th_name <- threshold_names[i]

    tmp <- m[, .(
      dataset = dataset_label,
      sample = sample,
      input_label = if ("input_label" %in% names(m)) input_label else NA_character_,
      predicted_source = predicted_source,
      predicted_source_from_model = if ("predicted_source_from_model" %in% names(m)) predicted_source_from_model else predicted_source,
      predicted_label_conf05 = if ("predicted_label_conf05" %in% names(m)) predicted_label_conf05 else NA_character_,
      confidence_pass_0.5 = `confidence_pass_0.5`,
      max_prob = max_prob,
      second_prob = if ("second_prob" %in% names(m)) second_prob else NA_real_,
      margin = if ("margin" %in% names(m)) margin else NA_real_,
      distance_metric = metric_name,
      threshold_name = th_name,
      distance_value = get(dcol),
      threshold_value = get(th_col),
      d1_normalized = d1_normalized,
      mean_top11_normalized = mean_top11_normalized,
      d1_dismatch_count = if ("d1_dismatch_count" %in% names(m)) d1_dismatch_count else NA_real_,
      mean_top11_dismatch_count = if ("mean_top11_dismatch_count" %in% names(m)) mean_top11_dismatch_count else NA_real_,
      mean_top11_valid_loci = if ("mean_top11_valid_loci" %in% names(m)) mean_top11_valid_loci else NA_real_
    )]

    tmp[, distance_pass := !is.na(distance_value) & !is.na(threshold_value) & distance_value <= threshold_value]
    tmp[, pass_both_conf_and_distance := confidence_pass_0.5 == TRUE & distance_pass == TRUE]
    tmp[, final_status := fifelse(pass_both_conf_and_distance, "Predictable", "Unpredictable")]

    long_list[[idx]] <- tmp
    idx <- idx + 1L
  }
}

long_dt <- rbindlist(long_list, use.names = TRUE, fill = TRUE)

# ============================================================
# Summary tables
# ============================================================
pass_summary <- long_dt[
  ,
  .(
    n_total = .N,
    n_confident = sum(confidence_pass_0.5 == TRUE, na.rm = TRUE),
    n_unconfident = sum(confidence_pass_0.5 != TRUE | is.na(confidence_pass_0.5), na.rm = TRUE),
    n_distance_evaluable = sum(!is.na(distance_value) & !is.na(threshold_value)),
    n_distance_pass = sum(distance_pass == TRUE, na.rm = TRUE),
    n_pass_both = sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE),
    pass_rate_total = safe_rate(sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE), .N),
    pass_rate_among_confident = safe_rate(
      sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE),
      sum(confidence_pass_0.5 == TRUE, na.rm = TRUE)
    ),
    distance_pass_rate_total = safe_rate(sum(distance_pass == TRUE, na.rm = TRUE), .N)
  ),
  by = .(dataset, distance_metric, threshold_name)
][order(distance_metric, threshold_name)]

pass_by_source <- long_dt[
  ,
  .(
    n_total = .N,
    n_confident = sum(confidence_pass_0.5 == TRUE, na.rm = TRUE),
    n_distance_pass = sum(distance_pass == TRUE, na.rm = TRUE),
    n_pass_both = sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE),
    pass_rate_total = safe_rate(sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE), .N),
    pass_rate_among_confident = safe_rate(
      sum(pass_both_conf_and_distance == TRUE, na.rm = TRUE),
      sum(confidence_pass_0.5 == TRUE, na.rm = TRUE)
    )
  ),
  by = .(dataset, distance_metric, threshold_name, predicted_source)
][order(distance_metric, threshold_name, predicted_source)]

# ============================================================
# Wide flags
# ============================================================
wide <- copy(m)

for (metric_name in names(distance_metrics)) {
  dcol <- distance_metrics[[metric_name]]

  for (i in seq_along(threshold_cols)) {
    th_col <- paste0("threshold_", threshold_cols[i])
    th_name <- threshold_names[i]

    flag_col <- paste0("pass_", metric_name, "_normalized_le_", th_name)
    status_col <- paste0("final_", metric_name, "_normalized_", th_name)

    wide[, (flag_col) := (
      confidence_pass_0.5 == TRUE &
        !is.na(get(dcol)) &
        !is.na(get(th_col)) &
        get(dcol) <= get(th_col)
    )]

    wide[, (status_col) := fifelse(get(flag_col) == TRUE, "Predictable", "Unpredictable")]
  }
}

# ============================================================
# Save outputs
# ============================================================
prefix <- paste0(tag, "_human_conf05_normalized_threshold_filter")

fwrite(wide, file.path(out_dir, paste0(prefix, "_wide_flags.tsv")), sep = "\t")
fwrite(long_dt, file.path(out_dir, paste0(prefix, "_long.tsv")), sep = "\t")
fwrite(pass_summary, file.path(out_dir, paste0(prefix, "_pass_summary.tsv")), sep = "\t")
fwrite(pass_by_source, file.path(out_dir, paste0(prefix, "_pass_by_predicted_source.tsv")), sep = "\t")

# Save retained samples for each threshold and distance metric
for (metric_name in names(distance_metrics)) {
  for (th_name in threshold_names) {
    retained <- long_dt[
      distance_metric == metric_name &
        threshold_name == th_name &
        pass_both_conf_and_distance == TRUE
    ]

    fwrite(
      retained,
      file.path(out_dir, paste0(prefix, "_retained_", metric_name, "_", th_name, ".tsv")),
      sep = "\t"
    )
  }
}

# Main q95 analysis: save both D1 and meanTop11 results
q95_main <- long_dt[
  threshold_name == "q95" &
    distance_metric %in% c("d1", "meanTop11")
]

fwrite(
  q95_main,
  file.path(out_dir, paste0(prefix, "_q95_main_D1_and_meanTop11_long.tsv")),
  sep = "\t"
)

q95_d1 <- q95_main[distance_metric == "d1"]
q95_mean <- q95_main[distance_metric == "meanTop11"]

q95_d1_keep <- q95_d1[, .(
  sample,
  d1_normalized_value = distance_value,
  d1_normalized_threshold = threshold_value,
  d1_normalized_pass_distance = distance_pass,
  d1_normalized_pass_both = pass_both_conf_and_distance,
  d1_normalized_final_status = final_status
)]

q95_mean_keep <- q95_mean[, .(
  sample,
  meanTop11_normalized_value = distance_value,
  meanTop11_normalized_threshold = threshold_value,
  meanTop11_normalized_pass_distance = distance_pass,
  meanTop11_normalized_pass_both = pass_both_conf_and_distance,
  meanTop11_normalized_final_status = final_status
)]

q95_wide <- merge(q95_d1_keep, q95_mean_keep, by = "sample", all = TRUE)

meta_keep <- unique(q95_main[, .(
  sample,
  input_label,
  predicted_source,
  predicted_source_from_model,
  predicted_label_conf05,
  confidence_pass_0.5,
  max_prob,
  second_prob,
  margin
)])

q95_wide <- merge(meta_keep, q95_wide, by = "sample", all.x = TRUE)

fwrite(
  q95_wide,
  file.path(out_dir, paste0(prefix, "_q95_main_D1_and_meanTop11_wide.tsv")),
  sep = "\t"
)

# Main meanTop11 q95 retained samples: save separately
main_meanTop11_q95 <- long_dt[
  distance_metric == "meanTop11" &
    threshold_name == "q95" &
    pass_both_conf_and_distance == TRUE
]

fwrite(
  main_meanTop11_q95,
  file.path(out_dir, paste0(prefix, "_MAIN_retained_meanTop11_q95.tsv")),
  sep = "\t"
)

run_summary <- data.table(
  item = c(
    "dataset",
    "human_distance_file",
    "threshold_file",
    "n_human_samples",
    "n_threshold_sources",
    "distance_metrics",
    "thresholds_tested",
    "main_rule",
    "out_dir"
  ),
  value = c(
    dataset_label,
    human_distance_file,
    threshold_file,
    as.character(nrow(dt)),
    as.character(nrow(thr)),
    paste(names(distance_metrics), collapse = ","),
    paste(threshold_names, collapse = ","),
    "confidence_pass_0.5 == TRUE and mean_top11_normalized <= predicted-source-specific q95 threshold",
    out_dir
  )
)

fwrite(run_summary, file.path(out_dir, paste0(prefix, "_run_summary.tsv")), sep = "\t")

cat("\n===== Pass summary =====\n")
print(pass_summary)

cat("\n===== q95 pass by source =====\n")
print(pass_by_source[threshold_name == "q95"][order(distance_metric, predicted_source)])

cat("\n[DONE] Outputs saved to:\n")
cat(out_dir, "\n")
