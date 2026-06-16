#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(30)

base_dir <- "/home/minhao/knn/xinxinxinknn/data/wgmlst/xunlianji"
top_dir <- file.path(base_dir, "top")

infile <- file.path(
  base_dir,
  "wgmlst_9source_filtered_missing99_major995_redundancy99_from_raw.with_label.tsv"
)

ranking_file <- file.path(top_dir, "knn_hamming_feature_importance_R.tsv")
all_stats_file <- file.path(top_dir, "knn_hamming_feature_importance_R.all_loci_stats.tsv")
summary_file <- file.path(top_dir, "knn_hamming_feature_importance_R_summary.txt")

top_loci_dir <- file.path(top_dir, "top_loci")
top_matrix_dir <- file.path(top_dir, "top_matrices")

dir.create(top_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(top_loci_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(top_matrix_dir, recursive = TRUE, showWarnings = FALSE)

unlink(file.path(top_loci_dir, "Top_*_loci.txt"))
unlink(file.path(top_matrix_dir, "wgmlst_top*_knn_hamming_importance_R.with_label.tsv"))

topN_vec <- c(50, 100, 400, 800, 1600, 3586)

logf <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf(...), "\n", sep = "")
  flush.console()
}

entropy_fun <- function(p) {
  p <- p[p > 0 & is.finite(p)]
  if (length(p) == 0) return(0)
  -sum(p * log(p))
}

minmax <- function(x) {
  y <- x
  y[!is.finite(y)] <- NA_real_
  if (all(is.na(y))) return(rep(0, length(x)))
  mn <- min(y, na.rm = TRUE)
  mx <- max(y, na.rm = TRUE)
  if (!is.finite(mn) || !is.finite(mx) || mx == mn) {
    out <- rep(0, length(x))
  } else {
    out <- (y - mn) / (mx - mn)
    out[is.na(out)] <- 0
  }
  out
}

calc_one_locus <- function(v, y, label_levels, locus_name) {
  n_samples <- length(v)
  missing <- is.na(v) | v == -1L
  n_missing <- sum(missing)
  missing_rate <- n_missing / n_samples
  call_rate <- 1 - missing_rate
  valid_idx <- which(!missing)
  n_nonmissing <- length(valid_idx)

  if (n_nonmissing == 0) {
    return(data.table(
      locus = locus_name,
      n_samples = n_samples,
      n_missing = n_missing,
      missing_rate = missing_rate,
      call_rate = call_rate,
      n_nonmissing = 0L,
      n_alleles_nonmissing = 0L,
      major_allele = NA_integer_,
      major_allele_count = 0L,
      major_allele_rate = NA_real_,
      mi = 0,
      nmi = 0,
      within_mismatch_balanced = NA_real_,
      between_mismatch_balanced = NA_real_,
      hamming_separation = 0
    ))
  }

  vv <- v[valid_idx]
  yy <- y[valid_idx]
  allele_tab <- sort(table(vv), decreasing = TRUE)
  n_alleles <- length(allele_tab)
  major_allele <- as.integer(names(allele_tab)[1])
  major_allele_count <- as.integer(allele_tab[1])
  major_allele_rate <- major_allele_count / n_nonmissing

  tab <- table(
    factor(yy, levels = label_levels),
    factor(vv)
  )

  total_valid <- sum(tab)

  if (total_valid == 0) {
    mi <- 0
    nmi <- 0
  } else {
    pxy <- tab / total_valid
    px <- rowSums(pxy)
    py <- colSums(pxy)
    mi <- 0
    nz <- which(pxy > 0, arr.ind = TRUE)

    if (nrow(nz) > 0) {
      for (k in seq_len(nrow(nz))) {
        i <- nz[k, 1]
        j <- nz[k, 2]
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }

    hx <- entropy_fun(px)
    hy <- entropy_fun(py)
    denom <- min(hx, hy)
    nmi <- ifelse(denom > 0, mi / denom, 0)
  }

  label_n <- rowSums(tab)
  within_rates <- c()

  for (g in seq_along(label_levels)) {
    ng <- label_n[g]
    if (ng >= 2) {
      cg <- as.numeric(tab[g, ])
      same_pairs <- sum(cg * (cg - 1)) / 2
      total_pairs <- ng * (ng - 1) / 2
      mismatch_rate <- 1 - same_pairs / total_pairs
      within_rates <- c(within_rates, mismatch_rate)
    }
  }

  within_balanced <- ifelse(length(within_rates) > 0, mean(within_rates), NA_real_)
  between_rates <- c()

  if (length(label_levels) >= 2) {
    for (g1 in 1:(length(label_levels) - 1)) {
      for (g2 in (g1 + 1):length(label_levels)) {
        n1 <- label_n[g1]
        n2 <- label_n[g2]
        if (n1 > 0 && n2 > 0) {
          c1 <- as.numeric(tab[g1, ])
          c2 <- as.numeric(tab[g2, ])
          same_cross <- sum(c1 * c2)
          total_cross <- n1 * n2
          mismatch_rate <- 1 - same_cross / total_cross
          between_rates <- c(between_rates, mismatch_rate)
        }
      }
    }
  }

  between_balanced <- ifelse(length(between_rates) > 0, mean(between_rates), NA_real_)
  hamming_separation <- between_balanced - within_balanced
  if (!is.finite(hamming_separation)) hamming_separation <- 0

  data.table(
    locus = locus_name,
    n_samples = n_samples,
    n_missing = n_missing,
    missing_rate = missing_rate,
    call_rate = call_rate,
    n_nonmissing = n_nonmissing,
    n_alleles_nonmissing = n_alleles,
    major_allele = major_allele,
    major_allele_count = major_allele_count,
    major_allele_rate = major_allele_rate,
    mi = mi,
    nmi = nmi,
    within_mismatch_balanced = within_balanced,
    between_mismatch_balanced = between_balanced,
    hamming_separation = hamming_separation
  )
}

logf("Input file: %s", infile)

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

logf("Reading training wgMLST matrix")
dt <- fread(infile, sep = "\t", data.table = TRUE, showProgress = TRUE)

if (ncol(dt) < 3) {
  stop("Input file must contain accession, label, and at least one locus column.")
}

if (!identical(names(dt)[1:2], c("accession", "label"))) {
  stop("Input file must have accession and label as the first two columns.")
}

labels <- dt$label
label_levels <- sort(unique(labels))
loci <- names(dt)[-(1:2)]

logf("Samples: %d", nrow(dt))
logf("Loci: %d", length(loci))
logf("Labels: %s", paste(label_levels, collapse = ", "))

label_counts <- dt[, .N, by = label][order(-N)]
print(label_counts)

res_list <- vector("list", length(loci))
logf("Calculating locus-level statistics")

for (i in seq_along(loci)) {
  locus <- loci[i]
  v <- as.integer(dt[[locus]])
  res_list[[i]] <- calc_one_locus(
    v = v,
    y = labels,
    label_levels = label_levels,
    locus_name = locus
  )

  if (i %% 100 == 0 || i == length(loci)) {
    logf("Processed loci: %d/%d", i, length(loci))
  }
}

ranking <- rbindlist(res_list)
ranking[, eligible_for_ranking := TRUE]
ranking[, nmi_component := minmax(nmi)]
ranking[, hamming_separation_positive := pmax(hamming_separation, 0)]
ranking[, hamming_separation_component := minmax(hamming_separation_positive)]
ranking[, call_rate_component := call_rate]
ranking[, core_score :=
          0.45 * nmi_component +
          0.45 * hamming_separation_component +
          0.10 * call_rate_component]
ranking[, final_score := core_score * sqrt(call_rate_component)]
ranking[!is.finite(final_score), final_score := 0]
ranking[!is.finite(core_score), core_score := 0]

setorder(
  ranking,
  -final_score,
  -core_score,
  -hamming_separation_component,
  -nmi_component,
  -call_rate_component,
  locus
)

ranking[, rank := .I]

setcolorder(ranking, c(
  "rank",
  "locus",
  "final_score",
  "core_score",
  "nmi_component",
  "hamming_separation_component",
  "call_rate_component",
  "mi",
  "nmi",
  "within_mismatch_balanced",
  "between_mismatch_balanced",
  "hamming_separation",
  "hamming_separation_positive",
  "missing_rate",
  "call_rate",
  "n_samples",
  "n_missing",
  "n_nonmissing",
  "n_alleles_nonmissing",
  "major_allele",
  "major_allele_count",
  "major_allele_rate",
  "eligible_for_ranking"
))

all_stats <- copy(ranking)
all_stats[, all_loci_order := rank]

logf("Writing ranking file: %s", ranking_file)
fwrite(ranking, ranking_file, sep = "\t", quote = FALSE)

logf("Writing all-loci stats file: %s", all_stats_file)
fwrite(all_stats, all_stats_file, sep = "\t", quote = FALSE)

topN_vec <- topN_vec[topN_vec <= nrow(ranking)]
if (!nrow(ranking) %in% topN_vec) {
  topN_vec <- c(topN_vec, nrow(ranking))
}
topN_vec <- sort(unique(topN_vec))

for (N in topN_vec) {
  top_loci <- ranking[1:N, locus]

  loci_file <- file.path(top_loci_dir, sprintf("Top_%d_loci.txt", N))
  fwrite(data.table(locus = top_loci), loci_file, sep = "\t", col.names = FALSE, quote = FALSE)

  out_mat <- file.path(
    top_matrix_dir,
    sprintf("wgmlst_top%d_knn_hamming_importance_R.with_label.tsv", N)
  )

  sub_dt <- dt[, c("accession", "label", top_loci), with = FALSE]
  fwrite(sub_dt, out_mat, sep = "\t", quote = FALSE)

  logf("Top-%d loci and matrix written", N)
}

summary_lines <- c(
  "===== KNN-Hamming feature importance ranking summary =====",
  sprintf("Input file: %s", infile),
  sprintf("Output directory: %s", top_dir),
  sprintf("Samples: %d", nrow(dt)),
  sprintf("Loci: %d", length(loci)),
  sprintf("Ranked loci: %d", nrow(ranking)),
  "Minimum call-rate hard filter: not applied",
  sprintf("Labels: %s", paste(label_levels, collapse = ", ")),
  "",
  "Feature importance rule:",
  "core_score = 0.45*NMI_component + 0.45*Hamming_separation_component + 0.10*call_rate",
  "final_score = core_score * sqrt(call_rate)",
  "",
  "Definitions:",
  "NMI_component: min-max normalized mutual information between allele state and source label.",
  "Hamming_separation_component: min-max normalized positive value of balanced between-source mismatch rate minus balanced within-source mismatch rate.",
  "call_rate: 1 - missing_rate. Missing calls are encoded as -1.",
  "All loci are retained for ranking; call_rate is used only as a score component and missingness penalty.",
  "",
  sprintf("Ranking file: %s", ranking_file),
  sprintf("All-loci stats file: %s", all_stats_file),
  sprintf("Top loci directory: %s", top_loci_dir),
  sprintf("Top matrices directory: %s", top_matrix_dir),
  "",
  "Generated Top-N sets:",
  paste(topN_vec, collapse = ", "),
  "",
  "Label counts:"
)

label_lines <- paste(label_counts$label, label_counts$N, sep = "\t")
writeLines(c(summary_lines, label_lines), summary_file)

cat("\n===== Top 30 loci preview =====\n")
print(ranking[1:min(30, .N), .(
  rank,
  locus,
  final_score,
  core_score,
  nmi,
  hamming_separation,
  missing_rate,
  call_rate,
  n_alleles_nonmissing,
  major_allele_rate
)])

cat("\n===== DONE =====\n")
cat("Ranking file: ", ranking_file, "\n", sep = "")
cat("All-loci stats file: ", all_stats_file, "\n", sep = "")
cat("Summary file: ", summary_file, "\n", sep = "")
cat("Top loci dir: ", top_loci_dir, "\n", sep = "")
cat("Top matrices dir: ", top_matrix_dir, "\n", sep = "")
