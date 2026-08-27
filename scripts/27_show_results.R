source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

suppressMessages({library(dplyr)})

R <- "results"
get <- function(f) {
  p <- file.path(R, f)
  if (!file.exists(p)) { cat(sprintf("  [missing %s, rerun the script that makes it]\n", f)); return(NULL) }
  read.csv(p, stringsAsFactors = FALSE)
}
hdr <- function(n, s) cat(sprintf("\n\n=========== %s. %s ===========\n", n, s))
tname <- function(x) ifelse(x == "10_9", "ICD-9 to ICD-10-CA", "ICD-9 to ICDA-8")

# ---------------------------------------------------------------- 1
hdr(1, "where the project started and where it is now")
d <- get("cv_rerank_results.csv")
if (!is.null(d)) {
  # the f1 column is the 3rd of each 4-column metric block
  f1 <- function(row, pre) {
    cols <- grep(pre, names(d), fixed = TRUE, value = TRUE)
    as.numeric(row[cols[3]])
  }
  for (tr in c("10_9", "8_9")) {
    x <- d[d$track == tr, ]
    cat(sprintf("\n%s\n", tname(tr)))
    cat(sprintf("  %-42s %.3f\n", "original rules, scored on its own data",
                f1(x[x$method == "baseline_in_sample", ], "m_bl_insample")))
    cat(sprintf("  %-42s %.3f\n", "original rules, scored on unseen codes",
                f1(x[x$method == "baseline_heldout_cv", ], "m_bl_all")))
    cat(sprintf("  %-42s %.3f\n", "two stage system, scored on unseen codes",
                f1(x[x$method == "rerank_heldout_cv", ], "m_rr_all")))
  }
  cat("\n  f1 is a score out of 1. higher is better.\n")
}

# ---------------------------------------------------------------- 2
hdr(2, "the main finding, where correct answers were being lost")
d <- get("error_analysis_summary.csv")
if (!is.null(d)) {
  x <- d[d$model == "SapBERT", ]
  cat("\n  of every correct answer in the manual crosswalk, the share that\n")
  cat("  survives each stage of the original pipeline:\n\n")
  cat(sprintf("  %-36s %10s %10s\n", "", "ICD-10-CA", "ICDA-8"))
  # pct_in_universe is 100 by construction, the matrix scores every code
  # against every target. it is a check that no validation pair names a code
  # missing from the label files, not a stage anything can be lost at
  rows <- list(
    c("survives the similarity cutoff",  "pct_pass_threshold"),
    c("in the co-occurrence list",       "pct_in_cooc"),
    c("kept as a candidate",             "pct_in_pool"),
    c("actually output",                 "pct_emitted"))
  for (r in rows)
    cat(sprintf("  %-36s %9.1f%% %9.1f%%\n", r[1],
                x[x$track == "10_9", r[2]], x[x$track == "8_9", r[2]]))
  cat(sprintf("\n  %-36s %10.3f %10.3f\n", "best f1 possible from what survives",
              x[x$track == "10_9", "oracle_f1"], x[x$track == "8_9", "oracle_f1"]))
  cat(sprintf("\n  the cutoff drops %.0f%% of correct ICD-10-CA pairs before anything is\n",
              100 - x[x$track == "10_9", "pct_pass_threshold"]))
  cat("  ranked, and a dropped pair cannot come back. that capped the pipeline.\n")
  cat("  widening the candidate step recovers them, which is what says the\n")
  cat("  retrieval step was the problem and not the model.\n")
}

# ---------------------------------------------------------------- 3
hdr(3, "embedding model comparison")
d <- get("full_grid_best.csv")
if (!is.null(d)) {
  for (tr in c("10_9", "8_9")) {
    cat(sprintf("\n%s\n", tname(tr)))
    x <- d[d$track == tr, ] %>% arrange(desc(f1))
    for (i in seq_len(nrow(x)))
      cat(sprintf("  %-38s f1 %.3f\n", x$model[i], x$f1[i]))
  }
  cat("\n  mpnet has no medical training and still matches sapbert on ICD-10-CA.\n")
}

# ---------------------------------------------------------------- 4
hdr(4, "is the right answer reachable at all")
d <- get("top1_accuracy.csv")
if (!is.null(d)) {
  cat(sprintf("\n  %-22s %8s %8s %8s\n", "", "top-1", "top-3", "top-5"))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-22s %7.1f%% %7.1f%% %7.1f%%\n", tname(d$track[i]),
                100*d$top1_accuracy[i], 100*d$top3_accuracy[i], 100*d$top5_accuracy[i]))
  cat("\n  this is whether the correct target is ranked first, separate from\n")
  cat("  whether the output rule decides to keep it.\n")
}

# ---------------------------------------------------------------- 5
hdr(5, "what happens at a 95% precision target")
d <- get("precision_coverage_triage.csv")
if (!is.null(d)) {
  for (tr in c("10_9", "8_9")) {
    x <- d[d$track == tr, ]
    cat(sprintf("\n%s\n", tname(tr)))
    for (i in seq_len(nrow(x)))
      cat(sprintf("  %-42s %4d codes  %5.1f%%\n", x$bucket[i], x$n[i], x$pct[i]))
  }
}

# ---------------------------------------------------------------- 6
hdr(6, "which parts of the system earn their place")
for (tr in c("10_9", "8_9")) {
  d <- get(sprintf("ablation_%s.csv", tr))
  if (is.null(d)) next
  full <- d$f1[d$variant == "full"]
  cat(sprintf("\n%s   (full system = %.3f)\n", tname(tr), full))
  # the rank_* rows swap the scoring method, they do not remove a feature
  x <- d[grepl("^no_", d$variant), ] %>% mutate(delta = f1 - full) %>% arrange(delta)
  for (i in seq_len(nrow(x)))
    cat(sprintf("  remove %-20s f1 %.3f   %+.3f\n",
                sub("^no_", "", x$variant[i]), x$f1[i], x$delta[i]))
  y <- d[grepl("^rank_", d$variant), ] %>% mutate(delta = f1 - full)
  for (i in seq_len(nrow(y)))
    cat(sprintf("  (scoring method %-11s f1 %.3f   %+.3f)\n",
                sub("^rank_", "", y$variant[i]), y$f1[i], y$delta[i]))
}
cat("\n  most negative = most important. the reverse direction features matter\n")
cat("  most on ICD-10-CA. clinicalbert earns nothing once the others are in.\n")

# ---------------------------------------------------------------- 7
hdr(7, "would more training data help")
for (tr in c("10_9", "8_9")) {
  d <- get(sprintf("learning_curve_%s.csv", tr))
  if (is.null(d)) next
  cat(sprintf("\n%s\n", tname(tr)))
  # the largest size can be empty when the track has fewer codes than that
  d <- d[!is.na(d$f1_mean), ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %4d training codes   f1 %.3f  (sd %.3f)\n",
                d$n_train[i], d$f1_mean[i], d$f1_sd[i]))
}
cat("\n  it flattens out. more manually mapped codes of the same kind would\n")
cat("  not move the result.\n")

# ---------------------------------------------------------------- 8
hdr(8, "does it work on clinical areas it has never seen")
for (tr in c("10_9", "8_9")) {
  d <- get(sprintf("category_holdout_%s.csv", tr))
  if (is.null(d)) next
  cat(sprintf("\n%s\n", tname(tr)))
  for (i in seq_len(nrow(d)))
    cat(sprintf("  split by %-10s f1 %.3f  (sd %.3f)\n",
                d$split[i], d$f1_mean[i], d$f1_sd[i]))
}
cat("\n  holding out whole CCS categories costs almost nothing, so it is not\n")
cat("  memorising categories.\n")

# ---------------------------------------------------------------- 9
hdr(9, "how much would these numbers move on a different set of codes")
d <- get("bootstrap_variability.csv")
if (!is.null(d)) {
  d <- d[d$metric == "f1", ]
  for (tr in c("10_9", "8_9")) {
    x <- d[d$track == tr, ]
    cat(sprintf("\n%s\n", tname(tr)))
    for (i in seq_len(nrow(x)))
      cat(sprintf("  %-22s f1 %.3f  (sd %.3f, 95%% %.3f to %.3f)\n",
                  x$model[i], x$estimate[i], x$sd[i], x$ci_low[i], x$ci_high[i]))
  }
  cat("\n  2000 bootstrap draws, resampling whole icd-9 codes. the sd is about\n")
  cat("  0.015 and 0.022, so the ~0.10 gain from sapbert is far bigger than noise.\n")
}

# ---------------------------------------------------------------- 10
hdr(10, "how evenly it works across the 130 ccs categories")
d <- get("ccs_category_summary.csv")
if (!is.null(d)) {
  for (tr in c("10_9", "8_9")) {
    x <- d[d$track == tr, ]
    cat(sprintf("\n%s\n", tname(tr)))
    for (i in seq_len(nrow(x)))
      cat(sprintf("  %-22s mean %.3f  sd %.3f  median %.3f  perfect %d  zero %d\n",
                  x$model[i], x$mean_f1[i], x$sd_f1[i], x$median_f1[i],
                  x$n_perfect[i], x$n_zero[i]))
  }
  cat("\n  spread between categories (sd ~0.25) dwarfs the noise on the overall\n")
  cat("  number. 59 of 130 categories hold a single code, so treat those gently.\n")
}

cat("\n\n=========== related scripts ===========\n")
cat("  24_show_similarity_matrix.R   what a similarity table looks like\n")
cat("  25_frequency_distributions.R  max similarity and co-occurrence spread\n")
cat("  26_stopword_choice.R          which stopword dictionary to use\n")
cat("  28_stopwords_and_codes.R      stop words in/out x code number in/out\n")
cat("  21_error_by_code_type.R       single vs multi target breakdown\n")
