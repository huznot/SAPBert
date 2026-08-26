# how much do the reported numbers move if the code set had been slightly
# different. the headline f1 is one number off a few hundred icd-9 codes with no
# indication of how stable it is, so this puts a standard deviation on it and on
# the gap between the two stage system and the original rules
#
# method: cluster bootstrap. resample the icd-9 codes with replacement, keeping
# all of a code's pairs together, recompute the metric, repeat B times. the sd
# of those B values is the standard error of the reported number. codes are the
# sampling unit because pairs from the same code are not independent
#
# everything is rebuilt from the out of fold predictions saved by 12_, so the
# numbers here are the same held out numbers, not a refit
#
# usage: Rscript 31_variability.R
source("pipeline_lib.R")
set.seed(42)

OUT_DIR <- "../results"
B <- 2000
TRACKS <- c("10_9", "8_9")

preds <- readRDS(file.path(OUT_DIR, "cv_rerank_predictions.rds"))
folds <- read.csv(file.path(OUT_DIR, "cv_rerank_folds.csv"), stringsAsFactors = FALSE)
reported <- read.csv(file.path(OUT_DIR, "cv_rerank_results.csv"), stringsAsFactors = FALSE)

# same rule as 12_: emit a candidate if it is confident on its own, or close
# enough to the best candidate for that code, and always emit the best one
emit_from_scores <- function(df, tau, rho) {
  df %>% group_by(ICD_9_CM) %>%
    mutate(.mx = max(.p_score)) %>% ungroup() %>%
    filter(.p_score >= tau | .p_score >= rho * .mx | .p_score == .mx) %>%
    distinct(ICD_9_CM, target)
}

counts_by_code <- function(emitted, truth, codes) {
  e <- emitted %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  t <- truth   %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  tp <- dplyr::semi_join(e, t, by = c("ICD_9_CM", "target"))
  bind_rows(
    tp %>% transmute(code = ICD_9_CM, TP = 1L, FP = 0L, FN = 0L),
    dplyr::anti_join(e, t, by = c("ICD_9_CM", "target")) %>%
      transmute(code = ICD_9_CM, TP = 0L, FP = 1L, FN = 0L),
    dplyr::anti_join(t, e, by = c("ICD_9_CM", "target")) %>%
      transmute(code = ICD_9_CM, TP = 0L, FP = 0L, FN = 1L)
  ) %>%
    group_by(code) %>%
    summarise(across(c(TP, FP, FN), sum), .groups = "drop")
}

metrics_from_totals <- function(tp, fp, fn) {
  p <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  r <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  list(precision = p, recall = r,
       f1 = ifelse(p + r > 0, 2 * p * r / (p + r), 0),
       accuracy = ifelse(tp + fp + fn > 0, tp / (tp + fp + fn), 0))
}

boot_metrics <- function(counts, idx) {
  tp <- colSums(matrix(counts$TP[idx], nrow = nrow(idx)))
  fp <- colSums(matrix(counts$FP[idx], nrow = nrow(idx)))
  fn <- colSums(matrix(counts$FN[idx], nrow = nrow(idx)))
  as.data.frame(metrics_from_totals(tp, fp, fn))
}

boot_summary <- function(point, boot) {
  tibble(estimate = round(point, 4), sd = round(sd(boot), 4),
         ci_low = round(quantile(boot, 0.025, names = FALSE), 4),
         ci_high = round(quantile(boot, 0.975, names = FALSE), 4))
}

# rebuilds the held out output of both systems, one row per icd-9 code
held_out_counts <- function(tr) {
  feat  <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr))) %>%
    filter(has_truth == 1)
  truth <- feat %>% filter(y == 1) %>% distinct(ICD_9_CM, target)
  codes <- sort(unique(feat$ICD_9_CM))
  base_emit <- readRDS(file.path(OUT_DIR, sprintf("base_emit_%s.rds", tr)))

  p  <- preds %>% filter(track == tr)
  fd <- folds %>% filter(track == tr)

  rr <- list(); bl <- list()
  for (fi in sort(unique(p$fold))) {
    row <- fd %>% filter(fold == fi)
    stopifnot(nrow(row) == 1, row$baseline_config %in% names(base_emit))
    d <- p %>% filter(fold == fi)
    test_codes <- unique(d$ICD_9_CM)
    rr[[length(rr) + 1]] <- emit_from_scores(d, row$tau, row$rho)
    bl[[length(bl) + 1]] <- base_emit[[row$baseline_config]] %>%
      filter(ICD_9_CM %in% test_codes) %>% distinct(ICD_9_CM, target)
  }

  list(codes = codes, truth = truth,
       counts = list(`two stage system` = counts_by_code(bind_rows(rr), truth, codes),
                     `original rules`   = counts_by_code(bind_rows(bl), truth, codes)),
       ccs = feat %>% distinct(ICD_9_CM, CCS_ID))
}

var_rows <- list(); delta_rows <- list()

for (tr in TRACKS) {
  cat(sprintf("\n=== track %s ===\n", tr))
  h <- held_out_counts(tr)

  # every code gets a row in both systems so the deltas are paired
  aligned <- lapply(h$counts, function(x)
    tibble(code = h$codes) %>% left_join(x, by = "code") %>%
      mutate(across(c(TP, FP, FN), ~ tidyr::replace_na(.x, 0L))))
  saveRDS(aligned, file.path(OUT_DIR, sprintf("heldout_counts_%s.rds", tr)))

  idx <- matrix(sample.int(length(h$codes), length(h$codes) * B, replace = TRUE),
                nrow = length(h$codes), ncol = B)
  cat(sprintf("  %d codes, %d bootstrap replicates\n", length(h$codes), B))

  point <- list(); boot <- list()
  for (sys in names(aligned)) {
    cnt <- aligned[[sys]]
    point[[sys]] <- metrics_from_totals(sum(cnt$TP), sum(cnt$FP), sum(cnt$FN))
    boot[[sys]]  <- boot_metrics(cnt, idx)
    cat(sprintf("  %-18s f1 %.4f\n", sys, point[[sys]]$f1))
    for (m in names(point[[sys]]))
      var_rows[[length(var_rows) + 1]] <- tibble(track = tr, system = sys, metric = m) %>%
        bind_cols(boot_summary(point[[sys]][[m]], boot[[sys]][[m]]))
  }

  # the number that matters, is the two stage system really ahead of the rules
  for (m in c("precision", "recall", "f1", "accuracy")) {
    d_boot <- boot[["two stage system"]][[m]] - boot[["original rules"]][[m]]
    delta_rows[[length(delta_rows) + 1]] <- tibble(
      track = tr, system = "two stage system", reference = "original rules", metric = m,
      delta = round(point[["two stage system"]][[m]] - point[["original rules"]][[m]], 4),
      sd = round(sd(d_boot), 4),
      ci_low = round(quantile(d_boot, 0.025, names = FALSE), 4),
      ci_high = round(quantile(d_boot, 0.975, names = FALSE), 4),
      prop_favouring_system = round(mean(d_boot > 0), 4),
      crosses_zero = quantile(d_boot, 0.025, names = FALSE) < 0 &
                     quantile(d_boot, 0.975, names = FALSE) > 0)
  }

  # the rebuild has to land on the numbers 12_ reported, otherwise the emit rule
  # or the fold assignment has drifted
  rep_f1 <- reported %>% filter(track == tr, method == "rerank_heldout_cv")
  rep_f1 <- as.numeric(rep_f1[[grep("m_rr_all", names(reported), fixed = TRUE, value = TRUE)[3]]])
  if (abs(rep_f1 - point[["two stage system"]]$f1) > 5e-4)
    stop(sprintf("rebuilt f1 %.4f does not match reported %.4f on %s",
                 point[["two stage system"]]$f1, rep_f1, tr))
  cat(sprintf("  matches the reported held out f1 of %.3f\n", rep_f1))
}

variability <- bind_rows(var_rows)
deltas      <- bind_rows(delta_rows)
write.csv(variability, file.path(OUT_DIR, "bootstrap_variability.csv"), row.names = FALSE)
write.csv(deltas, file.path(OUT_DIR, "bootstrap_deltas.csv"), row.names = FALSE)

cat("\n=== held out f1 with bootstrap standard deviation ===\n")
print(as.data.frame(variability %>% filter(metric == "f1") %>% select(-metric)))
cat("\n=== two stage system minus the original rules ===\n")
print(as.data.frame(deltas %>% select(-system, -reference)))

# the folds are a second, independent read on the same spread
fold_var <- folds %>%
  select(track, fold, baseline_f1, rerank_f1, baseline_accuracy, rerank_accuracy) %>%
  tidyr::pivot_longer(-c(track, fold), names_to = "metric", values_to = "value") %>%
  group_by(track, metric) %>%
  summarise(n_folds = n(), mean = round(mean(value), 4), sd = round(sd(value), 4),
            min = round(min(value), 4), max = round(max(value), 4), .groups = "drop")
write.csv(fold_var, file.path(OUT_DIR, "cv_fold_variability.csv"), row.names = FALSE)
cat("\n=== across the five cross validation folds ===\n")
print(as.data.frame(fold_var))

cat("\nDone. Wrote bootstrap_variability.csv, bootstrap_deltas.csv,",
    "cv_fold_variability.csv, heldout_counts_<track>.rds\n")
