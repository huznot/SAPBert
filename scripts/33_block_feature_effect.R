# did the block features actually change anything
#
# compares the current run against a saved earlier run, pairing the two inside
# each bootstrap draw. the folds are seeded identically in 12_ so both runs test
# on the same codes, which is what makes the pairing valid
#
# needs results/baseline_prefeat/ holding the earlier cv_rerank_predictions.rds
# and cv_rerank_folds.csv
#
# usage: Rscript 33_block_feature_effect.R
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
set.seed(42)

OUT_DIR <- "results"
OLD_DIR <- file.path(OUT_DIR, "baseline_prefeat")
B <- 2000
TRACKS <- c("10_9", "8_9")

RUNS <- list(
  `with block features` = list(pred = file.path(OUT_DIR, "cv_rerank_predictions.rds"),
                               folds = file.path(OUT_DIR, "cv_rerank_folds.csv")),
  `before block features` = list(pred = file.path(OLD_DIR, "cv_rerank_predictions.rds"),
                                 folds = file.path(OLD_DIR, "cv_rerank_folds.csv"))
)
for (r in RUNS) for (p in r) if (!file.exists(p)) stop("missing ", p)

emit_from_scores <- function(df, tau, rho) {
  df %>% group_by(ICD_9_CM) %>%
    mutate(.mx = max(.p_score)) %>% ungroup() %>%
    filter(.p_score >= tau | .p_score >= rho * .mx | .p_score == .mx) %>%
    distinct(ICD_9_CM, target)
}

counts_by_code <- function(emitted, truth, codes) {
  e <- emitted %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  t <- truth   %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  bind_rows(
    dplyr::semi_join(e, t, by = c("ICD_9_CM","target")) %>%
      transmute(code = ICD_9_CM, TP = 1L, FP = 0L, FN = 0L),
    dplyr::anti_join(e, t, by = c("ICD_9_CM","target")) %>%
      transmute(code = ICD_9_CM, TP = 0L, FP = 1L, FN = 0L),
    dplyr::anti_join(t, e, by = c("ICD_9_CM","target")) %>%
      transmute(code = ICD_9_CM, TP = 0L, FP = 0L, FN = 1L)
  ) %>% group_by(code) %>% summarise(across(c(TP, FP, FN), sum), .groups = "drop")
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

rows <- list(); drows <- list()

for (tr in TRACKS) {
  feat  <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr))) %>%
    filter(has_truth == 1)
  truth <- feat %>% filter(y == 1) %>% distinct(ICD_9_CM, target)
  codes <- sort(unique(feat$ICD_9_CM))

  cnts <- list()
  for (nm in names(RUNS)) {
    p  <- readRDS(RUNS[[nm]]$pred) %>% filter(track == tr)
    fd <- read.csv(RUNS[[nm]]$folds, stringsAsFactors = FALSE) %>% filter(track == tr)
    em <- bind_rows(lapply(sort(unique(p$fold)), function(fi) {
      row <- fd[fd$fold == fi, ]
      emit_from_scores(p %>% filter(fold == fi), row$tau, row$rho)
    }))
    cnts[[nm]] <- tibble(code = codes) %>%
      left_join(counts_by_code(em, truth, codes), by = "code") %>%
      mutate(across(c(TP, FP, FN), ~ tidyr::replace_na(.x, 0L)))
  }

  # the two runs must have tested on the same codes for the pairing to mean
  # anything. 12_ seeds the folds so they do
  same_codes <- identical(cnts[[1]]$code, cnts[[2]]$code)
  if (!same_codes) stop("the two runs do not cover the same codes on ", tr)

  idx <- matrix(sample.int(length(codes), length(codes) * B, replace = TRUE),
                nrow = length(codes), ncol = B)

  pt <- list(); bt <- list()
  for (nm in names(RUNS)) {
    c1 <- cnts[[nm]]
    pt[[nm]] <- metrics_from_totals(sum(c1$TP), sum(c1$FP), sum(c1$FN))
    bt[[nm]] <- boot_metrics(c1, idx)
    rows[[length(rows) + 1]] <- tibble(track = tr, run = nm,
      precision = round(pt[[nm]]$precision, 4), recall = round(pt[[nm]]$recall, 4),
      f1 = round(pt[[nm]]$f1, 4), accuracy = round(pt[[nm]]$accuracy, 4))
  }

  for (m in c("precision", "recall", "f1", "accuracy")) {
    d <- bt[[1]][[m]] - bt[[2]][[m]]
    drows[[length(drows) + 1]] <- tibble(
      track = tr, metric = m,
      delta = round(pt[[1]][[m]] - pt[[2]][[m]], 4),
      sd = round(sd(d), 4),
      ci_low = round(quantile(d, 0.025, names = FALSE), 4),
      ci_high = round(quantile(d, 0.975, names = FALSE), 4),
      prop_favouring_blocks = round(mean(d > 0), 4),
      crosses_zero = quantile(d, 0.025, names = FALSE) < 0 &
                     quantile(d, 0.975, names = FALSE) > 0)
  }
}

both <- bind_rows(rows); deltas <- bind_rows(drows)
write.csv(both, file.path(OUT_DIR, "block_feature_runs.csv"), row.names = FALSE)
write.csv(deltas, file.path(OUT_DIR, "block_feature_effect.csv"), row.names = FALSE)

cat("\n=== both runs, held out ===\n")
print(as.data.frame(both))
cat("\n=== with block features minus without, paired ===\n")
print(as.data.frame(deltas))
cat("\nDone. Wrote block_feature_runs.csv, block_feature_effect.csv\n")
