# produces the actual crosswalk
#
# 12_ measures how good the system is, this is the system. trains on every code
# with a known answer then maps the ones without, which is what youd run on a
# real legacy code set
#
# each mapping gets a confidence and a triage label. the threshold comes from
# 13_ which derives it from out of fold predictions, so it isnt tuned on
# anything this script sees

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "results"
TRACKS <- c("10_9", "8_9")

curve <- read.csv(file.path(OUT_DIR, "precision_coverage_curve.csv"))
folds <- read.csv(file.path(OUT_DIR, "cv_rerank_folds.csv"))

out <- list(); summ <- list()

for (tr in TRACKS) {
  feat_all <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  train <- feat_all %>% filter(has_truth == 1)
  topred <- feat_all %>% filter(has_truth == 0)

  cat(sprintf("\n=== %s: training on %d codes, predicting %d codes ===\n",
              tr, n_distinct(train$ICD_9_CM), n_distinct(topred$ICD_9_CM)))
  if (nrow(topred) == 0) { cat("  nothing to predict\n"); next }

  fc <- folds %>% filter(track == tr)
  cfg <- list(K = round(median(fc$K)), N = round(median(fc$N)),
              chapter = mean(fc$chapter) > 0.5)
  rho <- median(fc$rho)
  cat(sprintf("  retrieval K=%d N=%d chapter=%s (median across CV folds)\n",
              cfg$K, cfg$N, cfg$chapter))

  ok <- curve %>% filter(track == tr, !is.na(precision), precision >= 0.95)
  tau <- if (nrow(ok)) (ok %>% slice_max(n_emitted, n = 1, with_ties = FALSE))$tau else 0.9
  cat(sprintf("  auto-accept threshold %.2f (95%% precision point from CV)\n", tau))

  FEATURE_COLS <- setdiff(names(train),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair","CCS_ID",
      "chapter_icd9","chapter_target"))

  prior <- mean(train$y)
  enc <- train %>% group_by(chapter_pair) %>%
    summarise(n = n(), pos = sum(y), .groups = "drop") %>%
    mutate(chapter_rate = (pos + 20 * prior) / (n + 20)) %>%
    select(chapter_pair, chapter_rate)
  add_enc <- function(d) d %>% left_join(enc, by = "chapter_pair") %>%
    mutate(chapter_rate = ifelse(is.na(chapter_rate), prior, chapter_rate))
  train <- add_enc(train); topred <- add_enc(topred)

  fcols <- c(FEATURE_COLS, "chapter_rate")
  bst <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "logloss",
                  max_depth = 6, eta = 0.1, subsample = 0.8,
                  colsample_bytree = 0.8, min_child_weight = 5,
                  nthread = parallel::detectCores()),
    data = xgb.DMatrix(as.matrix(train[, fcols]), label = train$y),
    nrounds = 200, verbose = 0)

  topred$confidence <- predict(bst, as.matrix(topred[, fcols]))

  rank_cols <- grep("^simrank_", names(topred), value = TRUE)
  keep <- Reduce(`|`, lapply(rank_cols, function(c) topred[[c]] <= cfg$K)) |
          (topred$cooc_rank <= cfg$N)
  if (cfg$chapter) keep <- keep & (topred$chapter_ok == 1)
  elig <- topred[keep, , drop = FALSE]

  pred <- elig %>% group_by(ICD_9_CM) %>% mutate(.mx = max(confidence)) %>% ungroup() %>%
    filter(confidence >= tau | confidence >= rho * .mx | confidence == .mx) %>%
    mutate(triage = ifelse(confidence >= tau, "auto-accept", "review")) %>%
    arrange(ICD_9_CM, desc(confidence)) %>%
    transmute(track = tr, ICD_9_CM, icd9_label, target, target_label,
              confidence = round(confidence, 4), triage)

  missing <- setdiff(unique(topred$ICD_9_CM), unique(pred$ICD_9_CM))
  if (length(missing)) {
    pred <- bind_rows(pred, tibble(
      track = tr, ICD_9_CM = missing, icd9_label = NA_character_,
      target = NA_character_, target_label = NA_character_,
      confidence = NA_real_, triage = "no-candidate"))
  }

  bycode <- pred %>% group_by(ICD_9_CM) %>%
    summarise(best = ifelse(all(is.na(confidence)), "no-candidate",
                            ifelse(any(triage == "auto-accept"), "auto-accept", "review")),
              .groups = "drop")
  tb <- bycode %>% count(best)
  cat(sprintf("  %d mappings for %d codes\n", sum(!is.na(pred$confidence)), n_distinct(pred$ICD_9_CM)))
  for (i in seq_len(nrow(tb))) cat(sprintf("    %-14s %d codes\n", tb$best[i], tb$n[i]))

  out[[length(out)+1]] <- pred
  summ[[length(summ)+1]] <- tb %>% mutate(track = tr)
}

if (length(out)) {
  res <- bind_rows(out)
  write.csv(res, file.path(OUT_DIR, "predicted_crosswalk.csv"), row.names = FALSE)
  write.csv(bind_rows(summ), file.path(OUT_DIR, "predicted_crosswalk_summary.csv"), row.names = FALSE)
  cat(sprintf("\nWrote results/predicted_crosswalk.csv (%d rows)\n", nrow(res)))
} else {
  cat("\nNo codes without a known answer -- nothing written.\n")
}
