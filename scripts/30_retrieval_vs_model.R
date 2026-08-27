# separates the two changes that shipped together in the revised pipeline
#
# the revised pipeline widened the candidate step and added a trained scoring
# model at the same time, so its gain over the original rules cannot be
# attributed to either one. this holds the scoring model fixed and varies only
# how the candidates were selected.
#
# narrow_* reproduces the original steps 1 to 3 on the feature frame:
#   simrel is each pair's similarity divided by the best similarity for that
#   ICD-9-CM code, which is exactly what the original similarity cutoff
#   thresholds on. cooc_rank <= N is step 2. chapter_ok is step 3.
# wide is what the revised pipeline actually uses.
#
# usage: Rscript 30_retrieval_vs_model.R [track]

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "results"
N_OUTER <- 5
N_INNER <- 3

WIDE <- list(`10_9` = list(K = 10, N = 50), `8_9` = list(K = 10, N = 10))

SETTINGS <- list(
  list(name = "narrow_thr0.995_N30", kind = "narrow", thr = 0.995, N = 30, chapter = TRUE),
  list(name = "narrow_thr0.99_N10",  kind = "narrow", thr = 0.990, N = 10, chapter = TRUE),
  list(name = "narrow_thr0.95_N30",  kind = "narrow", thr = 0.950, N = 30, chapter = TRUE),
  list(name = "wide",                kind = "wide")
)

prf <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp/(tp+fp) else 0
  r <- if (tp + fn > 0) tp/(tp+fn) else 0
  c(precision = p, recall = r,
    f1 = if (p + r > 0) 2*p*r/(p+r) else 0,
    accuracy = if (tp+fp+fn > 0) tp/(tp+fp+fn) else 0)
}
score_pairs <- function(emitted, truth, codes) {
  e <- emitted %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  t <- truth   %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  tp <- nrow(dplyr::semi_join(e, t, by = c("ICD_9_CM","target")))
  prf(tp, nrow(e) - tp, nrow(t) - tp)
}
make_folds <- function(codes, k) {
  codes <- sample(codes); split(codes, rep(seq_len(k), length.out = length(codes)))
}
emit_from_scores <- function(df, tau, rho) {
  df %>% group_by(ICD_9_CM) %>% mutate(.mx = max(.p_score)) %>% ungroup() %>%
    filter(.p_score >= tau | .p_score >= rho * .mx | .p_score == .mx) %>%
    select(ICD_9_CM, target)
}
train_ranker <- function(train_df, fcols) {
  d <- xgb.DMatrix(as.matrix(train_df[, fcols]), label = train_df$y)
  xgb.train(params = list(objective = "binary:logistic",
                          max_depth = 6, eta = 0.1, subsample = 0.8,
                          colsample_bytree = 0.8, min_child_weight = 5,
                          nthread = parallel::detectCores()),
            data = d, nrounds = 200, verbose = 0)
}
fit_chapter_encoding <- function(train_df) {
  prior <- mean(train_df$y)
  enc <- train_df %>% group_by(chapter_pair) %>%
    summarise(n = n(), pos = sum(y), .groups = "drop") %>%
    mutate(chapter_rate = (pos + 20 * prior) / (n + 20)) %>%
    select(chapter_pair, chapter_rate)
  list(enc = enc, prior = prior)
}
apply_chapter_encoding <- function(df, ce) {
  df %>% left_join(ce$enc, by = "chapter_pair") %>%
    mutate(chapter_rate = ifelse(is.na(chapter_rate), ce$prior, chapter_rate))
}

select_candidates <- function(df, st, wide) {
  if (st$kind == "wide") {
    rank_cols <- grep("^simrank_", names(df), value = TRUE)
    keep <- Reduce(`|`, lapply(rank_cols, function(c) df[[c]] <= wide$K)) |
            (df$cooc_rank <= wide$N)
  } else {
    rel_cols <- grep("^simrel_", names(df), value = TRUE)
    keep <- Reduce(`|`, lapply(rel_cols, function(c) df[[c]] >= st$thr)) |
            (df$cooc_rank <= st$N)
    if (isTRUE(st$chapter)) keep <- keep & (df$chapter_ok == 1)
  }
  df[keep, , drop = FALSE]
}

args <- commandArgs(trailingOnly = TRUE)
tracks <- if (length(args)) args else c("10_9", "8_9")
rows <- list()

for (tr in tracks) {
  cat(sprintf("\n################ TRACK %s ################\n", tr))
  base <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr))) %>%
    filter(has_truth == 1)
  # denominator is the full pool, so pairs that candidate selection never
  # retrieved still count as misses. without this a narrower setting would
  # look better simply by being scored on fewer answers
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  codes <- sort(unique(base$ICD_9_CM))
  n_true <- nrow(truth)

  for (st in SETTINGS) {
    feat <- select_candidates(base, st, WIDE[[tr]])
    reachable <- nrow(dplyr::semi_join(truth, feat, by = c("ICD_9_CM","target")))
    fcols <- setdiff(names(feat),
      c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair","CCS_ID",
        "chapter_icd9","chapter_target"))
    fcols <- c(fcols, "chapter_rate")

    set.seed(20260819)
    outer_folds <- make_folds(codes, N_OUTER)
    emitted_all <- list()
    for (fi in seq_along(outer_folds)) {
      test_codes  <- outer_folds[[fi]]
      train_codes <- setdiff(codes, test_codes)
      train_full  <- feat %>% filter(ICD_9_CM %in% train_codes)
      if (!nrow(train_full) || sum(train_full$y) == 0) next

      inner <- make_folds(train_codes, N_INNER)
      isc <- list()
      for (ii in seq_along(inner)) {
        va <- inner[[ii]]; trn <- setdiff(train_codes, va)
        a <- train_full %>% filter(ICD_9_CM %in% trn)
        b <- train_full %>% filter(ICD_9_CM %in% va)
        if (!nrow(a) || !nrow(b) || sum(a$y) == 0) next
        ce <- fit_chapter_encoding(a)
        a <- apply_chapter_encoding(a, ce); b <- apply_chapter_encoding(b, ce)
        m <- train_ranker(a, fcols)
        b$.p_score <- predict(m, as.matrix(b[, fcols]))
        isc[[length(isc)+1]] <- b
      }
      isc <- bind_rows(isc)
      if (!nrow(isc)) next
      best <- list(f1 = -1, tau = 0.5, rho = 0.7)
      for (tau in seq(0.05, 0.95, by = 0.05)) for (rho in c(0.3, 0.5, 0.7, 0.85, 1.0)) {
        m <- score_pairs(emit_from_scores(isc, tau, rho), truth, train_codes)
        if (m["f1"] > best$f1) best <- list(f1 = m["f1"], tau = tau, rho = rho)
      }

      ce <- fit_chapter_encoding(train_full)
      trn <- apply_chapter_encoding(train_full, ce)
      te  <- apply_chapter_encoding(feat %>% filter(ICD_9_CM %in% test_codes), ce)
      if (!nrow(te)) next
      m <- train_ranker(trn, fcols)
      te$.p_score <- predict(m, as.matrix(te[, fcols]))
      emitted_all[[fi]] <- emit_from_scores(te, best$tau, best$rho)
    }

    mm <- score_pairs(bind_rows(emitted_all), truth, codes)
    cat(sprintf("  %-22s candidates %6d | correct pairs kept %5.1f%% | F1 %.3f\n",
                st$name, nrow(feat), 100 * reachable / n_true, mm["f1"]))
    rows[[length(rows)+1]] <- tibble(track = tr, setting = st$name,
      n_candidates = nrow(feat), pct_true_reachable = round(100 * reachable / n_true, 1),
      precision = mm["precision"], recall = mm["recall"],
      f1 = mm["f1"], accuracy = mm["accuracy"])
  }
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, "retrieval_vs_model.csv"), row.names = FALSE)

cat("\n\n===== scoring model held fixed, candidate selection varied =====\n\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr)
  w <- d$f1[d$setting == "wide"]
  cat(sprintf("-- %s (wide = %.3f) --\n", tr, w))
  print(as.data.frame(d %>% mutate(delta_vs_wide = round(f1 - w, 4)) %>%
    select(setting, n_candidates, pct_true_reachable, precision, recall, f1, delta_vs_wide) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))))
  cat("\n")
}
