# is fixing the retrieval config better than tuning it per fold, or did one lucky
# choice make it look that way
#
# 12_ picks K and N inside each fold and gets 0.668 / 0.840. 16_ fixes them and
# gets more, which looks like per fold tuning overfitting on ~276 training codes.
# but the fixed values came from medians of what 12_ picked so theyre mildly
# optimistic
#
# this scores the same model over a spread of fixed configs. if most land in the
# same place the conclusion holds, if only one does it was a lucky pick

source("pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "../results"
N_OUTER <- 5
N_INNER <- 3

CONFIGS <- list()
for (K in c(5, 10, 25, 50)) for (N in c(10, 30, 50)) for (ch in c(FALSE, TRUE))
  CONFIGS[[length(CONFIGS)+1]] <- list(K = K, N = N, chapter = ch)

prf <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp/(tp+fp) else 0
  r <- if (tp + fn > 0) tp/(tp+fn) else 0
  c(precision = p, recall = r, f1 = if (p + r > 0) 2*p*r/(p+r) else 0,
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
apply_retrieval <- function(df, cfg) {
  rc <- grep("^simrank_", names(df), value = TRUE)
  keep <- Reduce(`|`, lapply(rc, function(c) df[[c]] <= cfg$K)) | (df$cooc_rank <= cfg$N)
  if (cfg$chapter) keep <- keep & (df$chapter_ok == 1)
  df[keep, , drop = FALSE]
}
fit_ce <- function(d) {
  pr <- mean(d$y)
  list(enc = d %>% group_by(chapter_pair) %>%
         summarise(n = n(), pos = sum(y), .groups = "drop") %>%
         mutate(chapter_rate = (pos + 20*pr)/(n + 20)) %>% select(chapter_pair, chapter_rate),
       prior = pr)
}
add_ce <- function(d, ce) d %>% left_join(ce$enc, by = "chapter_pair") %>%
  mutate(chapter_rate = ifelse(is.na(chapter_rate), ce$prior, chapter_rate))
train_ranker <- function(d, fc) xgb.train(
  params = list(objective = "binary:logistic", max_depth = 6, eta = 0.1,
                subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5,
                nthread = parallel::detectCores()),
  data = xgb.DMatrix(as.matrix(d[, fc]), label = d$y), nrounds = 200, verbose = 0)

args <- commandArgs(trailingOnly = TRUE)
tracks <- if (length(args)) args else c("10_9", "8_9")
rows <- list()

for (tr in tracks) {
  cat(sprintf("\n######## %s ########\n", tr))
  fa <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  base <- fa %>% filter(has_truth == 1)
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  codes <- sort(unique(base$ICD_9_CM))
  FC <- c(setdiff(names(base),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair",
      "CCS_ID","chapter_icd9","chapter_target")), "chapter_rate")

  set.seed(20260819)
  outer <- make_folds(codes, N_OUTER)

  for (cfg in CONFIGS) {
    feat <- apply_retrieval(base, cfg)
    if (sum(feat$y) < 50) next
    em <- list()
    for (fi in seq_along(outer)) {
      te_codes <- outer[[fi]]; tr_codes <- setdiff(codes, te_codes)
      trf <- feat %>% filter(ICD_9_CM %in% tr_codes)
      inner <- make_folds(tr_codes, N_INNER)
      isc <- list()
      for (ii in seq_along(inner)) {
        va <- inner[[ii]]; tn <- setdiff(tr_codes, va)
        a <- trf %>% filter(ICD_9_CM %in% tn); b <- trf %>% filter(ICD_9_CM %in% va)
        if (!nrow(a) || !nrow(b) || sum(a$y) == 0) next
        ce <- fit_ce(a); a <- add_ce(a, ce); b <- add_ce(b, ce)
        b$.p_score <- predict(train_ranker(a, FC), as.matrix(b[, FC]))
        isc[[length(isc)+1]] <- b
      }
      isc <- bind_rows(isc)
      best <- list(f1 = -1, tau = 0.5, rho = 0.7)
      for (tau in seq(0.05, 0.95, by = 0.05)) for (rho in c(0.3,0.5,0.7,0.85,1.0)) {
        m <- score_pairs(emit_from_scores(isc, tau, rho), truth, tr_codes)
        if (m["f1"] > best$f1) best <- list(f1 = m["f1"], tau = tau, rho = rho)
      }
      ce <- fit_ce(trf)
      a <- add_ce(trf, ce)
      b <- add_ce(feat %>% filter(ICD_9_CM %in% te_codes), ce)
      b$.p_score <- predict(train_ranker(a, FC), as.matrix(b[, FC]))
      em[[fi]] <- emit_from_scores(b, best$tau, best$rho)
    }
    mm <- score_pairs(bind_rows(em), truth, codes)
    cat(sprintf("  K=%-2d N=%-2d chapter=%-5s  P %.3f R %.3f F1 %.3f\n",
                cfg$K, cfg$N, cfg$chapter, mm["precision"], mm["recall"], mm["f1"]))
    rows[[length(rows)+1]] <- tibble(track = tr, K = cfg$K, N = cfg$N,
      chapter = cfg$chapter, precision = mm["precision"], recall = mm["recall"],
      f1 = mm["f1"], accuracy = mm["accuracy"])
  }
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, sprintf("retrieval_sensitivity_%s.csv",
          paste(tracks, collapse = "_"))), row.names = FALSE)

cat("\n===== summary =====\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr)
  cat(sprintf("\n%s: %d configs, F1 median %.3f, range %.3f - %.3f\n",
              tr, nrow(d), median(d$f1), min(d$f1), max(d$f1)))
  cat(sprintf("  configs beating 12_'s per-fold-tuned result: %d of %d\n",
              sum(d$f1 > if (tr == "10_9") 0.668 else 0.840), nrow(d)))
  print(as.data.frame(d %>% arrange(desc(f1)) %>% head(5) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))))
}
