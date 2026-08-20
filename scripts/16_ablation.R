# ablation, which features earn their place and whether a ranking objective
# beats scoring candidates one at a time
#
# same folds and same retrieval as 12_ so the only thing that changes is the
# thing being tested. retrieval is fixed at what 12_ picked most often rather
# than re-searched, which keeps it clean and short
#
# two questions. which feature groups matter, since the mutual features and the
# training change shipped together and couldnt be attributed. and whether
# rank:pairwise / rank:ndcg beat binary:logistic, since the task is picking a
# set per code and most codes have several correct targets
#
# usage: Rscript 16_ablation.R [track]

source("pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "../results"
N_OUTER <- 5
N_INNER <- 3

RETRIEVAL <- list(`10_9` = list(K = 10, N = 50, chapter = FALSE),
                  `8_9`  = list(K = 10, N = 10, chapter = FALSE))

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

# feature groups by name pattern
GROUPS <- list(
  mutual  = "^simrankrev_|^simrelrev_|^ens_mean_relrev|^ens_best_rankrev|^ens_rrf_rev|^ens_mutual|^is_mutual_top1|^ens_direction_gap|^target_n_suitors|^target_rank_here|^target_is_best_here",
  lexical = "^lex_",
  cooc    = "^cooc_|^has_cooc",
  chapter = "^chapter_",
  ensemble= "^ens_",
  mpnet   = "_mpnet$",
  clinicalbert = "_clinicalbert$"
)

VARIANTS <- list(
  list(name = "full",              drop = NULL,          objective = "binary:logistic"),
  list(name = "no_mutual",         drop = "mutual",      objective = "binary:logistic"),
  list(name = "no_lexical",        drop = "lexical",     objective = "binary:logistic"),
  list(name = "no_cooccurrence",   drop = "cooc",        objective = "binary:logistic"),
  list(name = "no_chapter",        drop = "chapter",     objective = "binary:logistic"),
  list(name = "no_ensemble",       drop = "ensemble",    objective = "binary:logistic"),
  list(name = "no_mpnet",          drop = "mpnet",       objective = "binary:logistic"),
  list(name = "no_clinicalbert",   drop = "clinicalbert",objective = "binary:logistic"),
  list(name = "rank_pairwise",     drop = NULL,          objective = "rank:pairwise"),
  list(name = "rank_ndcg",         drop = NULL,          objective = "rank:ndcg")
)

train_ranker <- function(train_df, fcols, objective) {
  # ranking objectives need rows grouped by code, contiguous, with group sizes
  if (grepl("^rank:", objective)) {
    train_df <- train_df %>% arrange(ICD_9_CM)
    grp <- train_df %>% count(ICD_9_CM) %>% pull(n)
  }
  d <- xgb.DMatrix(as.matrix(train_df[, fcols]), label = train_df$y)
  if (grepl("^rank:", objective)) setinfo(d, "group", grp)
  xgb.train(params = list(objective = objective,
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

apply_retrieval <- function(df, cfg) {
  rank_cols <- grep("^simrank_", names(df), value = TRUE)
  keep <- Reduce(`|`, lapply(rank_cols, function(c) df[[c]] <= cfg$K)) | (df$cooc_rank <= cfg$N)
  if (cfg$chapter) keep <- keep & (df$chapter_ok == 1)
  df[keep, , drop = FALSE]
}

args <- commandArgs(trailingOnly = TRUE)
tracks <- if (length(args)) args else c("10_9", "8_9")
rows <- list()

for (tr in tracks) {
  cat(sprintf("\n################ TRACK %s ################\n", tr))
  feat_all <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  base <- feat_all %>% filter(has_truth == 1)
  # truth comes from the FULL pool, before retrieval. Taking it after would
  # drop the positives retrieval missed out of the recall denominator and
  # inflate F1, which is not comparable to what 12_ reports.
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  feat <- base %>% apply_retrieval(RETRIEVAL[[tr]])
  codes <- sort(unique(feat$ICD_9_CM))
  cat(sprintf("%d candidates, %d codes, %d positives\n", nrow(feat), length(codes), sum(feat$y)))

  ALL_FEATS <- setdiff(names(feat),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair","CCS_ID",
      "chapter_icd9","chapter_target"))

  set.seed(20260819)
  outer_folds <- make_folds(codes, N_OUTER)

  for (v in VARIANTS) {
    fcols <- ALL_FEATS
    if (!is.null(v$drop)) fcols <- setdiff(fcols, grep(GROUPS[[v$drop]], fcols, value = TRUE))
    fcols <- c(fcols, "chapter_rate")
    if (identical(v$drop, "chapter")) fcols <- setdiff(fcols, "chapter_rate")

    emitted_all <- list()
    for (fi in seq_along(outer_folds)) {
      test_codes  <- outer_folds[[fi]]
      train_codes <- setdiff(codes, test_codes)
      train_full  <- feat %>% filter(ICD_9_CM %in% train_codes)

      # inner cv only to pick tau/rho
      inner <- make_folds(train_codes, N_INNER)
      isc <- list()
      for (ii in seq_along(inner)) {
        va <- inner[[ii]]; trn <- setdiff(train_codes, va)
        a <- train_full %>% filter(ICD_9_CM %in% trn)
        b <- train_full %>% filter(ICD_9_CM %in% va)
        if (!nrow(a) || !nrow(b) || sum(a$y) == 0) next
        ce <- fit_chapter_encoding(a)
        a <- apply_chapter_encoding(a, ce); b <- apply_chapter_encoding(b, ce)
        m <- train_ranker(a, fcols, v$objective)
        b$.p_score <- predict(m, as.matrix(b[, fcols]))
        isc[[length(isc)+1]] <- b
      }
      isc <- bind_rows(isc)
      best <- list(f1 = -1, tau = 0.5, rho = 0.7)
      # ranking objectives give unbounded scores so sweep quantiles of the score
      # distribution rather than fixed probabilities
      taus <- if (grepl("^rank:", v$objective))
                unname(quantile(isc$.p_score, seq(0.5, 0.995, length.out = 19))) else
                seq(0.05, 0.95, by = 0.05)
      for (tau in taus) for (rho in c(0.3, 0.5, 0.7, 0.85, 1.0)) {
        m <- score_pairs(emit_from_scores(isc, tau, rho), truth, train_codes)
        if (m["f1"] > best$f1) best <- list(f1 = m["f1"], tau = tau, rho = rho)
      }

      ce <- fit_chapter_encoding(train_full)
      trn <- apply_chapter_encoding(train_full, ce)
      te  <- apply_chapter_encoding(feat %>% filter(ICD_9_CM %in% test_codes), ce)
      m <- train_ranker(trn, fcols, v$objective)
      te$.p_score <- predict(m, as.matrix(te[, fcols]))
      emitted_all[[fi]] <- emit_from_scores(te, best$tau, best$rho)
    }

    mm <- score_pairs(bind_rows(emitted_all), truth, codes)
    cat(sprintf("  %-18s P %.3f R %.3f F1 %.3f Acc %.3f  (%d features)\n",
                v$name, mm["precision"], mm["recall"], mm["f1"], mm["accuracy"], length(fcols)))
    rows[[length(rows)+1]] <- tibble(track = tr, variant = v$name,
      objective = v$objective, n_features = length(fcols),
      precision = mm["precision"], recall = mm["recall"],
      f1 = mm["f1"], accuracy = mm["accuracy"])
  }
}

res <- bind_rows(rows)
out <- file.path(OUT_DIR, sprintf("ablation_%s.csv", paste(tracks, collapse = "_")))
write.csv(res, out, row.names = FALSE)

cat("\n\n===== ablation, held-out =====\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr)
  full <- d$f1[d$variant == "full"]
  cat(sprintf("\n-- %s (full = %.3f) --\n", tr, full))
  print(as.data.frame(d %>% mutate(delta_vs_full = round(f1 - full, 4)) %>%
    select(variant, n_features, precision, recall, f1, delta_vs_full) %>%
    arrange(desc(f1)) %>% mutate(across(where(is.numeric), ~round(.x, 4)))))
}
cat(sprintf("\nWritten to %s\n", out))
