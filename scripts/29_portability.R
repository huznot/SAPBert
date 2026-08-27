# what does this method give someone who does not have everything we have
#
# 16_ablation.R only drops a feature group from the scoring model. the candidate
# list still keeps anything with cooc_rank <= N, so its no_cooccurrence row is
# not a "no co-occurrence data" number. this script drops co-occurrence from
# retrieval as well, which is the situation of a group that has code labels but
# no linked health records.
#
# three settings, same folds and same scoring model each time:
#   full          everything we have
#   no_cooc       code labels and chapter table only, no health records
#   labels_only   code labels only, no health records and no chapter table
#
# usage: Rscript 29_portability.R [track]

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "results"
N_OUTER <- 5
N_INNER <- 3

RETRIEVAL <- list(`10_9` = list(K = 10, N = 50, chapter = FALSE),
                  `8_9`  = list(K = 10, N = 10, chapter = FALSE))

# use_cooc FALSE means retrieval ignores the co-occurrence list entirely, so the
# candidates come from the embedding models alone
SETTINGS <- list(
  list(name = "full",        use_cooc = TRUE,  drop = character(0)),
  list(name = "no_cooc",     use_cooc = FALSE, drop = c("cooc")),
  list(name = "labels_only", use_cooc = FALSE, drop = c("cooc", "chapter"))
)

GROUPS <- list(cooc    = "^cooc_|^has_cooc",
               chapter = "^chapter_")

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

# the co-occurrence branch of the keep rule is dropped when use_cooc is FALSE
apply_retrieval <- function(df, cfg, use_cooc) {
  rank_cols <- grep("^simrank_", names(df), value = TRUE)
  keep <- Reduce(`|`, lapply(rank_cols, function(c) df[[c]] <= cfg$K))
  if (use_cooc) keep <- keep | (df$cooc_rank <= cfg$N)
  df[keep, , drop = FALSE]
}

args <- commandArgs(trailingOnly = TRUE)
tracks <- if (length(args)) args else c("10_9", "8_9")
rows <- list()

for (tr in tracks) {
  cat(sprintf("\n################ TRACK %s ################\n", tr))
  feat_all <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  base <- feat_all %>% filter(has_truth == 1)
  # truth comes from the full pool, before retrieval, so a positive that
  # retrieval missed still counts against recall
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  codes <- sort(unique(base$ICD_9_CM))

  for (st in SETTINGS) {
    feat <- base %>% apply_retrieval(RETRIEVAL[[tr]], st$use_cooc)
    ALL_FEATS <- setdiff(names(feat),
      c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair","CCS_ID",
        "chapter_icd9","chapter_target"))
    fcols <- ALL_FEATS
    for (g in st$drop) fcols <- setdiff(fcols, grep(GROUPS[[g]], fcols, value = TRUE))
    use_chapter <- !("chapter" %in% st$drop)
    if (use_chapter) fcols <- c(fcols, "chapter_rate")

    set.seed(20260819)
    outer_folds <- make_folds(codes, N_OUTER)
    emitted_all <- list()
    for (fi in seq_along(outer_folds)) {
      test_codes  <- outer_folds[[fi]]
      train_codes <- setdiff(codes, test_codes)
      train_full  <- feat %>% filter(ICD_9_CM %in% train_codes)

      inner <- make_folds(train_codes, N_INNER)
      isc <- list()
      for (ii in seq_along(inner)) {
        va <- inner[[ii]]; trn <- setdiff(train_codes, va)
        a <- train_full %>% filter(ICD_9_CM %in% trn)
        b <- train_full %>% filter(ICD_9_CM %in% va)
        if (!nrow(a) || !nrow(b) || sum(a$y) == 0) next
        if (use_chapter) {
          ce <- fit_chapter_encoding(a)
          a <- apply_chapter_encoding(a, ce); b <- apply_chapter_encoding(b, ce)
        }
        m <- train_ranker(a, fcols)
        b$.p_score <- predict(m, as.matrix(b[, fcols]))
        isc[[length(isc)+1]] <- b
      }
      isc <- bind_rows(isc)
      best <- list(f1 = -1, tau = 0.5, rho = 0.7)
      for (tau in seq(0.05, 0.95, by = 0.05)) for (rho in c(0.3, 0.5, 0.7, 0.85, 1.0)) {
        m <- score_pairs(emit_from_scores(isc, tau, rho), truth, train_codes)
        if (m["f1"] > best$f1) best <- list(f1 = m["f1"], tau = tau, rho = rho)
      }

      trn <- train_full
      te  <- feat %>% filter(ICD_9_CM %in% test_codes)
      if (use_chapter) {
        ce <- fit_chapter_encoding(train_full)
        trn <- apply_chapter_encoding(trn, ce); te <- apply_chapter_encoding(te, ce)
      }
      m <- train_ranker(trn, fcols)
      te$.p_score <- predict(m, as.matrix(te[, fcols]))
      emitted_all[[fi]] <- emit_from_scores(te, best$tau, best$rho)
    }

    mm <- score_pairs(bind_rows(emitted_all), truth, codes)
    cat(sprintf("  %-13s P %.3f R %.3f F1 %.3f Acc %.3f  (%d candidates, %d features)\n",
                st$name, mm["precision"], mm["recall"], mm["f1"], mm["accuracy"],
                nrow(feat), length(fcols)))
    rows[[length(rows)+1]] <- tibble(track = tr, setting = st$name,
      n_candidates = nrow(feat), n_features = length(fcols),
      precision = mm["precision"], recall = mm["recall"],
      f1 = mm["f1"], accuracy = mm["accuracy"])
  }
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, "portability.csv"), row.names = FALSE)

cat("\n\n===== what each setting needs, held out =====\n")
cat("  full         code labels + linked health records + chapter table\n")
cat("  no_cooc      code labels + chapter table\n")
cat("  labels_only  code labels only\n\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr)
  fl <- d$f1[d$setting == "full"]
  cat(sprintf("-- %s (full = %.3f) --\n", tr, fl))
  print(as.data.frame(d %>% mutate(delta_vs_full = round(f1 - fl, 4)) %>%
    select(setting, n_candidates, n_features, precision, recall, f1, delta_vs_full) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))))
  cat("\n")
}
