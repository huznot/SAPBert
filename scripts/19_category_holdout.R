# does the model generalise to clinical areas it has never seen
#
# everything so far splits folds by icd-9 code at random which is probably
# optimistic. codes inside a ccs category are clinically adjacent (140, 141, 142
# are all head and neck cancers) so a random split often leaves a near neighbour
# in training
#
# deployment doesnt look like that, mapping a legacy code set means whole areas
# with no examples. so this groups folds by ccs category instead
#
# usage: Rscript 19_category_holdout.R [track]

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "results"
N_OUTER <- 5
N_REPEAT <- 3

CFG <- list(`10_9` = list(K = 10, N = 50, chapter = FALSE, tau = 0.35, rho = 0.85),
            `8_9`  = list(K = 10, N = 10, chapter = FALSE, tau = 0.90, rho = 1.00))

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

# folds by random code, or by whole ccs category
build_folds <- function(map, k, mode) {
  if (mode == "code") {
    cs <- sample(map$ICD_9_CM)
    return(split(cs, rep(seq_len(k), length.out = length(cs))))
  }
  cats <- sample(unique(map$CCS_ID))
  # greedily balance categories across folds by size
  sizes <- map %>% count(CCS_ID)
  cats <- sizes$CCS_ID[order(-sizes$n)]
  load <- rep(0, k); assign <- setNames(integer(length(cats)), cats)
  for (cc in cats) {
    j <- which.min(load); assign[cc] <- j
    load[j] <- load[j] + sizes$n[sizes$CCS_ID == cc]
  }
  lapply(seq_len(k), function(j) map$ICD_9_CM[assign[as.character(map$CCS_ID)] == j])
}

args <- commandArgs(trailingOnly = TRUE)
tracks <- if (length(args)) args else c("10_9", "8_9")
rows <- list()

for (tr in tracks) {
  cfg <- CFG[[tr]]
  fa <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  base <- fa %>% filter(has_truth == 1)
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  feat <- apply_retrieval(base, cfg)
  codes <- sort(unique(base$ICD_9_CM))
  map <- base %>% distinct(ICD_9_CM, CCS_ID)
  FC <- c(setdiff(names(base),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair",
      "CCS_ID","chapter_icd9","chapter_target")), "chapter_rate")

  cat(sprintf("\n######## %s: %d codes, %d CCS categories ########\n",
              tr, length(codes), n_distinct(map$CCS_ID)))

  for (mode in c("code", "category")) {
    f1s <- c(); ps <- c(); rs <- c()
    for (rep in seq_len(N_REPEAT)) {
      set.seed(20260819 + rep)
      folds <- build_folds(map, N_OUTER, mode)
      em <- list()
      for (fi in seq_along(folds)) {
        te_codes <- folds[[fi]]; tr_codes <- setdiff(codes, te_codes)
        if (!length(te_codes) || !length(tr_codes)) next
        trf <- feat %>% filter(ICD_9_CM %in% tr_codes)
        tef <- feat %>% filter(ICD_9_CM %in% te_codes)
        if (!nrow(trf) || !nrow(tef) || sum(trf$y) < 10) next
        ce <- fit_ce(trf)
        a <- add_ce(trf, ce); b <- add_ce(tef, ce)
        b$.p_score <- predict(train_ranker(a, FC), as.matrix(b[, FC]))
        em[[length(em)+1]] <- emit_from_scores(b, cfg$tau, cfg$rho)
      }
      m <- score_pairs(bind_rows(em), truth, codes)
      f1s <- c(f1s, m["f1"]); ps <- c(ps, m["precision"]); rs <- c(rs, m["recall"])
    }
    cat(sprintf("  split by %-8s : P %.3f R %.3f F1 %.4f (sd %.4f over %d repeats)\n",
                mode, mean(ps), mean(rs), mean(f1s), sd(f1s), length(f1s)))
    rows[[length(rows)+1]] <- tibble(track = tr, split = mode,
      precision = mean(ps), recall = mean(rs), f1_mean = mean(f1s),
      f1_sd = sd(f1s), n_repeats = length(f1s))
  }
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, sprintf("category_holdout_%s.csv",
          paste(tracks, collapse = "_"))), row.names = FALSE)

cat("\n===== cost of holding out whole clinical areas =====\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr)
  a <- d$f1_mean[d$split == "code"]; b <- d$f1_mean[d$split == "category"]
  cat(sprintf("  %s: by code %.4f -> by category %.4f  (%+.4f)\n", tr, a, b, b - a))
}
