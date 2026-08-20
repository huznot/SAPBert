# How much would more training data actually buy?
#
# The validation set is small (345 and 302 codes) and that is the limitation
# most often raised about this work. "Small" on its own is not actionable
# though. What matters is whether performance is still climbing at the size we
# have: if it is, collecting more pairs is worth someone's time and this says
# roughly how many; if it has flattened, more data is not the bottleneck and
# effort should go elsewhere.
#
# Holds everything else fixed and varies only the number of training codes.
# Retrieval and the emission rule are fixed at the values 12_ settled on rather
# than re-tuned, so the only thing moving is training set size. Each size is
# repeated over several random subsamples to average out which codes are drawn.
#
# Usage:  Rscript 18_learning_curve.R [track]

source("pipeline_lib.R")
suppressMessages(library(xgboost))
set.seed(20260819)

OUT_DIR <- "../results"
N_OUTER <- 5
N_REPEAT <- 4
SIZES <- c(40, 70, 100, 140, 180, 220, 260)

CFG <- list(`10_9` = list(K = 10, N = 50, chapter = FALSE, tau = 0.35, rho = 0.85),
            `8_9`  = list(K = 10, N = 10, chapter = FALSE, tau = 0.90, rho = 1.00))

prf <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp/(tp+fp) else 0
  r <- if (tp + fn > 0) tp/(tp+fn) else 0
  c(precision = p, recall = r, f1 = if (p + r > 0) 2*p*r/(p+r) else 0)
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
  cfg <- CFG[[tr]]
  fa <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  base <- fa %>% filter(has_truth == 1)
  truth <- base %>% filter(y == 1) %>% select(ICD_9_CM, target)
  feat <- apply_retrieval(base, cfg)
  codes <- sort(unique(base$ICD_9_CM))
  FC <- c(setdiff(names(base),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair",
      "CCS_ID","chapter_icd9","chapter_target")), "chapter_rate")

  cat(sprintf("\n######## %s: %d codes total ########\n", tr, length(codes)))
  set.seed(20260819)
  outer <- make_folds(codes, N_OUTER)

  for (n_train in SIZES) {
    f1s <- c()
    for (rep in seq_len(N_REPEAT)) {
      em <- list()
      for (fi in seq_along(outer)) {
        te_codes <- outer[[fi]]
        pool_codes <- setdiff(codes, te_codes)
        if (n_train > length(pool_codes)) next
        tr_codes <- sample(pool_codes, n_train)
        trf <- feat %>% filter(ICD_9_CM %in% tr_codes)
        if (sum(trf$y) < 10) next
        ce <- fit_ce(trf)
        a <- add_ce(trf, ce)
        b <- add_ce(feat %>% filter(ICD_9_CM %in% te_codes), ce)
        b$.p_score <- predict(train_ranker(a, FC), as.matrix(b[, FC]))
        em[[length(em)+1]] <- emit_from_scores(b, cfg$tau, cfg$rho)
      }
      if (!length(em)) next
      f1s <- c(f1s, score_pairs(bind_rows(em), truth, codes)["f1"])
    }
    cat(sprintf("  n_train %3d : F1 %.4f (sd %.4f over %d repeats)\n",
                n_train, mean(f1s), sd(f1s), length(f1s)))
    rows[[length(rows)+1]] <- tibble(track = tr, n_train = n_train,
      f1_mean = mean(f1s), f1_sd = sd(f1s), n_repeats = length(f1s))
  }
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, sprintf("learning_curve_%s.csv",
          paste(tracks, collapse = "_"))), row.names = FALSE)

cat("\n===== is it still climbing? =====\n")
for (tr in unique(res$track)) {
  d <- res %>% filter(track == tr) %>% arrange(n_train)
  # slope over the last third of the curve, per 100 extra training codes
  tail_d <- tail(d, 3)
  slope <- coef(lm(f1_mean ~ n_train, data = tail_d))[2] * 100
  cat(sprintf("\n%s: F1 %.3f at n=%d, %.3f at n=%d\n", tr,
              d$f1_mean[1], d$n_train[1], tail(d$f1_mean,1), tail(d$n_train,1)))
  cat(sprintf("  slope over the last three points: %+.4f F1 per +100 training codes\n", slope))
  # log fit, for a rough sense of where it goes. Extrapolation beyond the
  # observed range is a guess, not a measurement.
  fit <- lm(f1_mean ~ log(n_train), data = d)
  for (n in c(500, 1000, 2000)) {
    cat(sprintf("  log-fit projection at n=%-5d: %.3f\n", n,
                predict(fit, newdata = data.frame(n_train = n))))
  }
}
cat("\nProjections assume the curve keeps its shape and are not measurements.\n")
