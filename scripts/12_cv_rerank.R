# trains the reranker and scores it with grouped cross validation
# 5 folds grouped by icd-9 code so a code and all its candidates stay together,
# splitting on pairs would leak some of a codes targets into training
# everything is fitted inside the training folds, test folds are only predicted on
# the old pipeline runs through the same folds so the comparison is fair
#
# usage: Rscript 12_cv_rerank.R          both tracks
#        Rscript 12_cv_rerank.R 10_9     one track, run both in parallel then 12b_

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(xgboost))

set.seed(20260819)

ORIG_BASE    <- "data/original"
SAPBERT_BASE <- "data/sapbert"
GEN_BASE     <- "data/generated"
OUT_DIR      <- "results"
N_OUTER <- 5
N_INNER <- 3

# abstain thresholds the inner cv picks between. 0 emits for every code, which
# is what the system did before, so it stays in the grid as the null option
FLOORS <- c(0, 0.15, 0.3, 0.45, 0.6)

VAL_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

TRACKS <- list(
  `10_9` = list(tcn = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                find_fn = find_icd10ca_chapter, align = chapter_alignment_10,
                manual_target_col = "ICD-10-CA", manual_sheet = "Validation_ICD9_ICD10",
                excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                sim = list(SapBERT = file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                           ClinicalBERT = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))),
  `8_9`  = list(tcn = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                find_fn = find_icda8_chapter, align = chapter_alignment_8,
                manual_target_col = "ICDA-8", manual_sheet = "Validaion_ICD9_ICD8",
                excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                sim = list(SapBERT = file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx"),
                           ClinicalBERT = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")))
)

prf <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp/(tp+fp) else 0
  r <- if (tp + fn > 0) tp/(tp+fn) else 0
  c(precision = p, recall = r,
    f1 = if (p + r > 0) 2*p*r/(p+r) else 0,
    accuracy = if (tp+fp+fn > 0) tp/(tp+fp+fn) else 0)
}

# score emitted pairs against the truth for a given set of codes
score_pairs <- function(emitted, truth, codes) {
  e <- emitted %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  t <- truth   %>% filter(ICD_9_CM %in% codes) %>% distinct(ICD_9_CM, target)
  tp <- nrow(dplyr::semi_join(e, t, by = c("ICD_9_CM","target")))
  prf(tp, nrow(e) - tp, nrow(t) - tp)
}

make_folds <- function(codes, k) {
  codes <- sample(codes)
  split(codes, rep(seq_len(k), length.out = length(codes)))
}

all_rows <- list(); fold_rows <- list(); imp_rows <- list(); pred_rows <- list()

# one track at a time so both can run in parallel, 12b_ merges them
.args <- commandArgs(trailingOnly = TRUE)
SELECTED_TRACKS <- if (length(.args)) .args else names(TRACKS)
stopifnot(all(SELECTED_TRACKS %in% names(TRACKS)))

for (tr in SELECTED_TRACKS) {
  tk <- TRACKS[[tr]]
  cat(sprintf("\n################ TRACK %s ################\n", tr))

  feat_all <- readRDS(file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  n_true_pairs <- attr(feat_all, "n_true_pairs")

  # only codes with a known answer can be cross validated. the rest go to 14_,
  # leaving them in would look like all negative rows and drag the model down
  feat <- feat_all %>% filter(has_truth == 1)
  truth <- feat %>% filter(y == 1) %>% select(ICD_9_CM, target)
  codes <- sort(unique(feat$ICD_9_CM))

  # truth restricted to pairs the candidate pool can reach, which is what the
  # recall denominator has always been here. a pair the pool never retrieved is
  # a miss either way, so the full crosswalk is scored alongside it and both are
  # reported. the pool one is the favourable choice and should not stand alone
  manual_full <- read_excel(VAL_XLSX, sheet = tk$manual_sheet)
  excl_full   <- as.character(read_excel(VAL_XLSX, sheet = tk$excl_sheet)$`ICD-9-CM`)
  # no filter on excl_full here. those codes have no rows in the manual
  # crosswalk to begin with, so they add nothing to the truth set, but they are
  # in codes and any pair emitted for them counts against precision
  truth_full <- manual_full %>%
    transmute(ICD_9_CM = as.character(`ICD-9-CM`),
              target = as.character(.data[[tk$manual_target_col]])) %>%
    filter(ICD_9_CM %in% codes) %>%
    distinct()
  cat(sprintf("truth pairs: %d reachable in the pool, %d in the full crosswalk\n",
              nrow(truth), nrow(truth_full)))
  cat(sprintf("%d candidates, %d codes, %d positives (pool ceiling %.4f of %d true pairs)\n",
              nrow(feat), length(codes), sum(feat$y), sum(feat$y)/n_true_pairs, n_true_pairs))

  FEATURE_COLS <- setdiff(names(feat),
    c("ICD_9_CM","target","y","has_truth","icd9_label","target_label","chapter_pair","CCS_ID",
      "chapter_icd9","chapter_target"))

  # baseline emissions per thr/top_n/flag. same every run so cache them,
  # rm results/base_emit_<track>.rds to rebuild
  base_cache <- file.path(OUT_DIR, sprintf("base_emit_%s.rds", tr))
  if (file.exists(base_cache)) {
    cat("Loading cached baseline emissions...\n")
    base_emit <- readRDS(base_cache)
  } else {
    cat("Precomputing baseline pipeline emissions (first run only)...\n")
    cooc_df <- load_cooccurrence_df(file.path(ORIG_BASE, tk$cooc_file))
    base_emit <- list()
    for (mdl in names(tk$sim)) {
      sheets <- load_similarity_sheets(tk$sim[[mdl]])
      for (thr in c(0.95, 0.99, 0.995, 0.999)) {
        sim_df <- get_similarity_scores_from_sheets(sheets, thr, tk$tcn)
        for (tn in c(3,5,10,15,20,25,30)) {
          cdf <- get_cooccurrence_codes_from_df(cooc_df, tn, tk$icd9_col, tk$target_col, tk$tcn)
          merged <- merge_and_flag(sim_df, cdf, tk$tcn, tk$find_fn, tk$align)
          for (fc in 1:4) {
            key <- sprintf("%s|%s|%s|%s", mdl, thr, tn, fc)
            base_emit[[key]] <- select_rows_by_flags(merged, fc) %>%
              transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(.data[[tk$tcn]])) %>%
              distinct()
          }
        }
      }
      cat(sprintf("  %s done\n", mdl))
    }
    saveRDS(base_emit, base_cache)
  }
  base_keys <- names(base_emit)

  # retrieval configs the inner cv picks between
  RETRIEVAL <- list()
  for (K in c(10, 25, 50)) for (N in c(10, 30, 50)) for (chap in c(TRUE, FALSE))
    RETRIEVAL[[length(RETRIEVAL)+1]] <- list(K = K, N = N, chapter = chap)

  apply_retrieval <- function(df, cfg) {
    rank_cols <- grep("^simrank_", names(df), value = TRUE)
    keep <- Reduce(`|`, lapply(rank_cols, function(c) df[[c]] <= cfg$K)) | (df$cooc_rank <= cfg$N)
    if (cfg$chapter) keep <- keep & (df$chapter_ok == 1)
    df[keep, , drop = FALSE]
  }

  # learn how often each chapter pair is real, from train only. the hand written
  # table misses cases, icd-9 ch 3 is "endocrine ... and immunity" but never maps
  # to D80-D89. patching it by hand would just be fitting the answers
  fit_chapter_encoding <- function(train_df) {
    prior <- mean(train_df$y)
    enc <- train_df %>% group_by(chapter_pair) %>%
      summarise(n = n(), pos = sum(y), .groups = "drop") %>%
      mutate(chapter_rate = (pos + 20 * prior) / (n + 20))
    list(enc = enc %>% select(chapter_pair, chapter_rate), prior = prior)
  }
  apply_chapter_encoding <- function(df, ce) {
    df %>% left_join(ce$enc, by = "chapter_pair") %>%
      mutate(chapter_rate = ifelse(is.na(chapter_rate), ce$prior, chapter_rate))
  }

  train_ranker <- function(train_df, feature_cols) {
    dtrain <- xgb.DMatrix(as.matrix(train_df[, feature_cols]), label = train_df$y)
    xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss",
                            max_depth = 6, eta = 0.1, subsample = 0.8,
                            colsample_bytree = 0.8, min_child_weight = 5,
                            nthread = parallel::detectCores()),
              data = dtrain, nrounds = 200, verbose = 0)
  }

  # emit if confident on its own (tau) or nearly as good as the best one for that
  # code (rho). most codes have several valid targets so one global threshold does
  # badly. rho = 1 means top 1 only.
  #
  # floor lets a code emit nothing. 61 of the 354 codes have no correct answer,
  # and without this the best candidate is always emitted, so those codes were
  # guaranteed a false positive. a code is skipped entirely when its best score
  # is under floor. floor = 0 is the old behaviour and is in the grid, so the
  # inner cv keeps abstention only if it earns its place
  emit_from_scores <- function(df, tau, rho, floor = 0) {
    df %>% group_by(ICD_9_CM) %>%
      mutate(.mx = max(.p_score)) %>% ungroup() %>%
      filter(.mx >= floor, .p_score >= tau | .p_score >= rho * .mx) %>%
      select(ICD_9_CM, target)
  }

  outer_folds <- make_folds(codes, N_OUTER)
  oof_pred <- list()
  fold_base_key <- character(N_OUTER)   # baseline config chosen per fold

  for (fi in seq_along(outer_folds)) {
    test_codes  <- outer_folds[[fi]]
    train_codes <- setdiff(codes, test_codes)
    cat(sprintf("\n--- outer fold %d/%d (%d train codes, %d test codes) ---\n",
                fi, N_OUTER, length(train_codes), length(test_codes)))

    train_full <- feat %>% filter(ICD_9_CM %in% train_codes)

    # inner cv picks the retrieval config and decision rule
    #
    # train once per inner fold on the full pool then test every config by
    # filtering the scored rows. retrieval only decides what is eligible, not what
    # a candidate looks like, so retraining per config was wasted work. 54 fits
    # per outer fold down to 3.
    #
    # note this changes the procedure not just the speed, the ranker now trains on
    # the full pool and retrieval becomes an inference time filter
    inner_folds <- make_folds(train_codes, N_INNER)
    fcols <- c(FEATURE_COLS, "chapter_rate")

    inner_scored <- list()
    for (ii in seq_along(inner_folds)) {
      va_codes <- inner_folds[[ii]]
      tr_codes <- setdiff(train_codes, va_codes)
      tr_df <- train_full %>% filter(ICD_9_CM %in% tr_codes)
      va_df <- train_full %>% filter(ICD_9_CM %in% va_codes)
      if (nrow(tr_df) == 0 || nrow(va_df) == 0 || sum(tr_df$y) == 0) next
      ce <- fit_chapter_encoding(tr_df)
      tr_df <- apply_chapter_encoding(tr_df, ce); va_df <- apply_chapter_encoding(va_df, ce)
      bst_i <- train_ranker(tr_df, fcols)
      va_df$.p_score <- predict(bst_i, as.matrix(va_df[, fcols]))
      inner_scored[[length(inner_scored)+1]] <- va_df
    }
    isc_all <- bind_rows(inner_scored)

    best <- list(f1 = -1)
    for (cfg in RETRIEVAL) {
      isc <- apply_retrieval(isc_all, cfg)
      if (nrow(isc) == 0) next
      for (tau in seq(0.05, 0.95, by = 0.05)) {
        for (rho in c(0.3, 0.5, 0.7, 0.85, 1.0)) {
          for (fl in FLOORS) {
            em <- emit_from_scores(isc, tau, rho, fl)
            m <- score_pairs(em, truth, train_codes)
            if (m["f1"] > best$f1)
              best <- list(f1 = m["f1"], cfg = cfg, tau = tau, rho = rho, floor = fl)
          }
        }
      }
    }
    cat(sprintf("  inner-CV pick: K=%d N=%d chapter=%s tau=%.2f rho=%.2f floor=%.2f (inner F1 %.3f)\n",
                best$cfg$K, best$cfg$N, best$cfg$chapter, best$tau, best$rho, best$floor, best$f1))

    # fit on the full training set, predict the test fold
    tr_df <- train_full
    te_df <- feat %>% filter(ICD_9_CM %in% test_codes)
    ce <- fit_chapter_encoding(tr_df)
    tr_df <- apply_chapter_encoding(tr_df, ce); te_df <- apply_chapter_encoding(te_df, ce)
    bst <- train_ranker(tr_df, fcols)
    te_df$.p_score <- predict(bst, as.matrix(te_df[, fcols]))
    # retrieval applied at inference as an eligibility filter
    te_df <- apply_retrieval(te_df, best$cfg)

    em <- emit_from_scores(te_df, best$tau, best$rho, best$floor)
    m_rr <- score_pairs(em, truth, test_codes)
    oof_pred[[fi]] <- te_df %>% mutate(.tau = best$tau, .rho = best$rho, .floor = best$floor, fold = fi)

    # baseline, same protocol
    best_key <- NULL; best_bf1 <- -1
    for (k in base_keys) {
      mm <- score_pairs(base_emit[[k]], truth, train_codes)
      if (mm["f1"] > best_bf1) { best_bf1 <- mm["f1"]; best_key <- k }
    }
    fold_base_key[fi] <- best_key
    m_bl <- score_pairs(base_emit[[best_key]], truth, test_codes)
    cat(sprintf("  baseline pick on train: %s (train F1 %.3f)\n", best_key, best_bf1))
    cat(sprintf("  TEST FOLD  baseline F1 %.3f | reranker F1 %.3f\n", m_bl["f1"], m_rr["f1"]))

    fold_rows[[length(fold_rows)+1]] <- tibble(
      track = tr, fold = fi, n_test_codes = length(test_codes),
      baseline_f1 = m_bl["f1"], baseline_precision = m_bl["precision"],
      baseline_recall = m_bl["recall"], baseline_accuracy = m_bl["accuracy"],
      rerank_f1 = m_rr["f1"], rerank_precision = m_rr["precision"],
      rerank_recall = m_rr["recall"], rerank_accuracy = m_rr["accuracy"],
      K = best$cfg$K, N = best$cfg$N, chapter = best$cfg$chapter, tau = best$tau, rho = best$rho, floor = best$floor,
      baseline_config = best_key)

    imp <- xgb.importance(model = bst)
    imp_rows[[length(imp_rows)+1]] <- as_tibble(imp) %>% mutate(track = tr, fold = fi)
  }

  # checkpoint before aggregating so a bug down there doesnt bin the whole run
  oof <- bind_rows(oof_pred)
  saveRDS(list(oof = oof, fold_base_key = fold_base_key, outer_folds = outer_folds),
          file.path(OUT_DIR, sprintf("cv_rerank_checkpoint_%s.rds", tr)))

  # pooled out-of-fold results
  em_all <- bind_rows(lapply(oof_pred, function(d) emit_from_scores(d, d$.tau[1], d$.rho[1], d$.floor[1])))
  m_rr_all <- score_pairs(em_all, truth, codes)

  # baseline pooled the same way, each folds test codes scored with the config
  # picked on that folds training codes
  stopifnot(all(nzchar(fold_base_key)), all(fold_base_key %in% names(base_emit)))
  bl_pairs <- bind_rows(lapply(seq_along(outer_folds), function(i) {
    base_emit[[fold_base_key[i]]] %>% filter(ICD_9_CM %in% outer_folds[[i]])
  }))
  m_bl_all <- score_pairs(bl_pairs, truth, codes)

  # in sample baseline, best config scored on everything
  # index by position. score_pairs returns a named vector so sapply names come out
  # as "SapBERT|0.95|3|1.f1" and the lookup silently returns NULL
  in_sample <- vapply(base_keys,
                      function(k) unname(score_pairs(base_emit[[k]], truth, codes)["f1"]),
                      numeric(1))
  m_bl_insample <- score_pairs(base_emit[[base_keys[which.max(in_sample)]]], truth, codes)

  cat(sprintf("\n===== TRACK %s POOLED RESULTS =====\n", tr))
  cat(sprintf("  baseline, IN-SAMPLE (old method)   : P %.3f R %.3f F1 %.3f Acc %.3f\n",
              m_bl_insample["precision"], m_bl_insample["recall"], m_bl_insample["f1"], m_bl_insample["accuracy"]))
  cat(sprintf("  baseline, held-out CV              : P %.3f R %.3f F1 %.3f Acc %.3f\n",
              m_bl_all["precision"], m_bl_all["recall"], m_bl_all["f1"], m_bl_all["accuracy"]))
  cat(sprintf("  RERANKER, held-out CV              : P %.3f R %.3f F1 %.3f Acc %.3f\n",
              m_rr_all["precision"], m_rr_all["recall"], m_rr_all["f1"], m_rr_all["accuracy"]))

  # the same two held-out outputs scored against the whole crosswalk
  m_bl_full <- score_pairs(bl_pairs, truth_full, codes)
  m_rr_full <- score_pairs(em_all,   truth_full, codes)
  cat(sprintf("  baseline, full crosswalk denominator: P %.3f R %.3f F1 %.3f Acc %.3f\n",
              m_bl_full["precision"], m_bl_full["recall"], m_bl_full["f1"], m_bl_full["accuracy"]))
  cat(sprintf("  RERANKER, full crosswalk denominator: P %.3f R %.3f F1 %.3f Acc %.3f\n",
              m_rr_full["precision"], m_rr_full["recall"], m_rr_full["f1"], m_rr_full["accuracy"]))

  all_rows[[length(all_rows)+1]] <- bind_rows(
    tibble(track = tr, method = "baseline_in_sample",  t(m_bl_insample)),
    tibble(track = tr, method = "baseline_heldout_cv", t(m_bl_all)),
    tibble(track = tr, method = "rerank_heldout_cv",   t(m_rr_all)),
    tibble(track = tr, method = "baseline_heldout_cv_fullden", t(m_bl_full)),
    tibble(track = tr, method = "rerank_heldout_cv_fullden",   t(m_rr_full)))
  pred_rows[[length(pred_rows)+1]] <- oof %>%
    select(ICD_9_CM, target, y, .p_score, fold) %>% mutate(track = tr)
}

results <- bind_rows(all_rows)

# per track files so parallel runs dont clobber each other
for (tr in SELECTED_TRACKS) {
  saveRDS(list(results = results %>% filter(track == tr),
               folds   = bind_rows(fold_rows) %>% filter(track == tr),
               imp     = bind_rows(imp_rows) %>% filter(track == tr),
               preds   = bind_rows(pred_rows) %>% filter(track == tr)),
          file.path(OUT_DIR, sprintf("cv_rerank_part_%s.rds", tr)))
}

write.csv(results, file.path(OUT_DIR, "cv_rerank_results.csv"), row.names = FALSE)
write.csv(bind_rows(fold_rows), file.path(OUT_DIR, "cv_rerank_folds.csv"), row.names = FALSE)
saveRDS(bind_rows(pred_rows), file.path(OUT_DIR, "cv_rerank_predictions.rds"))

importance <- bind_rows(imp_rows) %>% group_by(track, Feature) %>%
  summarise(gain = mean(Gain), .groups = "drop") %>% arrange(track, desc(gain))
write.csv(importance, file.path(OUT_DIR, "cv_rerank_importance.csv"), row.names = FALSE)

cat("\n\n================ FINAL ================\n")
print(as.data.frame(results %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("\n--- top features by mean gain ---\n")
print(as.data.frame(importance %>% group_by(track) %>% slice_head(n = 12) %>% ungroup()))
cat("\nWritten to results/cv_rerank_*.csv\n")
