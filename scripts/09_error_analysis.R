# works out where the pipeline loses points
#
# two stages fail differently. candidate generation decides what gets considered
# at all and anything it drops is gone for good. selection then picks what to
# emit and whatever it gets wrong costs precision
#
# measures both plus an oracle, the f1 a perfect selector would get on the
# current pool. thats the hard ceiling for any amount of reranking
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

ORIG_BASE    <- "data/original"
SAPBERT_BASE <- "data/sapbert"
GEN_BASE     <- "data/generated"
OUT_DIR      <- "results"

ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                     sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
VAL_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

TRACKS <- list(
  `10_9` = list(target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3",
                target_col = "ICD_10_CA_Code3", find_fn = find_icd10ca_chapter,
                align = chapter_alignment_10, manual_target_col = "ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                best = list(thr = 0.95, tn = 30, flag = 4)),
  `8_9`  = list(target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code",
                target_col = "ICDA_8_Code", find_fn = find_icda8_chapter,
                align = chapter_alignment_8, manual_target_col = "ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                best = list(thr = 0.99, tn = 3, flag = 2))
)
for (tr in names(TRACKS)) {
  TRACKS[[tr]]$manual <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$manual_sheet)
  TRACKS[[tr]]$excl   <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$excl_sheet)
  TRACKS[[tr]]$cooc   <- load_cooccurrence_df(file.path(ORIG_BASE, TRACKS[[tr]]$cooc_file))
}

MODELS <- list(
  SapBERT = list(`10_9` = file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                 `8_9`  = file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx")),
  mpnet   = list(`10_9` = file.path(GEN_BASE, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"),
                 `8_9`  = file.path(GEN_BASE, "cosine_similarity_matrices_8_9_mpnet_base.xlsx"))
)

# true pairs, after the same exclusion list the scorer uses
true_pairs <- function(tk) {
  excluded <- as.character(tk$excl$`ICD-9-CM`)
  tk$manual %>%
    transmute(ICD_9_CM = as.character(`ICD-9-CM`),
              target   = as.character(.data[[tk$manual_target_col]])) %>%
    filter(!(ICD_9_CM %in% excluded)) %>%
    distinct()
}

f1_of <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp / (tp + fp) else 0
  r <- if (tp + fn > 0) tp / (tp + fn) else 0
  list(precision = p, recall = r,
       f1 = if (p + r > 0) 2 * p * r / (p + r) else 0,
       accuracy = if (tp + fp + fn > 0) tp / (tp + fp + fn) else 0)
}

rows <- list()
loss_rows <- list()

for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  tcn <- tk$target_col_name
  tp_df <- true_pairs(tk)
  cat(sprintf("\n########## track %s : %d true pairs, %d distinct ICD-9 codes ##########\n",
              tr, nrow(tp_df), n_distinct(tp_df$ICD_9_CM)))

  for (mdl in names(MODELS)) {
    b <- tk$best
    sheets <- load_similarity_sheets(MODELS[[mdl]][[tr]])

    # --- stage 1: raw signals, before any filtering -------------------
    sim_raw <- get_similarity_scores_from_sheets(sheets, b$thr, tcn) %>%
      transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(.data[[tcn]]), Similarity)
    cooc_raw <- get_cooccurrence_codes_from_df(tk$cooc, b$tn, tk$icd9_col, tk$target_col, tcn) %>%
      transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(.data[[tcn]]),
                Co_Occurrence_Frequency)

    # everything the matrix could offer, unthresholded. separates "threshold cut
    # it" from "embedding never had it"
    universe <- get_similarity_scores_from_sheets(sheets, 0, tcn) %>%
      transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(.data[[tcn]])) %>%
      distinct()

    # --- stage 2: pool after merge + chapter filter --------------------
    pool <- merge_and_flag(
      sim_raw %>% rename(!!tcn := target),
      cooc_raw %>% rename(!!tcn := target),
      tcn, tk$find_fn, tk$align)
    pool_pairs <- pool %>% transmute(ICD_9_CM, target = .data[[tcn]]) %>% distinct()

    # --- stage 3: what the flag rule actually emits --------------------
    emitted <- select_rows_by_flags(pool, b$flag) %>%
      transmute(ICD_9_CM, target = .data[[tcn]]) %>% distinct()

    in_pool     <- semi_join(tp_df, pool_pairs, by = c("ICD_9_CM", "target"))
    in_emitted  <- semi_join(tp_df, emitted,    by = c("ICD_9_CM", "target"))
    in_universe <- semi_join(tp_df, universe,   by = c("ICD_9_CM", "target"))
    in_sim      <- semi_join(tp_df, sim_raw %>% select(ICD_9_CM, target), by = c("ICD_9_CM", "target"))
    in_cooc     <- semi_join(tp_df, cooc_raw %>% select(ICD_9_CM, target), by = c("ICD_9_CM", "target"))

    # oracle, perfect selection out of the pool
    oracle <- f1_of(tp = nrow(in_pool), fp = 0, fn = nrow(tp_df) - nrow(in_pool))
    # actual, recomputed here as a cross check on the main pipeline
    tp_n <- nrow(in_emitted)
    actual <- f1_of(tp = tp_n, fp = nrow(emitted) - tp_n, fn = nrow(tp_df) - tp_n)

    cat(sprintf("\n--- %s / %s  (thr %.3f, top_n %d, flag %d) ---\n", tr, mdl, b$thr, b$tn, b$flag))
    cat(sprintf("  true pairs                        : %d\n", nrow(tp_df)))
    cat(sprintf("  reachable in similarity universe  : %d (%.1f%%)\n", nrow(in_universe), 100*nrow(in_universe)/nrow(tp_df)))
    cat(sprintf("  survive similarity threshold      : %d (%.1f%%)\n", nrow(in_sim), 100*nrow(in_sim)/nrow(tp_df)))
    cat(sprintf("  present in top-N co-occurrence    : %d (%.1f%%)\n", nrow(in_cooc), 100*nrow(in_cooc)/nrow(tp_df)))
    cat(sprintf("  IN CANDIDATE POOL (recall ceiling): %d (%.1f%%)\n", nrow(in_pool), 100*nrow(in_pool)/nrow(tp_df)))
    cat(sprintf("  actually emitted (true positives) : %d (%.1f%%)\n", tp_n, 100*tp_n/nrow(tp_df)))
    cat(sprintf("  pool size / emitted size          : %d / %d pairs\n", nrow(pool_pairs), nrow(emitted)))
    cat(sprintf("  candidates per ICD-9 (pool/emit)  : %.1f / %.1f\n",
                nrow(pool_pairs)/n_distinct(pool_pairs$ICD_9_CM),
                nrow(emitted)/max(n_distinct(emitted$ICD_9_CM),1)))
    cat(sprintf("  ACTUAL  P %.3f R %.3f F1 %.3f\n", actual$precision, actual$recall, actual$f1))
    cat(sprintf("  ORACLE  P %.3f R %.3f F1 %.3f   <-- ceiling for any reranker on this pool\n",
                oracle$precision, oracle$recall, oracle$f1))

    rows[[length(rows)+1]] <- tibble(
      track = tr, model = mdl,
      n_true = nrow(tp_df),
      pct_in_universe = round(100*nrow(in_universe)/nrow(tp_df), 1),
      pct_pass_threshold = round(100*nrow(in_sim)/nrow(tp_df), 1),
      pct_in_cooc = round(100*nrow(in_cooc)/nrow(tp_df), 1),
      pct_in_pool = round(100*nrow(in_pool)/nrow(tp_df), 1),
      pct_emitted = round(100*tp_n/nrow(tp_df), 1),
      pool_pairs = nrow(pool_pairs), emitted_pairs = nrow(emitted),
      actual_f1 = round(actual$f1, 3), oracle_f1 = round(oracle$f1, 3),
      headroom_f1 = round(oracle$f1 - actual$f1, 3))

    # where did the LOST true pairs die?
    lost <- anti_join(tp_df, pool_pairs, by = c("ICD_9_CM", "target"))
    lost_diag <- lost %>%
      left_join(universe %>% mutate(in_universe = TRUE), by = c("ICD_9_CM", "target")) %>%
      left_join(sim_raw %>% select(ICD_9_CM, target) %>% mutate(in_sim = TRUE), by = c("ICD_9_CM", "target")) %>%
      left_join(cooc_raw %>% select(ICD_9_CM, target) %>% mutate(in_cooc = TRUE), by = c("ICD_9_CM", "target")) %>%
      mutate(
        chapter_ok = !is.na(compute_chapter_distance(
          find_icd9cm_chapter(ICD_9_CM), tk$find_fn(target), tk$align)) &
          compute_chapter_distance(find_icd9cm_chapter(ICD_9_CM), tk$find_fn(target), tk$align) < 1,
        reason = case_when(
          !chapter_ok                              ~ "killed by chapter filter",
          is.na(in_universe)                       ~ "target absent from similarity matrix",
          is.na(in_sim) & is.na(in_cooc)           ~ "below similarity threshold AND outside top-N co-occurrence",
          is.na(in_sim)                            ~ "below similarity threshold (co-occurrence had it)",
          TRUE                                     ~ "other"))
    cat("  lost true pairs by reason:\n")
    lr <- lost_diag %>% count(reason, sort = TRUE)
    for (i in seq_len(nrow(lr))) cat(sprintf("    %-62s %4d (%.1f%%)\n", lr$reason[i], lr$n[i],
                                             100*lr$n[i]/nrow(tp_df)))
    loss_rows[[length(loss_rows)+1]] <- lr %>% mutate(track = tr, model = mdl,
                                                      pct_of_true = round(100*n/nrow(tp_df), 1))
  }
}

summary_df <- bind_rows(rows)
write.csv(summary_df, out_path("error_analysis_summary.csv"), row.names = FALSE)
write.csv(bind_rows(loss_rows), out_path("error_analysis_losses.csv"), row.names = FALSE)

cat("\n\n================ SUMMARY ================\n")
print(as.data.frame(summary_df))
cat("\nWritten to results/error_analysis_summary.csv and results/error_analysis_losses.csv\n")
