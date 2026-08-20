# tries to raise the recall ceiling
#
# 09_ showed the embedding isnt the problem, every true target is in the
# similarity matrix but only ~63% (icd-10-ca) / ~85% (icda-8) make it into the
# pool, capping oracle f1 at 0.770 / 0.920. so vary candidate generation itself:
#   A. similarity cut, relative (current) vs topk
#   B. how many co-occurrence candidates
#   C. chapter filter on or off
# the chapter filter gets its own output since it drops 48 true pairs on icd-10-ca
source("pipeline_lib.R")

ORIG_BASE    <- "../data/original"
SAPBERT_BASE <- "../data/sapbert"
GEN_BASE     <- "../data/generated"
OUT_DIR      <- "../results"

ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                     sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
VAL_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

TRACKS <- list(
  `10_9` = list(tcn = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                find_fn = find_icd10ca_chapter, align = chapter_alignment_10,
                manual_target_col = "ICD-10-CA", manual_sheet = "Validation_ICD9_ICD10",
                excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                sim = list(SapBERT = file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                           mpnet   = file.path(GEN_BASE, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"))),
  `8_9`  = list(tcn = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                find_fn = find_icda8_chapter, align = chapter_alignment_8,
                manual_target_col = "ICDA-8", manual_sheet = "Validaion_ICD9_ICD8",
                excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                sim = list(SapBERT = file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx"),
                           mpnet   = file.path(GEN_BASE, "cosine_similarity_matrices_8_9_mpnet_base.xlsx")))
)
for (tr in names(TRACKS)) {
  TRACKS[[tr]]$manual <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$manual_sheet)
  TRACKS[[tr]]$excl   <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$excl_sheet)
  TRACKS[[tr]]$cooc   <- load_cooccurrence_df(file.path(ORIG_BASE, TRACKS[[tr]]$cooc_file))
}

true_pairs <- function(tk) {
  excluded <- as.character(tk$excl$`ICD-9-CM`)
  tk$manual %>%
    transmute(ICD_9_CM = as.character(`ICD-9-CM`),
              target   = as.character(.data[[tk$manual_target_col]])) %>%
    filter(!(ICD_9_CM %in% excluded)) %>% distinct()
}

# full long format similarity, computed once per model x track and reused for
# every design. this is the expensive part
long_similarity <- function(path, tcn) {
  sheets <- load_similarity_sheets(path)
  purrr::map_dfr(sheets, function(df) {
    id_col <- names(df)[1]
    tgt <- as.character(df[[id_col]])
    purrr::map_dfr(setdiff(names(df), id_col), function(col) {
      tibble(ICD_9_CM = as.character(col), target = tgt, Similarity = df[[col]])
    })
  }) %>% filter(!is.na(Similarity))
}

cat("Loading similarity matrices (once)...\n")
LONG <- list()
for (tr in names(TRACKS)) {
  LONG[[tr]] <- list()
  for (mdl in names(TRACKS[[tr]]$sim)) {
    LONG[[tr]][[mdl]] <- long_similarity(TRACKS[[tr]]$sim[[mdl]], TRACKS[[tr]]$tcn)
    cat(sprintf("  %s / %s: %d cells\n", tr, mdl, nrow(LONG[[tr]][[mdl]])))
  }
}

# chapter distance lookup for an arbitrary pair set
chapter_ok <- function(icd9, target, tk) {
  d <- compute_chapter_distance(find_icd9cm_chapter(icd9), tk$find_fn(target), tk$align)
  !is.na(d) & d < 1
}

sim_candidates <- function(long_df, mode, param) {
  if (mode == "relative") {
    long_df %>% group_by(ICD_9_CM) %>%
      filter(Similarity >= param * max(Similarity, na.rm = TRUE)) %>% ungroup()
  } else {
    long_df %>% group_by(ICD_9_CM) %>%
      slice_max(order_by = Similarity, n = param, with_ties = FALSE) %>% ungroup()
  }
}

results <- list()
for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  tp_df <- true_pairs(tk)
  n_true <- nrow(tp_df)
  cat(sprintf("\n######## track %s (%d true pairs) ########\n", tr, n_true))

  cooc_cache <- list()
  for (tn in c(3, 5, 10, 20, 30, 50)) {
    cooc_cache[[as.character(tn)]] <- get_cooccurrence_codes_from_df(
      tk$cooc, tn, tk$icd9_col, tk$target_col, tk$tcn) %>%
      transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(.data[[tk$tcn]]))
  }

  for (mdl in names(tk$sim)) {
    long_df <- LONG[[tr]][[mdl]]
    designs <- c(
      lapply(c(0.999, 0.995, 0.99, 0.95, 0.90, 0.85, 0.80), function(p) list(mode = "relative", param = p)),
      lapply(c(1, 3, 5, 10, 20, 30, 50), function(p) list(mode = "topk", param = p))
    )
    for (d in designs) {
      simc <- sim_candidates(long_df, d$mode, d$param) %>% select(ICD_9_CM, target) %>% distinct()
      for (tn in c(3, 5, 10, 20, 30, 50)) {
        coocc <- cooc_cache[[as.character(tn)]]
        for (chap in c(TRUE, FALSE)) {
          pool <- bind_rows(simc, coocc) %>% distinct()
          if (chap) pool <- pool %>% filter(chapter_ok(ICD_9_CM, target, tk))
          hit <- nrow(semi_join(tp_df, pool, by = c("ICD_9_CM", "target")))
          results[[length(results)+1]] <- tibble(
            track = tr, model = mdl,
            sim_mode = d$mode, sim_param = d$param, top_n = tn, chapter_filter = chap,
            pool_pairs = nrow(pool),
            pool_per_code = round(nrow(pool) / n_distinct(pool$ICD_9_CM), 2),
            recall_ceiling = round(hit / n_true, 4),
            # oracle f1 assumes a perfect selector, precision 1 and recall is
            # the ceiling
            oracle_f1 = round(2 * (hit/n_true) / (1 + hit/n_true), 4))
        }
      }
    }
    cat(sprintf("  %s done\n", mdl))
  }
}

study <- bind_rows(results)
write.csv(study, file.path(OUT_DIR, "candidate_generation_study.csv"), row.names = FALSE)

cat("\n=== Best recall ceiling per track x model x chapter_filter ===\n")
print(as.data.frame(study %>% group_by(track, model, chapter_filter) %>%
  slice_max(recall_ceiling, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(track, model, chapter_filter, sim_mode, sim_param, top_n, recall_ceiling, pool_per_code)))

cat("\n=== Cost of the chapter filter (matched designs) ===\n")
chap_cost <- study %>%
  select(track, model, sim_mode, sim_param, top_n, chapter_filter, recall_ceiling, pool_pairs) %>%
  tidyr::pivot_wider(names_from = chapter_filter,
                     values_from = c(recall_ceiling, pool_pairs),
                     names_prefix = "chap_") %>%
  mutate(recall_lost = round(recall_ceiling_chap_FALSE - recall_ceiling_chap_TRUE, 4),
         pool_saved_pct = round(100 * (1 - pool_pairs_chap_TRUE / pool_pairs_chap_FALSE), 1)) %>%
  group_by(track, model) %>%
  summarise(mean_recall_lost = round(mean(recall_lost), 4),
            mean_pool_saved_pct = round(mean(pool_saved_pct), 1), .groups = "drop")
print(as.data.frame(chap_cost))

# which true pairs the chapter filter destroys
cas <- list()
for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  tp_df <- true_pairs(tk)
  bad <- tp_df %>% filter(!chapter_ok(ICD_9_CM, target, tk)) %>%
    mutate(track = tr,
           chapter_icd9 = find_icd9cm_chapter(ICD_9_CM),
           chapter_target = tk$find_fn(target),
           allowed_target_chapters = vapply(chapter_icd9, function(c9) {
             a <- tk$align[[as.character(c9)]]
             if (is.null(a)) NA_character_ else paste(a, collapse = ",")
           }, character(1)))
  cas[[length(cas)+1]] <- bad
}
casualties <- bind_rows(cas)
write.csv(casualties, file.path(OUT_DIR, "chapter_filter_casualties.csv"), row.names = FALSE)
cat(sprintf("\n=== True pairs destroyed by the chapter filter: %d ===\n", nrow(casualties)))
print(as.data.frame(casualties %>% count(track, chapter_icd9, chapter_target,
                                         allowed_target_chapters, sort = TRUE) %>% head(20)))
cat("\nWritten to results/candidate_generation_study.csv, results/chapter_filter_casualties.csv\n")
