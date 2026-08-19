source("pipeline_lib.R")

ORIG_BASE    <- "../data/original"
SAPBERT_BASE <- "../data/sapbert"
OUT_DIR      <- "../results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

thresholds  <- c(0.95, 0.99, 0.995, 0.999)
# Both tracks sweep top-N {3,5,10,15,20,25,30}, a strict SUPERSET of the two
# ranges this script has used historically ({3,5,10,15} originally, and
# {5,10,15,20,25,30} in the first attempt at a fix). Using the superset for
# both tracks means no earlier result can be lost, and the widened tail is
# actually searched on both tracks rather than just one.
#
# WHY THIS MATTERS: the original range {3,5,10,15} was too narrow and cut off
# the real optimum on the ICD-9-CM -> ICD-10-CA track. Under it, the best
# reachable numbers were F1 0.423 (ClinicalBERT) and 0.523 (SapBERT), which is
# why the reference numbers F1 0.427/Acc 0.271 and F1 0.530/Acc 0.360 looked
# unreproducible and were flagged as a discrepancy. They are NOT wrong: with
# the widened range they reproduce exactly, at top_n=30 (ClinicalBERT,
# threshold 0.995) and top_n=30 (SapBERT, threshold 0.95). See report.Rmd.
#
# Sweeping the same range on both tracks is affordable because
# merge_and_flag() is now ~130x faster (vectorized chapter lookup in
# pipeline_lib.R, verified bit-identical by verify_vectorized_equivalence.R).
top_ns_8_9  <- c(3, 5, 10, 15, 20, 25, 30)
top_ns_10_9 <- c(3, 5, 10, 15, 20, 25, 30)
flags       <- 1:4

cat("Loading shared reference data...\n")
ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                      sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)

manual_10_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10")
excl_10_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10_Excld")
cooc_10_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))

manual_8_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validaion_ICD9_ICD8")
excl_8_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD8_Excld")
cooc_8_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))

cat("Loading ClinicalBERT similarity sheets...\n")
sheets_clin_10_9 <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))
sheets_clin_8_9  <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"))

cat("Loading SapBERT similarity sheets...\n")
sheets_sap_10_9 <- load_similarity_sheets(file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"))
sheets_sap_8_9  <- load_similarity_sheets(file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx"))

run_grid <- function(track, model_name, sheets, cooc_df, manual_df, valid_excluding_df,
                      target_col_name, icd9_col, target_col,
                      find_target_chapter_fn, chapter_alignment, manual_target_col, top_ns) {
  n_combinations <- length(thresholds) * length(top_ns) * length(flags)
  cat(sprintf("\n--- Running grid: track=%s model=%s (%d combinations) ---\n",
              track, model_name, n_combinations))
  rows <- list()
  i <- 0L

  for (thr in thresholds) {
    similarity_df <- get_similarity_scores_from_sheets(sheets, thr, target_col_name)

    for (tn in top_ns) {
      cooccurrence_df <- get_cooccurrence_codes_from_df(cooc_df, tn, icd9_col, target_col, target_col_name)
      merged_df <- merge_and_flag(similarity_df, cooccurrence_df, target_col_name,
                                   find_target_chapter_fn, chapter_alignment)

      for (fc in flags) {
        i <- i + 1L
        auto_df <- select_rows_by_flags(merged_df, fc)
        final_valid_df <- validate_mapping(manual_df, auto_df, ccs_df, valid_excluding_df,
                                            manual_target_col = manual_target_col,
                                            target_col_name = target_col_name)
        metrics <- calculate_performance_metrics(final_valid_df)

        rows[[length(rows) + 1]] <- tibble(
          track = track, model = model_name,
          similarity_threshold = thr, top_n = tn, flag_combination = fc,
          precision = metrics$overall$overall_precision,
          recall    = metrics$overall$overall_recall,
          f1        = metrics$overall$overall_f1_score,
          accuracy  = metrics$overall$overall_accuracy
        )
      }
    }
    cat(sprintf("  threshold %.3f done (%d/%d combinations)\n", thr, i, n_combinations))
  }
  bind_rows(rows)
}

res_clin_10_9 <- run_grid("10_9", "ClinicalBERT", sheets_clin_10_9, cooc_10_9, manual_10_9, excl_10_9,
                           target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                           find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
                           manual_target_col = "ICD-10-CA", top_ns = top_ns_10_9)
res_sap_10_9  <- run_grid("10_9", "SapBERT", sheets_sap_10_9, cooc_10_9, manual_10_9, excl_10_9,
                           target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                           find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
                           manual_target_col = "ICD-10-CA", top_ns = top_ns_10_9)
res_clin_8_9  <- run_grid("8_9", "ClinicalBERT", sheets_clin_8_9, cooc_8_9, manual_8_9, excl_8_9,
                           target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                           find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
                           manual_target_col = "ICDA-8", top_ns = top_ns_8_9)
res_sap_8_9   <- run_grid("8_9", "SapBERT", sheets_sap_8_9, cooc_8_9, manual_8_9, excl_8_9,
                           target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                           find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
                           manual_target_col = "ICDA-8", top_ns = top_ns_8_9)

grid_10_9 <- bind_rows(res_clin_10_9, res_sap_10_9)
grid_8_9  <- bind_rows(res_clin_8_9, res_sap_8_9)

write.csv(grid_10_9, file.path(OUT_DIR, "grid_results_10_9.csv"), row.names = FALSE)
write.csv(grid_8_9,  file.path(OUT_DIR, "grid_results_8_9.csv"), row.names = FALSE)

all_grid <- bind_rows(grid_10_9, grid_8_9)

best_by_model <- all_grid %>%
  group_by(track, model) %>%
  slice_max(order_by = f1, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(best_by_model, file.path(OUT_DIR, "best_by_model.csv"), row.names = FALSE)

cat("\n=== Best combination per track x model (by F1) ===\n")
print(best_by_model)

get_best_grouped <- function(track, model_name, sheets, cooc_df, manual_df, valid_excluding_df, pipeline_fn) {
  b <- best_by_model %>% filter(track == !!track, model == !!model_name)
  res <- pipeline_fn(sheets, cooc_df, manual_df, ccs_df, valid_excluding_df,
                      similarity_threshold = b$similarity_threshold,
                      top_n = b$top_n,
                      flag_combination = b$flag_combination)
  res$grouped %>% mutate(track = track, model = model_name)
}

ccs_10_9 <- bind_rows(
  get_best_grouped("10_9", "ClinicalBERT", sheets_clin_10_9, cooc_10_9, manual_10_9, excl_10_9, run_pipeline_10_9_cached),
  get_best_grouped("10_9", "SapBERT",      sheets_sap_10_9,  cooc_10_9, manual_10_9, excl_10_9, run_pipeline_10_9_cached)
)
ccs_8_9 <- bind_rows(
  get_best_grouped("8_9", "ClinicalBERT", sheets_clin_8_9, cooc_8_9, manual_8_9, excl_8_9, run_pipeline_8_9_cached),
  get_best_grouped("8_9", "SapBERT",      sheets_sap_8_9,  cooc_8_9, manual_8_9, excl_8_9, run_pipeline_8_9_cached)
)

write.csv(ccs_10_9, file.path(OUT_DIR, "ccs_breakdown_10_9.csv"), row.names = FALSE)
write.csv(ccs_8_9,  file.path(OUT_DIR, "ccs_breakdown_8_9.csv"), row.names = FALSE)

cat("\nDone. Results written to", OUT_DIR, "\n")
