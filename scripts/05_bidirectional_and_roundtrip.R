source("pipeline_lib.R")

# Reverse mapping (target -> ICD-9-CM) and round-trip consistency.
# Reuses the forward parameters from best_by_model.csv instead of running a
# separate reverse grid search.

ORIG_BASE    <- "../data/original"
SAPBERT_BASE <- "../data/sapbert"
OUT_DIR      <- "../results"

best_by_model <- read.csv(file.path(OUT_DIR, "best_by_model.csv"), stringsAsFactors = FALSE)

ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                      sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)

manual_10_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10")
excl_10_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10_Excld")
cooc_10_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))

manual_8_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validaion_ICD9_ICD8")
excl_8_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD8_Excld")
cooc_8_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))

cat("Loading similarity sheets...\n")
sheets_clin_10_9 <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))
sheets_clin_8_9  <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"))
sheets_sap_10_9  <- load_similarity_sheets(file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"))
sheets_sap_8_9   <- load_similarity_sheets(file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx"))

get_params <- function(track, model_name) {
  b <- best_by_model %>% filter(track == !!track, model == !!model_name)
  if (nrow(b) == 0) stop(sprintf("no best_by_model row for %s / %s", track, model_name))
  list(similarity_threshold = b$similarity_threshold[1], top_n = b$top_n[1], flag_combination = b$flag_combination[1])
}

run_bidir <- function(track, model_name, sheets, cooc_df, manual_df, excl_df,
                       fwd_fn, rev_fn, target_col_name) {
  p <- get_params(track, model_name)
  cat(sprintf("\n=== %s / %s (threshold=%.3f, top_n=%d, flag=%d) ===\n",
              track, model_name, p$similarity_threshold, p$top_n, p$flag_combination))

  fwd <- fwd_fn(sheets, cooc_df, manual_df, ccs_df, excl_df,
                similarity_threshold = p$similarity_threshold, top_n = p$top_n,
                flag_combination = p$flag_combination)
  rev <- rev_fn(sheets, cooc_df, manual_df, ccs_df, excl_df,
                similarity_threshold = p$similarity_threshold, top_n = p$top_n,
                flag_combination = p$flag_combination)

  cat(sprintf("Forward (ICD-9 -> %s):  P%.3f R%.3f F1 %.3f Acc %.3f\n",
              target_col_name, fwd$overall$overall_precision, fwd$overall$overall_recall,
              fwd$overall$overall_f1_score, fwd$overall$overall_accuracy))
  cat(sprintf("Reverse (%s -> ICD-9):  P%.3f R%.3f F1 %.3f Acc %.3f\n",
              target_col_name, rev$overall$overall_precision, rev$overall$overall_recall,
              rev$overall$overall_f1_score, rev$overall$overall_accuracy))

  rt <- compute_roundtrip_consistency(fwd$auto_df, rev$auto_df, target_col_name)
  cat(sprintf("Round-trip consistency: %d/%d = %.3f\n", rt$n_consistent, rt$n_total, rt$rate))

  tibble(
    track = track, model = model_name,
    fwd_precision = fwd$overall$overall_precision, fwd_recall = fwd$overall$overall_recall,
    fwd_f1 = fwd$overall$overall_f1_score, fwd_accuracy = fwd$overall$overall_accuracy,
    rev_precision = rev$overall$overall_precision, rev_recall = rev$overall$overall_recall,
    rev_f1 = rev$overall$overall_f1_score, rev_accuracy = rev$overall$overall_accuracy,
    roundtrip_rate = rt$rate, roundtrip_n_consistent = rt$n_consistent, roundtrip_n_total = rt$n_total
  )
}

results <- bind_rows(
  run_bidir("10_9", "ClinicalBERT", sheets_clin_10_9, cooc_10_9, manual_10_9, excl_10_9,
            run_pipeline_10_9_cached, run_pipeline_10_9_reverse_cached, "ICD_10_CA"),
  run_bidir("10_9", "SapBERT", sheets_sap_10_9, cooc_10_9, manual_10_9, excl_10_9,
            run_pipeline_10_9_cached, run_pipeline_10_9_reverse_cached, "ICD_10_CA"),
  run_bidir("8_9", "ClinicalBERT", sheets_clin_8_9, cooc_8_9, manual_8_9, excl_8_9,
            run_pipeline_8_9_cached, run_pipeline_8_9_reverse_cached, "ICDA_8"),
  run_bidir("8_9", "SapBERT", sheets_sap_8_9, cooc_8_9, manual_8_9, excl_8_9,
            run_pipeline_8_9_cached, run_pipeline_8_9_reverse_cached, "ICDA_8")
)

write.csv(results, file.path(OUT_DIR, "bidirectional_roundtrip.csv"), row.names = FALSE)
cat("\n=== Summary ===\n")
print(results)
cat("\nWrote results/bidirectional_roundtrip.csv\n")
