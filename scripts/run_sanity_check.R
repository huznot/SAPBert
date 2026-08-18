source("pipeline_lib.R")

BASE <- "../data/original"

cat("Loading 10_9 ClinicalBERT similarity sheets...\n")
sheets_10_9 <- load_similarity_sheets(file.path(BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))
cooc_10_9 <- load_cooccurrence_df(file.path(BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))
manual_10_9 <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10")
ccs_df <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"), sheet = "CCS ICD-9-CM-3Level") %>%
  select(ICD_9_CM, CCS_ID)
excl_10_9 <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10_Excld")

result_10_9 <- run_pipeline_10_9_cached(sheets_10_9, cooc_10_9, manual_10_9, ccs_df, excl_10_9,
                                         similarity_threshold = 0.995, top_n = 10, flag_combination = 3)
cat("\n=== ICD-9->ICD-10-CA sanity check ===\n")
cat("Expect: Precision 0.686, Recall 0.277, F1 0.395, Accuracy 0.246\n")
cat("Got:    Precision", result_10_9$overall$overall_precision,
    ", Recall", result_10_9$overall$overall_recall,
    ", F1", result_10_9$overall$overall_f1_score,
    ", Accuracy", result_10_9$overall$overall_accuracy, "\n")

cat("\nLoading 8_9 ClinicalBERT similarity sheets...\n")
sheets_8_9 <- load_similarity_sheets(file.path(BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"))
cooc_8_9 <- load_cooccurrence_df(file.path(BASE, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))
manual_8_9 <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validaion_ICD9_ICD8")
excl_8_9 <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD8_Excld")

result_8_9 <- run_pipeline_8_9_cached(sheets_8_9, cooc_8_9, manual_8_9, ccs_df, excl_8_9,
                                       similarity_threshold = 0.995, top_n = 5, flag_combination = 3)
cat("\n=== ICD-9->ICDA-8 sanity check ===\n")
cat("Expect: Precision 0.753, Recall 0.680, F1 ~0.714-0.715, Accuracy 0.556\n")
cat("Got:    Precision", result_8_9$overall$overall_precision,
    ", Recall", result_8_9$overall$overall_recall,
    ", F1", result_8_9$overall$overall_f1_score,
    ", Accuracy", result_8_9$overall$overall_accuracy, "\n")
