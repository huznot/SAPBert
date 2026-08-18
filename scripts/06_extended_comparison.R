source("pipeline_lib.R")

# Tasks 1 + 2: filler-word stripping and the general-purpose embedding
# model (all-mpnet-base-v2) arm.
#
# Compares 6 embedding conditions against the original baseline:
#   ClinicalBERT / base (regenerated, for an apples-to-apples check against
#     the original externally-generated ClinicalBERT matrices -- see note
#     below), ClinicalBERT / stripped, SapBERT / base (regenerated, same
#     check against data/sapbert), SapBERT / stripped, mpnet / base,
#     mpnet / stripped.
#
# NOTE on "base" regeneration: the original data/original ClinicalBERT
# matrices were produced by an external pipeline not present in this repo
# (no model name/pooling method documented anywhere -- flagged separately
# for the PI). To isolate the effect of filler-word stripping cleanly, this
# script regenerates a ClinicalBERT "base" arm with the *same* generation
# script/methodology as the "stripped" arm (scripts/generate_embeddings.py,
# emilyalsentzer/Bio_ClinicalBERT, CLS-token pooling), so base vs stripped
# is a controlled comparison. The regenerated "base" arm is expected to
# land close to, but not necessarily exactly match, the original externally
# generated numbers -- any gap there is attributable to the unknown
# generation choices in the original pipeline, not to this task's changes.
# SapBERT / base uses the *existing* data/sapbert matrices directly (that
# generation method IS documented and already validated in this repo), so
# SapBERT / base == the already-reported baseline exactly.
#
# GRID NOTE: merge_and_flag() costs ~20-30s per (threshold, top_n) point on
# this machine (rowwise() over ~2-4k candidate rows), confirmed by timing a
# single reverse-direction call. The full parameter sweep in
# 02_run_comparison.R (~96 points for 10_9, ~64 for 8_9, per model) already
# takes ~15-20 min per model x track; running that same sweep for 7 more
# embedding conditions x 2 tracks would be several hours, which isn't
# practical here. So this script fixes (threshold, top_n) at the single
# point that won for the *unstripped* model on each track in the original
# full grid (10_9: threshold 0.995, top_n 30; 8_9: threshold 0.99, top_n 5
# -- picked from best_by_model.csv as a reasonable shared operating point
# rather than re-deriving one per new condition) and only sweeps the 4 flag
# rules, which is the lever the original grid showed mattered most for
# F1. This is a resource trade-off, not a change in methodology: if more
# compute time is available, widen thresholds_10_9/top_ns_10_9/etc. below
# and re-run for a fuller sweep per condition.

ORIG_BASE <- "../data/original"
GEN_BASE  <- "../data/generated"
OUT_DIR   <- "../results"

thresholds_10_9 <- c(0.995)
top_ns_10_9     <- c(30)
thresholds_8_9  <- c(0.99)
top_ns_8_9      <- c(5)
flags           <- 1:4

cat("Loading shared reference data...\n")
ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                      sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)

manual_10_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10")
excl_10_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10_Excld")
cooc_10_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))

manual_8_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validaion_ICD9_ICD8")
excl_8_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD8_Excld")
cooc_8_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))

run_grid_reduced <- function(track, model_name, sheets, cooc_df, manual_df, valid_excluding_df,
                              target_col_name, icd9_col, target_col,
                              find_target_chapter_fn, chapter_alignment, manual_target_col,
                              thresholds, top_ns) {
  rows <- list()
  cat(sprintf("--- %s / %s ---\n", track, model_name))
  for (thr in thresholds) {
    similarity_df <- get_similarity_scores_from_sheets(sheets, thr, target_col_name)
    for (tn in top_ns) {
      cooccurrence_df <- get_cooccurrence_codes_from_df(cooc_df, tn, icd9_col, target_col, target_col_name)
      merged_df <- merge_and_flag(similarity_df, cooccurrence_df, target_col_name,
                                   find_target_chapter_fn, chapter_alignment)
      for (fc in flags) {
        auto_df <- select_rows_by_flags(merged_df, fc)
        final_valid_df <- validate_mapping(manual_df, auto_df, ccs_df, valid_excluding_df,
                                            manual_target_col = manual_target_col,
                                            target_col_name = target_col_name)
        metrics <- calculate_performance_metrics(final_valid_df)
        rows[[length(rows) + 1]] <- tibble(
          track = track, model = model_name,
          similarity_threshold = thr, top_n = tn, flag_combination = fc,
          precision = metrics$overall$overall_precision, recall = metrics$overall$overall_recall,
          f1 = metrics$overall$overall_f1_score, accuracy = metrics$overall$overall_accuracy
        )
      }
    }
  }
  bind_rows(rows)
}

conditions <- list(
  list(model = "ClinicalBERT-base",     file_tag = "clinicalbert_base"),
  list(model = "ClinicalBERT-stripped", file_tag = "clinicalbert_stripped"),
  list(model = "SapBERT-stripped",      file_tag = "sapbert_stripped"),
  list(model = "mpnet-base",            file_tag = "mpnet_base"),
  list(model = "mpnet-stripped",        file_tag = "mpnet_stripped")
)

all_results <- list()

for (cond in conditions) {
  path_10_9 <- file.path(GEN_BASE, sprintf("cosine_similarity_matrices_10_9_%s.xlsx", cond$file_tag))
  path_8_9  <- file.path(GEN_BASE, sprintf("cosine_similarity_matrices_8_9_%s.xlsx", cond$file_tag))
  if (!file.exists(path_10_9) || !file.exists(path_8_9)) {
    cat(sprintf("SKIP %s: generated matrices not found yet (%s)\n", cond$model, path_10_9))
    next
  }
  sheets_10_9 <- load_similarity_sheets(path_10_9)
  sheets_8_9  <- load_similarity_sheets(path_8_9)

  all_results[[length(all_results) + 1]] <- run_grid_reduced(
    "10_9", cond$model, sheets_10_9, cooc_10_9, manual_10_9, excl_10_9,
    target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
    find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
    manual_target_col = "ICD-10-CA", thresholds = thresholds_10_9, top_ns = top_ns_10_9)

  all_results[[length(all_results) + 1]] <- run_grid_reduced(
    "8_9", cond$model, sheets_8_9, cooc_8_9, manual_8_9, excl_8_9,
    target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
    find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
    manual_target_col = "ICDA-8", thresholds = thresholds_8_9, top_ns = top_ns_8_9)
}

# SapBERT / base is the existing, already-validated data/sapbert matrices --
# no need to regenerate, use exactly what 02_run_comparison.R already used.
sheets_sap_10_9 <- load_similarity_sheets("../data/sapbert/cosine_similarity_matrices_10_9_SapBERT.xlsx")
sheets_sap_8_9  <- load_similarity_sheets("../data/sapbert/cosine_similarity_matrices_8_9_SapBERT.xlsx")
all_results[[length(all_results) + 1]] <- run_grid_reduced(
  "10_9", "SapBERT-base", sheets_sap_10_9, cooc_10_9, manual_10_9, excl_10_9,
  target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
  find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
  manual_target_col = "ICD-10-CA", thresholds = thresholds_10_9, top_ns = top_ns_10_9)
all_results[[length(all_results) + 1]] <- run_grid_reduced(
  "8_9", "SapBERT-base", sheets_sap_8_9, cooc_8_9, manual_8_9, excl_8_9,
  target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
  find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
  manual_target_col = "ICDA-8", thresholds = thresholds_8_9, top_ns = top_ns_8_9)

extended_grid <- bind_rows(all_results)
write.csv(extended_grid, file.path(OUT_DIR, "extended_grid_results.csv"), row.names = FALSE)

best_extended <- extended_grid %>%
  group_by(track, model) %>%
  slice_max(order_by = f1, n = 1, with_ties = FALSE) %>%
  ungroup()
write.csv(best_extended, file.path(OUT_DIR, "extended_best_by_model.csv"), row.names = FALSE)

cat("\n=== Best combination per track x model (extended, reduced grid) ===\n")
print(best_extended)

old_best <- read.csv(file.path(OUT_DIR, "best_by_model.csv"), stringsAsFactors = FALSE)
cat("\n=== Old (original grid) vs new (this run) best F1/Accuracy ===\n")
comparison <- bind_rows(
  old_best %>% mutate(source = "original_full_grid"),
  best_extended %>% select(track, model, similarity_threshold, top_n, flag_combination, precision, recall, f1, accuracy) %>%
    mutate(source = "extended_reduced_grid")
) %>% arrange(track, model, source)
write.csv(comparison, file.path(OUT_DIR, "old_vs_new_comparison.csv"), row.names = FALSE)
print(comparison)
