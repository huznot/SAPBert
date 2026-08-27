# checks the vectorized find_chapter / compute_chapter_distance give exactly the
# same output as the old scalar rowwise version
#
# the refactor is purely for speed, it exists so a full sweep per condition is
# affordable, so it must not move a single metric
#
# compares chapter lookups, the full merge_and_flag output, and all four metrics
# under all four flag rules, on real data for both tracks and both directions

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

OLD_LIB <- tempfile(fileext = ".R")
system2("git", c("show", "HEAD:scripts/pipeline_lib.R"), stdout = OLD_LIB)

source("pipeline_lib.R")          # NEW (vectorized) into globalenv
old <- new.env(parent = globalenv())
sys.source(OLD_LIB, envir = old)  # OLD (scalar/rowwise)

ORIG_BASE <- "data/original"
stopifnot(file.exists(file.path(ORIG_BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx")))

ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                     sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
manual_10_9 <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10")
excl_10_9   <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD10_Excld")
cooc_10_9   <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))
manual_8_9  <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validaion_ICD9_ICD8")
excl_8_9    <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"), sheet = "Validation_ICD9_ICD8_Excld")
cooc_8_9    <- load_cooccurrence_df(file.path(ORIG_BASE, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))

sheets_10_9 <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))
sheets_8_9  <- load_similarity_sheets(file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"))

failures <- 0L
check <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "OK  " else "FAIL", label))
  if (!isTRUE(ok)) failures <<- failures + 1L
  invisible(ok)
}

compare_track <- function(track, sheets, cooc_df, manual_df, excl_df, target_col_name,
                          icd9_col, target_col, new_chapter_fn, old_chapter_fn,
                          alignment, alignment_rev, manual_target_col, thr, tn) {
  cat(sprintf("\n=== %s (threshold %.3f, top_n %d) ===\n", track, thr, tn))

  sim_df  <- get_similarity_scores_from_sheets(sheets, thr, target_col_name)
  cooc    <- get_cooccurrence_codes_from_df(cooc_df, tn, icd9_col, target_col, target_col_name)

  # --- 1. chapter lookup over every distinct code in play ---
  codes9 <- unique(c(sim_df$ICD_9_CM, cooc$ICD_9_CM))
  codesT <- unique(c(sim_df[[target_col_name]], cooc[[target_col_name]]))
  check(sprintf("find_icd9cm_chapter over %d distinct ICD-9-CM codes", length(codes9)),
        identical(as.numeric(find_icd9cm_chapter(codes9)),
                  as.numeric(vapply(codes9, function(x) as.numeric(old$find_icd9cm_chapter(x)), numeric(1),
                                    USE.NAMES = FALSE))))
  check(sprintf("target chapter fn over %d distinct target codes", length(codesT)),
        identical(as.numeric(new_chapter_fn(codesT)),
                  as.numeric(vapply(codesT, function(x) as.numeric(old_chapter_fn(x)), numeric(1),
                                    USE.NAMES = FALSE))))

  # --- 2. full merge_and_flag output, forward direction ---
  t_new <- system.time(m_new <- merge_and_flag(sim_df, cooc, target_col_name, new_chapter_fn, alignment))[["elapsed"]]
  t_old <- system.time(m_old <- old$merge_and_flag(sim_df, cooc, target_col_name, old_chapter_fn, alignment))[["elapsed"]]
  key_cols <- c("ICD_9_CM", target_col_name, "Similarity", "Co_Occurrence_Frequency",
                "highest_similarity_flag", "highest_cooccurrence_flag", "exists_in_both_flag")
  check(sprintf("merge_and_flag output identical (%d rows)", nrow(m_new)),
        isTRUE(all.equal(as.data.frame(m_new[key_cols]), as.data.frame(m_old[key_cols]),
                         check.attributes = FALSE)))
  check("merge_and_flag chapter_distance identical",
        isTRUE(all.equal(as.numeric(m_new$chapter_distance), as.numeric(m_old$chapter_distance))))
  cat(sprintf("       merge_and_flag: old %.1fs -> new %.1fs (%.0fx faster)\n",
              t_old, t_new, t_old / max(t_new, 0.001)))

  # --- 3. metrics for all 4 flag rules ---
  for (fc in 1:4) {
    mn <- calculate_performance_metrics(validate_mapping(
      manual_df, select_rows_by_flags(m_new, fc), ccs_df, excl_df,
      manual_target_col = manual_target_col, target_col_name = target_col_name))$overall
    mo <- calculate_performance_metrics(old$validate_mapping(
      manual_df, old$select_rows_by_flags(m_old, fc), ccs_df, excl_df,
      manual_target_col = manual_target_col, target_col_name = target_col_name))$overall
    check(sprintf("flag %d metrics identical (P=%.3f R=%.3f F1=%.3f Acc=%.3f)",
                  fc, mn$overall_precision, mn$overall_recall, mn$overall_f1_score, mn$overall_accuracy),
          identical(unlist(mn), unlist(mo)))
  }

  # --- 4. reverse direction (task 3 bidirectional path uses the same helpers) ---
  sim_rev  <- get_similarity_scores_from_sheets_reverse(sheets, thr, target_col_name)
  cooc_rev <- get_cooccurrence_codes_from_df_reverse(cooc_df, tn, icd9_col, target_col, target_col_name)
  r_new <- merge_and_flag_reverse(sim_rev, cooc_rev, target_col_name, new_chapter_fn, alignment_rev)
  r_old <- old$merge_and_flag_reverse(sim_rev, cooc_rev, target_col_name, old_chapter_fn, alignment_rev)
  check(sprintf("merge_and_flag_reverse output identical (%d rows)", nrow(r_new)),
        isTRUE(all.equal(as.data.frame(r_new[key_cols]), as.data.frame(r_old[key_cols]),
                         check.attributes = FALSE)))
}

compare_track("10_9", sheets_10_9, cooc_10_9, manual_10_9, excl_10_9, "ICD_10_CA",
              "ICD_9_CM_Code3", "ICD_10_CA_Code3", find_icd10ca_chapter, old$find_icd10ca_chapter,
              chapter_alignment_10, chapter_alignment_10_rev, "ICD-10-CA", thr = 0.995, tn = 15)
compare_track("8_9", sheets_8_9, cooc_8_9, manual_8_9, excl_8_9, "ICDA_8",
              "ICD_9_CM_Code", "ICDA_8_Code", find_icda8_chapter, old$find_icda8_chapter,
              chapter_alignment_8, chapter_alignment_8_rev, "ICDA-8", thr = 0.99, tn = 5)

cat(sprintf("\n%s: %d check(s) failed\n", if (failures == 0) "PASS" else "FAIL", failures))
if (failures > 0) quit(status = 1)
