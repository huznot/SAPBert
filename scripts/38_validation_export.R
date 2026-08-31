# exports the validation data as csv, with labels and chapters attached, and
# builds a review sheet for the icd-9 codes recorded as having no match.
#
# two things prompted this. the validation workbook is the only place the
# "no match" decision is recorded and it was not readable outside excel, and
# code 338 turned out to have no reachable correct target: g89, chronic pain,
# is not in the icd-10-ca code list at all, and r52, which carries 338's own
# label, is chapter 18 against 338's aligned chapters 6, 7 and 8, so the
# chapter filter removes it. that is a property of the reference data, not of
# the algorithm, and the other 8 need the same check.
#
# the review sheet carries what the pipeline emits for each code next to the
# chapters, so the manual check has the candidate list in front of it. the
# last two columns are left blank to be filled in by hand.
#
# usage: Rscript 38_validation_export.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

ORIG <- "data/original"; OUT <- "results"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")

ccs <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level")
icd9_lab <- ccs %>% transmute(icd9 = as.character(ICD_9_CM), icd9_label = ICD_9_CM_LABEL,
                              ccs_id = CCS_ID, ccs_category = CCS_CATEGORY_DESCRIPTION)

TRACKS <- list(
  `10_9` = list(label = "ICD-9-CM to ICD-10-CA", tcn = "ICD_10_CA", manual_col = "ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                labels_sheet = "ICD-10-CA-3Level", chap_fn = find_icd10ca_chapter,
                align = chapter_alignment_10,
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                sim_file = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
                run = run_pipeline_10_9_cached, thr = 0.995, top_n = 25, algo = 4),
  `8_9`  = list(label = "ICD-9-CM to ICDA-8", tcn = "ICDA_8", manual_col = "ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                labels_sheet = "ICDA-8-3Level", chap_fn = find_icda8_chapter,
                align = chapter_alignment_8,
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                sim_file = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"),
                run = run_pipeline_8_9_cached, thr = 0.999, top_n = 5, algo = 3)
)

for (tr in names(TRACKS)) {
  tk  <- TRACKS[[tr]]
  man <- read_excel(VAL, sheet = tk$manual_sheet)
  exc <- read_excel(VAL, sheet = tk$excl_sheet)
  tgt <- read_excel(LAB, sheet = tk$labels_sheet) %>% setNames(c("target", "target_label")) %>%
    mutate(target = as.character(target), target_chapter = tk$chap_fn(target))

  unmatched <- unique(as.character(exc$`ICD-9-CM`))
  cat(sprintf("\n=== %s ===\n%d target codes in the code list, %d icd-9 with a match, %d without\n",
              tk$label, nrow(tgt), n_distinct(man$`ICD-9-CM`), length(unmatched)))

  # ---- the reference crosswalk itself, one row per manually mapped pair ----
  pairs <- man %>%
    transmute(icd9 = as.character(`ICD-9-CM`), target = as.character(.data[[tk$manual_col]])) %>%
    left_join(icd9_lab, by = "icd9") %>%
    left_join(tgt, by = "target") %>%
    mutate(icd9_chapter = find_icd9cm_chapter(icd9),
           aligned_target_chapters = vapply(icd9_chapter,
             function(ch) paste(tk$align[[as.character(ch)]], collapse = "/"), character(1)),
           chapter_distance = compute_chapter_distance(icd9_chapter, target_chapter, tk$align),
           survives_chapter_filter = chapter_distance < 1) %>%
    select(icd9, icd9_label, icd9_chapter, ccs_id, ccs_category,
           target, target_label, target_chapter, aligned_target_chapters,
           chapter_distance, survives_chapter_filter)
  write.csv(pairs, file.path(OUT, sprintf("validation_pairs_%s.csv", tr)), row.names = FALSE)

  # a reference pair the chapter filter would reject is a mapping the pipeline
  # cannot make however good its scores are, so it is worth counting
  cat(sprintf("%d reference pairs, %d of them fail the chapter filter\n",
              nrow(pairs), sum(!pairs$survives_chapter_filter, na.rm = TRUE)))

  # ---- what the pipeline emits for the codes with no reference match ----
  cooc <- load_cooccurrence_df(file.path(ORIG, tk$cooc_file))
  res  <- tk$run(load_similarity_sheets(tk$sim_file), cooc, man, ccs %>% select(ICD_9_CM, CCS_ID), exc,
                 similarity_threshold = tk$thr, top_n = tk$top_n, flag_combination = tk$algo)

  emitted <- res$auto_df %>% filter(ICD_9_CM %in% unmatched) %>%
    transmute(icd9 = as.character(ICD_9_CM), target = as.character(.data[[tk$tcn]]),
              similarity = round(Similarity, 3), cooccurrence = Co_Occurrence_Frequency) %>%
    left_join(tgt, by = "target") %>%
    arrange(icd9, desc(coalesce(similarity, 0)))

  # one row per code even when nothing was emitted, so the 52 icda-8 codes with
  # no output still appear on the sheet to be reviewed
  review <- tibble(icd9 = unmatched) %>%
    left_join(icd9_lab, by = "icd9") %>%
    mutate(icd9_chapter = find_icd9cm_chapter(icd9),
           aligned_target_chapters = vapply(icd9_chapter,
             function(ch) paste(tk$align[[as.character(ch)]], collapse = "/"), character(1))) %>%
    left_join(emitted, by = "icd9", relationship = "one-to-many") %>%
    mutate(pipeline_output = ifelse(is.na(target), "none", target),
           correct_target_should_be = "", present_in_code_list = "", notes = "") %>%
    select(icd9, icd9_label, icd9_chapter, aligned_target_chapters, ccs_category,
           pipeline_output, target_label, target_chapter, similarity, cooccurrence,
           correct_target_should_be, present_in_code_list, notes) %>%
    arrange(icd9)
  write.csv(review, file.path(OUT, sprintf("unmatched_review_%s.csv", tr)), row.names = FALSE)
  cat(sprintf("%d codes with no reference match produced %d mappings between them\n",
              length(unmatched), sum(review$pipeline_output != "none")))
}

cat("\nwrote results/validation_pairs_10_9.csv and validation_pairs_8_9.csv\n")
cat("wrote results/unmatched_review_10_9.csv and unmatched_review_8_9.csv\n")
