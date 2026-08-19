# Task 1 + 2, FULL parameter grid.
#
# scripts/06_extended_comparison.R could only afford ONE (threshold, top_n)
# operating point per track per embedding condition, because merge_and_flag()
# cost ~10-12s per point under the old rowwise() implementation. Picking a
# single arbitrary operating point is a weak basis for the headline claim
# ("does a general-purpose embedding model beat a domain-specific one?"):
# each condition has a different similarity distribution, so the point that
# is optimal for one model need not be optimal for another, and a comparison
# at one shared point can rank models by an artifact of that choice.
#
# find_chapter()/compute_chapter_distance() in pipeline_lib.R are now
# vectorized (bit-identical output, verified by
# scripts/verify_vectorized_equivalence.R), which makes merge_and_flag()
# ~130x faster and the full sweep affordable for every condition. This
# script therefore runs the SAME full grid for all of them:
#
#   thresholds {0.95, 0.99, 0.995, 0.999}
#   x top_n    {3, 5, 10, 15, 20, 25, 30}  (superset of every range used before)
#   x flags    {1, 2, 3, 4}                = 112 points per condition per track
#
# The top_n range is a strict superset of both ranges previously used
# ({3,5,10,15} originally, {5,10,15,20,25,30} after the widening in
# 02_run_comparison.R), so no earlier result can be lost by using it, and the
# widened tail is now searched for BOTH tracks rather than just one.
#
# Usage (from scripts/):
#   Rscript 07_full_grid_comparison.R                 # every condition, sequentially
#   Rscript 07_full_grid_comparison.R mpnet_base ...  # named conditions only
# Each condition writes its own results/full_grid/<tag>.csv, so several can be
# run as parallel background processes without colliding. Assemble the pieces
# with scripts/08_assemble_full_grid.R.

source("pipeline_lib.R")

ORIG_BASE    <- "../data/original"
SAPBERT_BASE <- "../data/sapbert"
GEN_BASE     <- "../data/generated"
OUT_DIR      <- "../results/full_grid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

thresholds <- c(0.95, 0.99, 0.995, 0.999)
top_ns     <- c(3, 5, 10, 15, 20, 25, 30)
flags      <- 1:4

# Each condition names the embedding matrices to score. "family" and
# "stripping" are carried into the output so the two questions this grid is
# meant to answer (which model family wins; does filler-word stripping help)
# can be read off directly.
# NOTE ON TAG NAMING: tags become output FILENAMES, and this pipeline runs on
# Windows, where the filesystem is case-insensitive. Two tags differing only
# by case ("SapBERT_base" vs "sapbert_base") are distinct R list names but the
# SAME file on disk, so parallel background runs silently overwrite each
# other's results -- which happened once and produced a results file whose
# "model" column did not match its filename. Tags must therefore be unique
# case-INSENSITIVELY; the check below enforces it.
CONDITIONS <- list(
  # The two arms the existing baseline was built on. clinicalbert_original
  # matrices came from an external pipeline not present in this repo (its
  # generation method is undocumented -- flagged for the PI); sapbert_original
  # generation IS documented and validated here.
  clinicalbert_original = list(
    model = "ClinicalBERT-original", family = "ClinicalBERT", stripping = "base",
    path_10_9 = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
    path_8_9  = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")),
  sapbert_original = list(
    model = "SapBERT-base", family = "SapBERT", stripping = "base",
    path_10_9 = file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
    path_8_9  = file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx")),
  # Regenerated arms: all produced by the same generation script and pooling
  # (scripts/generate_embeddings.py), so base-vs-stripped and
  # family-vs-family are controlled comparisons. SapBERT-base-regen doubles
  # as a control on the regeneration itself: it should land on top of
  # SapBERT-base if the regeneration reproduces data/sapbert.
  clinicalbert_base     = list(model = "ClinicalBERT-base",     family = "ClinicalBERT", stripping = "base"),
  clinicalbert_stripped = list(model = "ClinicalBERT-stripped", family = "ClinicalBERT", stripping = "stripped"),
  sapbert_base          = list(model = "SapBERT-base-regen",    family = "SapBERT",      stripping = "base"),
  sapbert_stripped      = list(model = "SapBERT-stripped",      family = "SapBERT",      stripping = "stripped"),
  mpnet_base            = list(model = "mpnet-base",            family = "mpnet",        stripping = "base"),
  mpnet_stripped        = list(model = "mpnet-stripped",        family = "mpnet",        stripping = "stripped")
)
dupe <- names(CONDITIONS)[duplicated(tolower(names(CONDITIONS)))]
if (length(dupe)) {
  stop("condition tags must be unique case-insensitively (they become filenames on a ",
       "case-insensitive filesystem); offending: ", paste(dupe, collapse = ", "))
}

for (tag in names(CONDITIONS)) {
  if (is.null(CONDITIONS[[tag]]$path_10_9)) {
    CONDITIONS[[tag]]$path_10_9 <- file.path(GEN_BASE, sprintf("cosine_similarity_matrices_10_9_%s.xlsx", tag))
    CONDITIONS[[tag]]$path_8_9  <- file.path(GEN_BASE, sprintf("cosine_similarity_matrices_8_9_%s.xlsx", tag))
  }
}

TRACKS <- list(
  `10_9` = list(target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
                manual_target_col = "ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"),
  `8_9`  = list(target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
                manual_target_col = "ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx")
)

args <- commandArgs(trailingOnly = TRUE)
selected <- if (length(args) == 0) names(CONDITIONS) else args
unknown <- setdiff(selected, names(CONDITIONS))
if (length(unknown)) {
  stop("unknown condition(s): ", paste(unknown, collapse = ", "),
       "\nknown: ", paste(names(CONDITIONS), collapse = ", "))
}

cat("Loading shared reference data...\n")
VAL_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
ccs_df <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                     sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
for (tr in names(TRACKS)) {
  TRACKS[[tr]]$manual <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$manual_sheet)
  TRACKS[[tr]]$excl   <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$excl_sheet)
  TRACKS[[tr]]$cooc   <- load_cooccurrence_df(file.path(ORIG_BASE, TRACKS[[tr]]$cooc_file))
}

run_full_grid <- function(cond, track_name, sheets) {
  tk  <- TRACKS[[track_name]]
  tcn <- tk$target_col_name
  n_total <- length(thresholds) * length(top_ns) * length(flags)
  cat(sprintf("--- %s / %s (%d points) ---\n", track_name, cond$model, n_total))

  # co-occurrence candidates depend only on top_n, so build each one once and
  # reuse it across all four thresholds
  cooc_by_top_n <- lapply(top_ns, function(tn)
    get_cooccurrence_codes_from_df(tk$cooc, tn, tk$icd9_col, tk$target_col, tcn))
  names(cooc_by_top_n) <- as.character(top_ns)

  rows <- list()
  for (thr in thresholds) {
    similarity_df <- get_similarity_scores_from_sheets(sheets, thr, tcn)
    for (tn in top_ns) {
      merged_df <- merge_and_flag(similarity_df, cooc_by_top_n[[as.character(tn)]], tcn,
                                  tk$find_target_chapter_fn, tk$chapter_alignment)
      for (fc in flags) {
        auto_df <- select_rows_by_flags(merged_df, fc)
        final_valid_df <- validate_mapping(tk$manual, auto_df, ccs_df, tk$excl,
                                           manual_target_col = tk$manual_target_col,
                                           target_col_name = tcn)
        m <- calculate_performance_metrics(final_valid_df)$overall
        rows[[length(rows) + 1]] <- tibble(
          track = track_name, model = cond$model,
          family = cond$family, stripping = cond$stripping,
          similarity_threshold = thr, top_n = tn, flag_combination = fc,
          precision = m$overall_precision, recall = m$overall_recall,
          f1 = m$overall_f1_score, accuracy = m$overall_accuracy,
          n_auto_mappings = nrow(auto_df)
        )
      }
    }
    cat(sprintf("    threshold %.3f done (%d/%d)\n", thr, length(rows), n_total))
  }
  bind_rows(rows)
}

for (tag in selected) {
  cond <- CONDITIONS[[tag]]
  if (!file.exists(cond$path_10_9) || !file.exists(cond$path_8_9)) {
    cat(sprintf("SKIP %s: missing similarity matrices (%s)\n", cond$model, cond$path_10_9))
    next
  }
  t0 <- Sys.time()
  cat(sprintf("\n=== condition %s (%s) ===\n", tag, cond$model))
  out <- bind_rows(
    run_full_grid(cond, "10_9", load_similarity_sheets(cond$path_10_9)),
    run_full_grid(cond, "8_9",  load_similarity_sheets(cond$path_8_9))
  )
  out_path <- file.path(OUT_DIR, sprintf("%s.csv", tag))
  write.csv(out, out_path, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows, %.1f min)\n", out_path, nrow(out),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  best <- out %>% group_by(track) %>% slice_max(order_by = f1, n = 1, with_ties = FALSE) %>% ungroup()
  print(best %>% select(track, model, similarity_threshold, top_n, flag_combination, f1, accuracy))
}

cat("\nDone.\n")
