# how the pipeline handles icd-9 codes that have no correct answer
#
# 9 icd-9-cm codes have no icd-10-ca reference match and 52 have no icda-8
# match. they used to be dropped before scoring, so the pipeline was never
# charged for mapping a code with nothing to map to. they are now scored.
#
# f1 cannot show whether that is handled well. a code with no correct answer
# has no true positive to earn, so staying silent on it scores nothing and
# emitting on it only costs precision. the number that does show it is the
# share of those codes the pipeline correctly returns nothing for, which is
# what this script reports alongside the f1 cost
#
# usage: Rscript 35_unmatched_codes.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

ORIG <- "data/original"; SAP <- "data/sapbert"; GEN <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
OUT  <- "results"

ccs <- read_excel(file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                  sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)

TRACKS <- list(
  `10_9` = list(tcn = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                run = run_pipeline_10_9_cached),
  `8_9`  = list(tcn = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                run = run_pipeline_8_9_cached)
)

SIM <- list(
  ClinicalBERT = list(`10_9` = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
                      `8_9`  = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")),
  SapBERT      = list(`10_9` = file.path(SAP, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                      `8_9`  = file.path(SAP, "cosine_similarity_matrices_8_9_SapBERT.xlsx")),
  mpnet        = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_mpnet_base.xlsx"))
)

# the operating point each model wins at, re-selected under the new scoring
best <- read.csv(file.path(OUT, "full_grid_best.csv"), stringsAsFactors = FALSE)
MODEL_KEY <- c(ClinicalBERT = "ClinicalBERT-original", SapBERT = "SapBERT-base", mpnet = "mpnet-base")

rows <- list()
for (tr in names(TRACKS)) {
  tk   <- TRACKS[[tr]]
  man  <- read_excel(VAL, sheet = tk$manual_sheet)
  exc  <- read_excel(VAL, sheet = tk$excl_sheet)
  cooc <- load_cooccurrence_df(file.path(ORIG, tk$cooc_file))
  unmatched <- unique(as.character(exc$`ICD-9-CM`))

  for (mdl in names(SIM)) {
    b <- best[best$track == tr & best$model == MODEL_KEY[mdl], ]
    if (!nrow(b)) { cat(sprintf("  [no grid row for %s / %s, skipped]\n", tr, mdl)); next }

    res <- tk$run(load_similarity_sheets(SIM[[mdl]][[tr]]), cooc, man, ccs, exc,
                  similarity_threshold = b$similarity_threshold[1],
                  top_n = b$top_n[1], flag_combination = b$flag_combination[1])

    emitted_for <- unique(as.character(res$auto_df$ICD_9_CM))
    silent      <- setdiff(unmatched, emitted_for)
    spurious    <- res$auto_df %>%
      filter(ICD_9_CM %in% unmatched) %>%
      distinct(ICD_9_CM, .data[[tk$tcn]])

    rows[[length(rows) + 1]] <- tibble(
      track = tr, model = mdl,
      similarity_threshold = b$similarity_threshold[1],
      top_n = b$top_n[1], flag_combination = b$flag_combination[1],
      n_unmatched = length(unmatched),
      n_correctly_silent = length(silent),
      correct_rejection_rate = round(length(silent) / length(unmatched), 3),
      n_spurious_pairs = nrow(spurious),
      spurious_per_code = round(nrow(spurious) / max(length(unmatched) - length(silent), 1), 2)
    )
  }
}

out <- bind_rows(rows)
write.csv(out, file.path(OUT, "unmatched_code_handling.csv"), row.names = FALSE)

cat("\n=== ICD-9 codes with no reference match ===\n")
cat("correct rejection rate = share of those codes the pipeline emits nothing for\n\n")
print(as.data.frame(out), row.names = FALSE)
cat("\nwrote results/unmatched_code_handling.csv\n")
