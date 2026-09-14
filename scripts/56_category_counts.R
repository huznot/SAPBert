source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# per ccs category counts at the plain cutoffs, for her to pick where to drill.
# 31_category_breakdown.R already broke the categories down, but at the old
# relative threshold settings out of full_grid_best.csv. this is the same shape
# at the cutoff each model actually does best at, and with the raw counts kept
# rather than the per category f1 on its own.
#
# the caution that has to travel with this table: most ccs categories are tiny,
# so a per category f1 is often computed off one or two pairs and swings between
# 0 and 1 for reasons that have nothing to do with the model. the code counts are
# in every sheet so that is always visible, and the single code categories are
# flagged.
#
# icd-9 to icd-10-ca only. the plain cutoffs were only worked out on that track.

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("4_counts_by_ccs_category.xlsx")

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
            mpnet        = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx")

SCENARIO <- c("1" = "top cosine or top co-occurrence",
              "2" = "top cosine or in both lists",
              "3" = "top co-occurrence or in both lists",
              "4" = "any of the three")

ccs_full <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>%
  mutate(ICD_9_CM = as.character(ICD_9_CM))
ccs_df <- ccs_full %>% select(ICD_9_CM, CCS_ID)
ccs_index <- ccs_full %>%
  group_by(CCS_ID) %>%
  summarise(description = first(CCS_CATEGORY_DESCRIPTION), n_codes = n(), .groups = "drop")

man  <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
excl <- read_excel(VAL, sheet = "Validation_ICD9_ICD10_Excld")
cooc <- load_cooccurrence_df(file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))

res <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

cat(sprintf("%d ccs categories over %d icd-9 codes, %d of them hold a single code\n\n",
            nrow(ccs_index), nrow(ccs_full), sum(ccs_index$n_codes == 1)))

rows <- list(); settings <- list()
for (m in names(MODELS)) {
  d <- res[res$model == m & res$mode == "absolute", ]
  b <- d[which.max(d$f1), ]
  cat(sprintf("%-13s cutoff %.2f, top %d, rule %d\n", m, b$threshold, b$top_n, b$flag))

  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[m]]))
  sim <- long_at_cutoff <- NULL
  sheets_long <- purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble(ICD_9_CM = as.character(c1), ICD_10_CA = as.character(df[[id]]), Similarity = df[[c1]]))
  })
  sim <- sheets_long %>% filter(!is.na(Similarity), Similarity >= b$threshold)
  co  <- get_cooccurrence_codes_from_df(cooc, b$top_n, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
  mg  <- merge_and_flag(sim, co, "ICD_10_CA", find_icd10ca_chapter, chapter_alignment_10)
  au  <- select_rows_by_flags(mg, b$flag)
  fin <- validate_mapping(man, au, ccs_df, excl, "ICD-10-CA", "ICD_10_CA")

  # the overall row has to match what 53_ already reported or something is off
  mt <- calculate_performance_metrics(fin)
  stopifnot(sum(fin$`True Positive`, na.rm = TRUE) == b$tp,
            sum(fin$`False Positive`, na.rm = TRUE) == b$fp,
            sum(fin$`False Negative`, na.rm = TRUE) == b$fn)

  per <- fin %>%
    group_by(CCS_ID) %>%
    summarise(`ICD-9 codes seen` = n_distinct(`ICD-9-CM`),
              TP = sum(`True Positive`, na.rm = TRUE),
              FP = sum(`False Positive`, na.rm = TRUE),
              FN = sum(`False Negative`, na.rm = TRUE), .groups = "drop")

  # left join off ccs_index so all 130 categories show up, including any the
  # pipeline emitted nothing for
  per <- ccs_index %>%
    left_join(per, by = "CCS_ID") %>%
    mutate(across(c(`ICD-9 codes seen`, TP, FP, FN), ~ tidyr::replace_na(.x, 0L)))

  rows[[length(rows) + 1]] <- data.frame(
    `CCS ID` = per$CCS_ID,
    Category = per$description,
    `ICD-9 codes in the category` = per$n_codes,
    `Single code category` = ifelse(per$n_codes == 1, "yes", ""),
    Model = m,
    `Cosine cutoff` = b$threshold,
    `Top N co-occurring` = b$top_n,
    `Correct pairs to find` = per$TP + per$FN,
    `Mappings produced` = per$TP + per$FP,
    `True positives` = per$TP,
    `False positives` = per$FP,
    `False negatives` = per$FN,
    Precision = round(ifelse(per$TP + per$FP > 0, per$TP / (per$TP + per$FP), NA), 3),
    Recall    = round(ifelse(per$TP + per$FN > 0, per$TP / (per$TP + per$FN), NA), 3),
    F1        = round(ifelse(2 * per$TP + per$FP + per$FN > 0,
                             2 * per$TP / (2 * per$TP + per$FP + per$FN), NA), 3),
    check.names = FALSE)

  settings[[length(settings) + 1]] <- data.frame(
    Model = m, `Cosine cutoff` = b$threshold, `Top N co-occurring` = b$top_n,
    Rule = unname(SCENARIO[as.character(b$flag)]),
    `True positives` = b$tp, `False positives` = b$fp, `False negatives` = b$fn,
    `Overall F1` = b$f1, check.names = FALSE)
}

all_cat <- do.call(rbind, rows)
all_cat <- all_cat[order(all_cat$`CCS ID`,
                         match(all_cat$Model, names(MODELS))), ]
settings <- do.call(rbind, settings)

stopifnot(nrow(all_cat) == nrow(ccs_index) * length(MODELS))
for (m in names(MODELS)) {
  s <- all_cat[all_cat$Model == m, ]
  stopifnot(sum(s$`True positives`) == settings$`True positives`[settings$Model == m])
}

cb <- all_cat[all_cat$Model == "ClinicalBERT",
              setdiff(names(all_cat), c("Model", "Cosine cutoff", "Top N co-occurring"))]

## how lopsided the categories are. this is the concrete version of "too much
## data", and it is what says which categories can carry a number at all
size <- all_cat %>%
  filter(Model == "ClinicalBERT") %>%
  mutate(band = cut(`ICD-9 codes in the category`, c(0, 1, 2, 4, 9, Inf),
                    labels = c("1 code", "2 codes", "3 to 4 codes",
                               "5 to 9 codes", "10 or more codes"))) %>%
  group_by(band) %>%
  summarise(Categories = n(),
            `ICD-9 codes` = sum(`ICD-9 codes in the category`),
            `Correct pairs to find` = sum(`Correct pairs to find`),
            .groups = "drop") %>%
  as.data.frame()
size$`Share of categories` <- sprintf("%.0f%%", 100 * size$Categories / sum(size$Categories))
size$`Share of correct pairs` <- sprintf("%.0f%%", 100 * size$`Correct pairs to find` /
                                           sum(size$`Correct pairs to find`))
names(size)[1] <- "Category size"

cat("\nhow the categories are sized\n")
print(size, row.names = FALSE)
cat("\nsettings each model is scored at\n")
print(settings, row.names = FALSE)


hdr <- createStyle(textDecoration = "bold", valign = "bottom")
wb  <- createWorkbook()

add <- function(name, df, widths) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}
add("counts by category", all_cat,
    c(8, 46, 24, 19, 14, 14, 19, 20, 18, 15, 15, 16, 10, 9, 8))
add("ClinicalBERT only", cb, c(8, 46, 24, 19, 20, 18, 15, 15, 16, 10, 9, 8))
add("category size", size, c(18, 12, 13, 21, 19, 21))
add("settings", settings, c(14, 14, 19, 32, 15, 15, 16, 11))
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
