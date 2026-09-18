source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# head to head, clinicalbert against sapbert, after the 18 sep meeting. the
# cutoffs were set from each model's score distribution, not tuned for f1, and
# they are the only thing allowed to differ between the two models. top n and
# the rule are held the same for both, even where that is not the best setting
# for one of them, so the comparison is fair and shows where each one does
# badly. mpnet is out.
#
# rule 4 for both. under rule 1 the cutoff and top n barely change what comes
# out (it only takes the single top cosine and top co-occurring code), so the
# grid would look flat for clinicalbert.

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("clinicalbert_vs_sapbert.xlsx")

MODELS <- list(
  ClinicalBERT = list(file = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
                      cutoffs = c(0.85, 0.90, 0.95, 0.99)),
  SapBERT      = list(file = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
                      cutoffs = c(0.55, 0.60, 0.70, 0.80, 0.90, 0.99)))
TOP_NS <- c(10, 20, 25)
RULE   <- 4
RULE_TEXT <- "any of the three"

ccs_full <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>%
  mutate(ICD_9_CM = as.character(ICD_9_CM))
ccs_df <- ccs_full %>% select(ICD_9_CM, CCS_ID)
ccs_index <- ccs_full %>%
  group_by(CCS_ID) %>%
  summarise(description = first(CCS_CATEGORY_DESCRIPTION), n_codes = n(), .groups = "drop")

man  <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
excl <- read_excel(VAL, sheet = "Validation_ICD9_ICD10_Excld")
cooc <- load_cooccurrence_df(file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))
co_by_n <- setNames(lapply(TOP_NS, function(tn)
  get_cooccurrence_codes_from_df(cooc, tn, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")),
  as.character(TOP_NS))

# the settings that were already in the old grid are checked against it below
grid <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

prf <- function(tp, fp, fn) list(
  p  = round(ifelse(tp + fp > 0, tp / (tp + fp), NA), 3),
  r  = round(ifelse(tp + fn > 0, tp / (tp + fn), NA), 3),
  f1 = round(ifelse(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn), NA), 3))

overall <- list(); bycat <- list()
for (m in names(MODELS)) {
  cat("reading", m, "\n")
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[m]]$file))
  long <- purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble(ICD_9_CM = as.character(c1), ICD_10_CA = as.character(df[[id]]), Similarity = df[[c1]]))
  }) %>% filter(!is.na(Similarity))
  n_pairs <- length(unique(long$ICD_9_CM)) * length(unique(long$ICD_10_CA))

  for (thr in MODELS[[m]]$cutoffs) {
    sim <- long %>% filter(Similarity >= thr)
    for (tn in TOP_NS) {
      mg  <- merge_and_flag(sim, co_by_n[[as.character(tn)]], "ICD_10_CA",
                            find_icd10ca_chapter, chapter_alignment_10)
      au  <- select_rows_by_flags(mg, RULE)
      fin <- validate_mapping(man, au, ccs_df, excl, "ICD-10-CA", "ICD_10_CA")
      tp <- sum(fin$`True Positive`, na.rm = TRUE)
      fp <- sum(fin$`False Positive`, na.rm = TRUE)
      fn <- sum(fin$`False Negative`, na.rm = TRUE)
      stopifnot(tp + fn == nrow(man), tp + fp == nrow(au))

      old <- grid[grid$model == m & grid$mode == "absolute" & grid$flag == RULE &
                  abs(grid$threshold - thr) < 1e-9 & grid$top_n == tn, ]
      if (nrow(old)) stopifnot(old$tp == tp, old$fp == fp, old$fn == fn)

      s <- prf(tp, fp, fn)
      cat(sprintf("  %-12s %.2f top %2d  tp %3d fp %3d fn %3d  f1 %.3f%s\n",
                  m, thr, tn, tp, fp, fn, s$f1, if (nrow(old)) "  matches old grid" else ""))
      overall[[length(overall) + 1]] <- data.frame(
        Model = m, `Cosine cutoff` = thr, `Top N co-occurring` = tn, Rule = RULE_TEXT,
        `Possible pairs` = n_pairs, `Pairs above the cutoff` = nrow(sim),
        `Mappings produced` = tp + fp, `Correct pairs to find` = tp + fn,
        `True positives` = tp, `False positives` = fp, `False negatives` = fn,
        Precision = s$p, Recall = s$r, F1 = s$f1, check.names = FALSE)

      per <- fin %>% group_by(CCS_ID) %>%
        summarise(TP = sum(`True Positive`, na.rm = TRUE),
                  FP = sum(`False Positive`, na.rm = TRUE),
                  FN = sum(`False Negative`, na.rm = TRUE), .groups = "drop")
      per <- ccs_index %>% left_join(per, by = "CCS_ID") %>%
        mutate(across(c(TP, FP, FN), ~ tidyr::replace_na(.x, 0L)))
      stopifnot(sum(per$TP) == tp, sum(per$FP) == fp, sum(per$FN) == fn)
      c3 <- prf(per$TP, per$FP, per$FN)
      bycat[[length(bycat) + 1]] <- data.frame(
        `CCS ID` = per$CCS_ID, Category = per$description,
        `ICD-9 codes in the category` = per$n_codes,
        Model = m, `Cosine cutoff` = thr, `Top N co-occurring` = tn,
        `Correct pairs to find` = per$TP + per$FN, `Mappings produced` = per$TP + per$FP,
        `True positives` = per$TP, `False positives` = per$FP, `False negatives` = per$FN,
        Precision = c3$p, Recall = c3$r, F1 = c3$f1, check.names = FALSE)
    }
  }
}

overall <- do.call(rbind, overall)
bycat   <- do.call(rbind, bycat)
bycat   <- bycat[order(as.numeric(bycat$`CCS ID`), match(bycat$Model, names(MODELS)),
                       bycat$`Cosine cutoff`, bycat$`Top N co-occurring`), ]
stopifnot(nrow(bycat) == nrow(ccs_index) * nrow(overall))

# correct pairs per category do not depend on the model, so this is one table
size <- bycat[bycat$Model == "ClinicalBERT" & bycat$`Cosine cutoff` == 0.90 &
              bycat$`Top N co-occurring` == 25, ] %>%
  mutate(band = cut(`ICD-9 codes in the category`, c(0, 1, 2, 4, 9, Inf),
                    labels = c("1 code", "2 codes", "3 to 4 codes",
                               "5 to 9 codes", "10 or more codes"))) %>%
  group_by(band) %>%
  summarise(Categories = n(), `ICD-9 codes` = sum(`ICD-9 codes in the category`),
            `Correct pairs to find` = sum(`Correct pairs to find`), .groups = "drop") %>%
  as.data.frame()
size$`Share of categories`    <- sprintf("%.0f%%", 100 * size$Categories / sum(size$Categories))
size$`Share of correct pairs` <- sprintf("%.0f%%", 100 * size$`Correct pairs to find` /
                                           sum(size$`Correct pairs to find`))
names(size)[1] <- "Category size"

hdr <- createStyle(textDecoration = "bold", valign = "bottom")
wb  <- createWorkbook()
add <- function(name, df, widths, filter = FALSE) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr, withFilter = filter)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}
add("overall", overall, c(14, 14, 19, 17, 15, 21, 18, 20, 15, 15, 16, 10, 9, 8))
add("by category", bycat, c(8, 46, 24, 14, 14, 19, 20, 18, 15, 15, 16, 10, 9, 8), filter = TRUE)
add("category size", size, c(18, 12, 13, 21, 19, 21))
saveWorkbook(wb, OUT, overwrite = TRUE)

cat("\n"); print(overall[, c("Model", "Cosine cutoff", "Top N co-occurring",
                             "True positives", "False positives", "False negatives", "F1")],
                 row.names = FALSE)
cat("\nwrote", OUT, "with", nrow(overall), "settings and", nrow(bycat), "category rows\n")
