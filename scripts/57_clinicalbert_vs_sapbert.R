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
RULES  <- 1:4
RULE_TEXT <- c("1" = "top cosine or top co-occurrence",
               "2" = "top cosine or in both lists",
               "3" = "top co-occurrence or in both lists",
               "4" = "any of the three")
# all four rules, in the overall tab and then one category tab each, so a
# category can be read under whichever rule gets picked

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
      # the merge is the slow part and does not depend on the rule, so it is
      # done once and the four rules are applied to it
      mg <- merge_and_flag(sim, co_by_n[[as.character(tn)]], "ICD_10_CA",
                           find_icd10ca_chapter, chapter_alignment_10)
      for (rl in RULES) {
        au  <- select_rows_by_flags(mg, rl)
        fin <- validate_mapping(man, au, ccs_df, excl, "ICD-10-CA", "ICD_10_CA")
        tp <- sum(fin$`True Positive`, na.rm = TRUE)
        fp <- sum(fin$`False Positive`, na.rm = TRUE)
        fn <- sum(fin$`False Negative`, na.rm = TRUE)
        stopifnot(tp + fn == nrow(man), tp + fp == nrow(au))

        old <- grid[grid$model == m & grid$mode == "absolute" & grid$flag == rl &
                    abs(grid$threshold - thr) < 1e-9 & grid$top_n == tn, ]
        if (nrow(old)) stopifnot(old$tp == tp, old$fp == fp, old$fn == fn)

        s <- prf(tp, fp, fn)
        cat(sprintf("  %-12s %.2f top %2d rule %d  tp %3d fp %4d fn %3d  f1 %.3f%s\n",
                    m, thr, tn, rl, tp, fp, fn, s$f1, if (nrow(old)) "  ok" else ""))
        overall[[length(overall) + 1]] <- data.frame(
          Model = m, `Cosine cutoff` = thr, `Top N co-occurring` = tn,
          `Rule number` = rl, Rule = unname(RULE_TEXT[as.character(rl)]),
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
          Rule = unname(RULE_TEXT[as.character(rl)]),
          `Correct pairs to find` = per$TP + per$FN, `Mappings produced` = per$TP + per$FP,
          `True positives` = per$TP, `False positives` = per$FP, `False negatives` = per$FN,
          Precision = c3$p, Recall = c3$r, F1 = c3$f1, check.names = FALSE)
      }
    }
  }
}

overall <- do.call(rbind, overall)
bycat   <- do.call(rbind, bycat)
bycat   <- bycat[order(as.numeric(bycat$`CCS ID`), match(bycat$Model, names(MODELS)),
                       bycat$`Cosine cutoff`, bycat$`Top N co-occurring`), ]
stopifnot(nrow(bycat) == nrow(overall) * nrow(ccs_index))
bycat_by_rule <- split(bycat, match(bycat$Rule, RULE_TEXT))

# correct pairs per category do not depend on the model, so this is one table
size <- bycat[bycat$Model == "ClinicalBERT" & bycat$`Cosine cutoff` == 0.90 &
              bycat$`Top N co-occurring` == 25 & bycat$Rule == RULE_TEXT[["4"]], ] %>%
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

## one row per category: the best each model manages anywhere in the grid, and
## the setting that produced it. this is the tab to present from, 130 rows
## instead of 3,900. she asked that every table say which cutoff and top n it
## is using, so the winning setting travels with each score
# where several settings tie on f1, pick one deterministically rather than
# whichever row happened to come first: fewest false positives, then the
# smallest cutoff, top n and rule. otherwise the setting shown moves around
win <- bycat %>%
  mutate(.rule = match(Rule, RULE_TEXT)) %>%
  arrange(`CCS ID`, Model, desc(F1), `False positives`,
          `Cosine cutoff`, `Top N co-occurring`, .rule) %>%
  group_by(`CCS ID`, Model) %>%
  slice(1) %>%
  ungroup()

base <- bycat %>%
  group_by(`CCS ID`, Category, `ICD-9 codes in the category`) %>%
  summarise(`Correct pairs to find` = first(`Correct pairs to find`), .groups = "drop") %>%
  as.data.frame()

grab <- function(mdl, col) {
  w <- win[win$Model == mdl, ]
  w[[col]][match(base$`CCS ID`, w$`CCS ID`)]
}
summ <- data.frame(
  `CCS ID` = base$`CCS ID`, Category = base$Category,
  `ICD-9 codes in the category` = base$`ICD-9 codes in the category`,
  `Correct pairs to find` = base$`Correct pairs to find`,
  `Best F1 ClinicalBERT`  = grab("ClinicalBERT", "F1"),
  `ClinicalBERT cutoff`   = grab("ClinicalBERT", "Cosine cutoff"),
  `ClinicalBERT top N`    = grab("ClinicalBERT", "Top N co-occurring"),
  `ClinicalBERT rule`     = grab("ClinicalBERT", "Rule"),
  `Best F1 SapBERT`       = grab("SapBERT", "F1"),
  `SapBERT cutoff`        = grab("SapBERT", "Cosine cutoff"),
  `SapBERT top N`         = grab("SapBERT", "Top N co-occurring"),
  `SapBERT rule`          = grab("SapBERT", "Rule"),
  check.names = FALSE)
summ$`Better model` <- ifelse(summ$`Best F1 SapBERT` > summ$`Best F1 ClinicalBERT`, "SapBERT",
                       ifelse(summ$`Best F1 SapBERT` < summ$`Best F1 ClinicalBERT`, "ClinicalBERT", "tied"))
summ$`Best either`  <- pmax(summ$`Best F1 ClinicalBERT`, summ$`Best F1 SapBERT`)
summ <- summ[order(summ$`Best either`, -summ$`ICD-9 codes in the category`), ]

hdr  <- createStyle(textDecoration = "bold", valign = "bottom", fgFill = "#D9D9D9",
                    border = "bottom", borderStyle = "medium")
CB   <- createStyle(fgFill = "#EAF1FB")   # clinicalbert rows
SB   <- createStyle(fgFill = "#FDF1E4")   # sapbert rows
SEP  <- createStyle(border = "top", borderColour = "#808080", borderStyle = "thin")
wb   <- createWorkbook()

# red to green across the F1 column, so the weak categories stand out while she
# scrolls rather than having to read every number
scale3 <- function(sheet, col, n)
  conditionalFormatting(wb, sheet, cols = col, rows = 2:(n + 1), type = "colourScale",
                        style = c("#F8696B", "#FFEB84", "#63BE7B"), rule = c(0, 0.4, 0.8))

add <- function(name, df, widths, filter = FALSE, freeze_col = 1, band_model = TRUE) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr, withFilter = filter)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = freeze_col)
  n <- nrow(df)
  if (band_model && "Model" %in% names(df)) {
    cb <- which(df$Model == "ClinicalBERT") + 1
    sb <- which(df$Model == "SapBERT") + 1
    if (length(cb)) addStyle(wb, name, CB, rows = cb, cols = seq_along(widths),
                             gridExpand = TRUE, stack = TRUE)
    if (length(sb)) addStyle(wb, name, SB, rows = sb, cols = seq_along(widths),
                             gridExpand = TRUE, stack = TRUE)
  }
  # a line between categories so the blocks are visible at a glance
  if ("CCS ID" %in% names(df)) {
    brk <- which(c(TRUE, df$`CCS ID`[-1] != df$`CCS ID`[-n])) + 1
    addStyle(wb, name, SEP, rows = brk, cols = seq_along(widths),
             gridExpand = TRUE, stack = TRUE)
  }
  for (cn in intersect(c("F1", "Best F1 ClinicalBERT", "Best F1 SapBERT", "Best either"), names(df)))
    scale3(name, which(names(df) == cn), n)
}

add("overall", overall, c(14, 12, 12, 8, 30, 14, 15, 14, 14, 12, 12, 13, 10, 9, 8),
    filter = TRUE, freeze_col = 2)
add("category summary", summ, c(7, 40, 13, 13, 17, 15, 14, 30, 15, 13, 12, 30, 12, 10),
    filter = TRUE, freeze_col = 3)
for (rl in RULES) {
  d <- bycat_by_rule[[as.character(rl)]]
  d$Rule <- NULL   # the tab name already says which rule it is
  add(sprintf("rule %d by category", rl), d,
      c(7, 42, 13, 14, 12, 12, 14, 14, 12, 12, 13, 10, 9, 8),
      filter = TRUE, freeze_col = 4)
}
add("category size", size, c(18, 12, 13, 21, 19, 21))
saveWorkbook(wb, OUT, overwrite = TRUE)

cat("\n"); best <- do.call(rbind, lapply(split(overall, list(overall$Model, overall$`Rule number`)),
  function(d) d[which.max(d$F1), ]))
print(best[order(best$Model, best$`Rule number`),
           c("Model", "Rule number", "Cosine cutoff", "Top N co-occurring",
             "True positives", "False positives", "False negatives", "F1")],
      row.names = FALSE)
cat("\nwrote", OUT, "with", nrow(overall), "settings and", nrow(bycat), "category rows\n")
