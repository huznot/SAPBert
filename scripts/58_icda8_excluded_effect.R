source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# what she asked for on 30 sep: the 52 icd-9 codes with no icda-8 match were
# reviewed (52_icda8_excluded_list.R, six look wrong), and now the icda-8
# estimates have to account for them the way the icd-9 to icd-10-ca ones
# already do.
#
# nothing in the repo ever dropped them: validate_mapping defaults to
# drop_unmatched = FALSE and no script overrides it. so they are already in
# every reported number. this quantifies what that is worth, and what the six
# corrections would change, by running three scenarios side by side:
#
#   codes dropped        the old behaviour, scored only on codes with a match
#   codes included       what the repo already does
#   six corrections      the six from the review treated as real matches

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("52_codes_accuracy_effect.xlsx")

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_8_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_8_9_sapbert_base_nocode.xlsx")
THRESHOLDS <- c(0.95, 0.99, 0.995, 0.999)   # the 8_9 track still uses the relative rule
TOP_NS     <- c(10, 20, 25)
RULES      <- 1:4
RULE_TEXT  <- c("1" = "top cosine or top co-occurrence",
                "2" = "top cosine or in both lists",
                "3" = "top co-occurrence or in both lists",
                "4" = "any of the three")

# the six the review marked as validation data wrong, where the icda-8 rubric is
# unambiguous and nothing competes. same six as 52_icda8_excluded_list.R
CORRECTIONS <- data.frame(
  icd9  = c("165", "175", "179", "555", "576", "745"),
  icda8 = c("163", "174", "182", "563", "576", "746"),
  note  = c("ICD-8 stops at 163, no other home for ill-defined respiratory sites",
            "ICD-8 did not split breast cancer by sex, 174 covers both",
            "uterus part unspecified has only one home, 182 other",
            "563 is the only chronic inflammatory bowel rubric in ICD-8",
            "same rubric under the same number",
            "746 is the only heart rubric, the numbers do not line up"),
  stringsAsFactors = FALSE)

man  <- read_excel(VAL, sheet = "Validaion_ICD9_ICD8")
excl <- read_excel(VAL, sheet = "Validation_ICD9_ICD8_Excld")
ccs  <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
cooc <- load_cooccurrence_df(file.path(ORIG, "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"))
a8   <- read_excel(LAB, sheet = "ICDA-8-3Level")
lab9 <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level")

e52 <- as.character(excl$`ICD-9-CM`)
stopifnot(length(e52) == 52, all(CORRECTIONS$icd9 %in% e52))

# scenario 3: the six move out of the excluded list and into the validated pairs
man_fix <- rbind(
  data.frame(`ICD-9-CM` = as.character(man$`ICD-9-CM`),
             `ICDA-8` = as.character(man$`ICDA-8`), check.names = FALSE),
  data.frame(`ICD-9-CM` = CORRECTIONS$icd9, `ICDA-8` = CORRECTIONS$icda8, check.names = FALSE))
excl_fix <- excl[!(as.character(excl$`ICD-9-CM`) %in% CORRECTIONS$icd9), ]

SCEN <- list(
  list(name = "codes dropped",   man = man,     excl = excl,     drop = TRUE),
  list(name = "codes included",  man = man,     excl = excl,     drop = FALSE),
  list(name = "six corrections", man = man_fix, excl = excl_fix, drop = FALSE))

cat(sprintf("%d validated pairs over %d icd-9 codes, %d codes with no match\n",
            nrow(man), length(unique(as.character(man$`ICD-9-CM`))), length(e52)))
cat(sprintf("with the six corrections: %d pairs, %d codes with no match\n\n",
            nrow(man_fix), nrow(excl_fix)))

co_by_n <- setNames(lapply(TOP_NS, function(tn)
  get_cooccurrence_codes_from_df(cooc, tn, "ICD_9_CM_Code", "ICDA_8_Code", "ICDA_8")),
  as.character(TOP_NS))

rows <- list()
for (m in names(MODELS)) {
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[m]]))
  for (thr in THRESHOLDS) {
    sim <- get_similarity_scores_from_sheets(sheets, thr, "ICDA_8")
    for (tn in TOP_NS) {
      mg <- merge_and_flag(sim, co_by_n[[as.character(tn)]], "ICDA_8",
                           find_icda8_chapter, chapter_alignment_8)
      for (rl in RULES) {
        au <- select_rows_by_flags(mg, rl)
        n52 <- sum(as.character(au$ICD_9_CM) %in% e52)
        for (sc in SCEN) {
          fin <- validate_mapping(sc$man, au, ccs, sc$excl, "ICDA-8", "ICDA_8",
                                  drop_unmatched = sc$drop)
          tp <- sum(fin$`True Positive`, na.rm = TRUE)
          fp <- sum(fin$`False Positive`, na.rm = TRUE)
          fn <- sum(fin$`False Negative`, na.rm = TRUE)
          rows[[length(rows) + 1]] <- data.frame(
            Model = m, Scenario = sc$name,
            `Similarity threshold` = thr, `Top N co-occurring` = tn,
            `Rule number` = rl, Rule = unname(RULE_TEXT[as.character(rl)]),
            `Mappings produced` = nrow(au),
            `Mappings for the 52 excluded codes` = n52,
            `Correct pairs to find` = tp + fn,
            `True positives` = tp, `False positives` = fp, `False negatives` = fn,
            Precision = round(ifelse(tp + fp > 0, tp / (tp + fp), NA), 3),
            Recall    = round(ifelse(tp + fn > 0, tp / (tp + fn), NA), 3),
            F1        = round(ifelse(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn), NA), 3),
            check.names = FALSE)
        }
      }
    }
  }
  cat("done", m, "\n")
}
res <- do.call(rbind, rows)
res <- res[order(res$Model, res$`Rule number`, res$`Similarity threshold`,
                 res$`Top N co-occurring`, match(res$Scenario, vapply(SCEN, `[[`, "", "name"))), ]

## what including them costs, at each model's best setting under each rule
eff <- do.call(rbind, lapply(split(res, list(res$Model, res$`Rule number`)), function(d) {
  inc <- d[d$Scenario == "codes included", ]
  b   <- inc[which.max(inc$F1), ]
  key <- d$`Similarity threshold` == b$`Similarity threshold` &
         d$`Top N co-occurring` == b$`Top N co-occurring`
  g <- function(s) d$F1[key & d$Scenario == s]
  data.frame(Model = b$Model, `Rule number` = b$`Rule number`, Rule = b$Rule,
             `Similarity threshold` = b$`Similarity threshold`,
             `Top N co-occurring` = b$`Top N co-occurring`,
             `Mappings for the 52 excluded codes` = b$`Mappings for the 52 excluded codes`,
             `F1 codes dropped` = g("codes dropped"),
             `F1 codes included` = g("codes included"),
             `F1 six corrections` = g("six corrections"),
             `Cost of including them` = round(g("codes included") - g("codes dropped"), 3),
             `Change from the corrections` = round(g("six corrections") - g("codes included"), 3),
             check.names = FALSE)
}))
eff <- eff[order(eff$Model, eff$`Rule number`), ]

cat("\neffect at each model's best setting, per rule\n")
print(eff[, c("Model", "Rule number", "Mappings for the 52 excluded codes",
              "F1 codes dropped", "F1 codes included", "F1 six corrections",
              "Cost of including them", "Change from the corrections")], row.names = FALSE)

a8lab <- setNames(as.character(a8$ICDA_8_LABEL), as.character(a8$ICDA_8))
n9    <- setNames(as.character(lab9$ICD_9_CM_LABEL), as.character(lab9$ICD_9_CM))
fixes <- data.frame(
  `ICD-9-CM` = CORRECTIONS$icd9,
  `ICD-9-CM label` = unname(n9[CORRECTIONS$icd9]),
  `Proposed ICDA-8` = CORRECTIONS$icda8,
  `ICDA-8 label` = unname(a8lab[CORRECTIONS$icda8]),
  `Why` = CORRECTIONS$note, check.names = FALSE)

hdr <- createStyle(textDecoration = "bold", valign = "bottom", fgFill = "#D9D9D9",
                   border = "bottom", borderStyle = "medium")
CB  <- createStyle(fgFill = "#EAF1FB"); SB <- createStyle(fgFill = "#FDF1E4")
wb  <- createWorkbook()
add <- function(name, df, widths, filter = FALSE, freeze_col = 1) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr, withFilter = filter)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = freeze_col)
  if ("Model" %in% names(df)) {
    for (p in list(list("ClinicalBERT", CB), list("SapBERT", SB))) {
      r <- which(df$Model == p[[1]]) + 1
      if (length(r)) addStyle(wb, name, p[[2]], rows = r, cols = seq_along(widths),
                              gridExpand = TRUE, stack = TRUE)
    }
  }
  for (cn in intersect(c("F1", "F1 codes included"), names(df)))
    conditionalFormatting(wb, name, cols = which(names(df) == cn), rows = 2:(nrow(df) + 1),
                          type = "colourScale", style = c("#F8696B", "#FFEB84", "#63BE7B"),
                          rule = c(0.5, 0.7, 0.85))
}
add("effect at best setting", eff, c(14, 8, 30, 12, 12, 16, 14, 15, 15, 14, 16), freeze_col = 2)
add("all settings", res, c(14, 16, 12, 12, 8, 30, 14, 16, 14, 12, 12, 13, 10, 9, 8),
    filter = TRUE, freeze_col = 3)
add("the six corrections", fixes, c(11, 42, 12, 46, 60))
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "with", nrow(res), "rows\n")
