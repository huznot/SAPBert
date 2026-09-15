source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# does pasting the code number onto the front of the label before embedding
# matter on the icd-9 to icda-8 track. it mattered a lot on the icd-9 to
# icd-10-ca track, and her expectation was that it would not matter here
# because both classifications are numeric so the numbers themselves look
# alike. both arms already ran in the full grid, so the grid comparison is
# read from results/grid/conditions/. the counts at the best setting are not in
# those files, so those six settings per track are rerun to get tp, fp and fn.

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("2_icda8_code_in_label.xlsx")

# the four pairs of arms that differ only by whether the code is in the text
PAIRS <- list(
  list(model = "ClinicalBERT", text = "none",
       code = "clinicalbert_base", nocode = "clinicalbert_base_nocode"),
  list(model = "SapBERT", text = "none",
       code = "sapbert_base", nocode = "sapbert_base_nocode"),
  list(model = "mpnet", text = "none",
       code = "mpnet_base", nocode = "mpnet_base_nocode"),
  list(model = "ClinicalBERT", text = "stopwords removed",
       code = "clinicalbert_stopwords", nocode = "clinicalbert_stopwords_nocode")
)

TRACKS <- list(
  `8_9` = list(
    label = "ICD-9-CM to ICDA-8", target_col_name = "ICDA_8",
    icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
    find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
    manual_target_col = "ICDA-8", manual_sheet = "Validaion_ICD9_ICD8",
    excl_sheet = "Validation_ICD9_ICD8_Excld",
    cooc = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx"),
  `10_9` = list(
    label = "ICD-9-CM to ICD-10-CA", target_col_name = "ICD_10_CA",
    icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
    find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
    manual_target_col = "ICD-10-CA", manual_sheet = "Validation_ICD9_ICD10",
    excl_sheet = "Validation_ICD9_ICD10_Excld",
    cooc = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx")
)

SCENARIO <- c("1" = "top cosine or top co-occurrence",
              "2" = "top cosine or in both lists",
              "3" = "top co-occurrence or in both lists",
              "4" = "any of the three")

ccs <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
for (tr in names(TRACKS)) {
  TRACKS[[tr]]$manual   <- read_excel(VAL, sheet = TRACKS[[tr]]$manual_sheet)
  TRACKS[[tr]]$excl     <- read_excel(VAL, sheet = TRACKS[[tr]]$excl_sheet)
  TRACKS[[tr]]$cooc_df  <- load_cooccurrence_df(file.path(ORIG, TRACKS[[tr]]$cooc))
}

# keyed by the condition tag, which is the filename. the model column cannot be
# used for this because sapbert_base is recorded as SapBERT-base-regen
gfiles <- list.files("results/grid/conditions", pattern = "[.]csv$", full.names = TRUE)
grid <- setNames(lapply(gfiles, read.csv, stringsAsFactors = FALSE),
                 sub("[.]csv$", "", basename(gfiles)))
stopifnot(all(unlist(lapply(PAIRS, `[`, c("code", "nocode"))) %in% names(grid)))

KEY <- c("track", "similarity_threshold", "top_n", "flag_combination")

## every grid point, the two arms side by side. 4 thresholds x 7 top_n x 4 rules
## is 112 points per model per track, so this is a paired comparison rather than
## a best against best, which one lucky setting could swing
side <- do.call(rbind, lapply(PAIRS, function(p) {
  a <- grid[[p$code]]; b <- grid[[p$nocode]]
  m <- merge(a[, c(KEY, "f1", "precision", "recall", "n_auto_mappings")],
             b[, c(KEY, "f1", "precision", "recall", "n_auto_mappings")],
             by = KEY, suffixes = c("_code", "_nocode"))
  data.frame(
    Track = vapply(m$track, function(t) TRACKS[[t]]$label, character(1)),
    Model = p$model, `Text cleaning` = p$text,
    `Cosine cutoff` = m$similarity_threshold,
    `Top N co-occurring` = m$top_n,
    Rule = unname(SCENARIO[as.character(m$flag_combination)]),
    `F1 with the code in the label` = m$f1_code,
    `F1 label only` = m$f1_nocode,
    `Difference` = round(m$f1_nocode - m$f1_code, 4),
    check.names = FALSE, row.names = NULL)
}))
side <- side[order(side$Track, side$Model, side$`Text cleaning`, side$`Cosine cutoff`,
                   side$`Top N co-occurring`), ]

## the summary of that: how often removing the code helped, and by how much
summ <- side %>%
  group_by(Track, Model, `Text cleaning`) %>%
  summarise(`Grid points` = n(),
            `Label only better` = sum(Difference > 0),
            Tied = sum(Difference == 0),
            `Code in label better` = sum(Difference < 0),
            `Mean difference in F1` = round(mean(Difference), 4),
            `Largest gain` = round(max(Difference), 4),
            `Largest loss` = round(min(Difference), 4),
            .groups = "drop") %>%
  as.data.frame()
summ <- summ[order(summ$Track, summ$Model, summ$`Text cleaning`), ]

cat("\npaired across every grid point, label only minus code in the label\n")
print(summ, row.names = FALSE)

## counts at the best setting for each arm. the grid files only carry precision,
## recall and f1, and she wants the raw numbers
run_one <- function(track, tag, thr, tn, fl) {
  tk <- TRACKS[[track]]
  tcn <- tk$target_col_name
  sheets <- load_similarity_sheets(
    file.path(GEN, sprintf("cosine_similarity_matrices_%s_%s.xlsx", track, tag)))
  long <- purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble(ICD_9_CM = as.character(c1), !!tcn := as.character(df[[id]]), Similarity = df[[c1]]))
  })
  long <- long[!is.na(long$Similarity), ]
  n_icd9   <- length(unique(long$ICD_9_CM))
  n_target <- length(unique(long[[tcn]]))

  sim <- get_similarity_scores_from_sheets(sheets, thr, tcn)
  co  <- get_cooccurrence_codes_from_df(tk$cooc_df, tn, tk$icd9_col, tk$target_col, tcn)
  mg  <- merge_and_flag(sim, co, tcn, tk$find_target_chapter_fn, tk$chapter_alignment)
  au  <- select_rows_by_flags(mg, fl)
  fin <- validate_mapping(tk$manual, au, ccs, tk$excl, tk$manual_target_col, tcn)
  mt  <- calculate_performance_metrics(fin)$overall
  list(possible = n_icd9 * n_target, n_icd9 = n_icd9, n_target = n_target,
       kept = nrow(sim), emitted = nrow(au),
       tp = sum(fin$`True Positive`, na.rm = TRUE),
       fp = sum(fin$`False Positive`, na.rm = TRUE),
       fn = sum(fin$`False Negative`, na.rm = TRUE),
       precision = mt$overall_precision, recall = mt$overall_recall, f1 = mt$overall_f1_score)
}

rows <- list()
for (tr in names(TRACKS)) {
  for (p in PAIRS) {
    for (arm in c("code", "nocode")) {
      tag <- p[[arm]]
      d <- grid[[tag]][grid[[tag]]$track == tr, ]
      b <- d[which.max(d$f1), ]
      cat(sprintf("%s / %s / %s  cutoff %.3f top %d rule %d\n", tr, p$model,
                  if (arm == "code") "code in label" else "label only",
                  b$similarity_threshold, b$top_n, b$flag_combination))
      r <- run_one(tr, tag, b$similarity_threshold, b$top_n, b$flag_combination)
      stopifnot(abs(r$f1 - b$f1) < 0.002)
      rows[[length(rows) + 1]] <- data.frame(
        Track = TRACKS[[tr]]$label, Model = p$model, `Text cleaning` = p$text,
        `Label content` = if (arm == "code") "code number then the label" else "the label only",
        `Cosine cutoff` = b$similarity_threshold,
        `Top N co-occurring` = b$top_n,
        Rule = unname(SCENARIO[as.character(b$flag_combination)]),
        `Possible pairs` = r$possible,
        `Pairs above the cutoff` = r$kept,
        `Mappings produced` = r$emitted,
        `Correct pairs to find` = r$tp + r$fn,
        `True positives` = r$tp, `False positives` = r$fp, `False negatives` = r$fn,
        Precision = r$precision, Recall = r$recall, F1 = r$f1,
        check.names = FALSE)
    }
  }
}
best <- do.call(rbind, rows)

cat("\ncounts at each arm's own best setting\n")
print(best[, c("Track", "Model", "Text cleaning", "Label content", "Mappings produced",
               "True positives", "False positives", "False negatives", "F1")],
      row.names = FALSE)

hdr <- createStyle(textDecoration = "bold", valign = "bottom")
wb  <- createWorkbook()

add <- function(name, df, widths) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}
add("counts at the best setting", best,
    c(22, 14, 18, 26, 14, 19, 32, 15, 21, 18, 20, 15, 15, 16, 11, 9, 8))
add("paired difference", summ, c(22, 14, 18, 12, 18, 7, 21, 22, 13, 13))
add("grid side by side", side, c(22, 14, 18, 14, 19, 32, 28, 14, 11))
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
