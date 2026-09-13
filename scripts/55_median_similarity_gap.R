source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# she asked in the meeting why SapBERT's median similarity is so much lower
# than ClinicalBERT's, and the answer given was a guess. this splits every one
# of the 721,452 possible pairs into the ones the validation data marks correct
# and everything else, and scores the two groups separately for each model.
# the median on its own is not the interesting number. the distance between the
# two groups is, and that is where ClinicalBERT is weakest.

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
OUT  <- out_path("3_similarity_gap_by_model.xlsx")

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
            mpnet        = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx")

CUTOFFS <- c(0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99)

man <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
key <- paste(as.character(man$`ICD-9-CM`), as.character(man$`ICD-10-CA`))

res <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)
best_cut <- vapply(names(MODELS), function(m) {
  d <- res[res$model == m & res$mode == "absolute", ]
  d$threshold[which.max(d$f1)]
}, numeric(1))

dist <- list(); adm <- list(); summ <- list()

for (m in names(MODELS)) {
  cat("reading", m, "\n")
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[m]]))
  long <- purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble(icd9 = as.character(c1), t10 = as.character(df[[id]]), s = df[[c1]]))
  })
  long <- long[!is.na(long$s), ]
  long$correct <- paste(long$icd9, long$t10) %in% key

  cor_s <- long$s[long$correct]
  oth_s <- long$s[!long$correct]
  cat(sprintf("  %s pairs, %d correct, %s other\n", format(nrow(long), big.mark = ","),
              length(cor_s), format(length(oth_s), big.mark = ",")))

  q <- function(x) as.numeric(quantile(x, c(0, .25, .5, .75, 1)))
  for (g in list(list("pairs the validation data marks correct", cor_s),
                 list("every other pair", oth_s))) {
    v <- q(g[[2]])
    dist[[length(dist) + 1]] <- data.frame(
      Model = m, Group = g[[1]], Pairs = length(g[[2]]),
      Lowest = round(v[1], 3), `25th percentile` = round(v[2], 3),
      Median = round(v[3], 3), `75th percentile` = round(v[4], 3),
      Highest = round(v[5], 3), Mean = round(mean(g[[2]]), 3),
      check.names = FALSE)
  }

  ## the table that makes the point. at every cutoff, how many correct pairs
  ## get through against how many of the ~720,000 that should not
  adm[[length(adm) + 1]] <- data.frame(
    Model = m, `Cosine cutoff` = CUTOFFS,
    `Correct pairs above it` = vapply(CUTOFFS, function(c) sum(cor_s >= c), integer(1)),
    `Share of correct pairs` = sprintf("%.1f%%", 100 * vapply(CUTOFFS,
        function(c) mean(cor_s >= c), numeric(1))),
    `Other pairs above it` = vapply(CUTOFFS, function(c) sum(oth_s >= c), integer(1)),
    `Share of other pairs` = sprintf("%.2f%%", 100 * vapply(CUTOFFS,
        function(c) mean(oth_s >= c), numeric(1))),
    `This model's best cutoff` = ifelse(CUTOFFS == best_cut[[m]], "yes", ""),
    check.names = FALSE)

  bc <- best_cut[[m]]
  summ[[length(summ) + 1]] <- data.frame(
    Model = m,
    `Median, correct pairs` = round(median(cor_s), 3),
    `Median, every other pair` = round(median(oth_s), 3),
    `Gap between the two medians` = round(median(cor_s) - median(oth_s), 3),
    `Best cutoff` = bc,
    `Correct pairs above it` = sum(cor_s >= bc),
    `Share of correct pairs` = sprintf("%.1f%%", 100 * mean(cor_s >= bc)),
    `Other pairs above it` = sum(oth_s >= bc),
    `Share of other pairs` = sprintf("%.2f%%", 100 * mean(oth_s >= bc)),
    check.names = FALSE)
}

dist <- do.call(rbind, dist)
adm  <- do.call(rbind, adm)
summ <- do.call(rbind, summ)

cat("\nscore distribution, the two groups apart\n")
print(dist, row.names = FALSE)
cat("\nthe gap, and what each model's own best cutoff lets through\n")
print(summ, row.names = FALSE)

hdr  <- createStyle(textDecoration = "bold", valign = "bottom")
note <- createStyle(fontColour = "#595959", textDecoration = "italic")
wb   <- createWorkbook()

# one line of context in A1, blank row, then the header on row 3. no notes tab
add <- function(name, df, widths, msg) {
  addWorksheet(wb, name)
  writeData(wb, name, msg, startRow = 1, startCol = 1)
  addStyle(wb, name, note, rows = 1, cols = 1)
  writeData(wb, name, df, startRow = 3, headerStyle = hdr)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 4)
}
add("the gap", summ, c(14, 20, 24, 26, 12, 21, 21, 19, 19),
    "Why SapBERT's median sits lower. The 937 correct pairs against the other 720,515. A low median means unrelated pairs are pushed down, which is what leaves room between the two groups.")
add("score distribution", dist, c(14, 38, 11, 16, 16, 9, 16, 10, 8),
    "Median and quartiles for each group, per model.")
add("what each cutoff admits", adm, c(14, 14, 21, 21, 19, 19, 22),
    "At each cutoff, correct pairs kept against wrong pairs kept. Cosine step only, before co-occurrence, so these counts are larger than the false positives elsewhere.")
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
