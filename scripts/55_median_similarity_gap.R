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
OUT  <- out_path("absolute_threshold_median_gap.xlsx")

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

readme <- data.frame(
  Item = c(
    "Question",
    "Short answer",
    "",
    "Why it happens",
    "",
    "",
    "Why a low median is the good sign",
    "",
    "",
    "The two groups",
    "pairs marked correct",
    "every other pair",
    "",
    "score distribution",
    "the gap",
    "what each cutoff admits",
    "",
    "Caveat",
    "",
    "Source"),
  Detail = c(
    "why is SapBERT's median cosine similarity so much lower than ClinicalBERT's",
    "because ClinicalBERT gives almost every pair of codes a high score, so its median is high and its cutoff has to be high too. SapBERT only scores a pair high when the two labels mean close to the same thing",
    "",
    "ClinicalBERT was trained on clinical notes, so it mostly learns that a piece of text is medical. Every ICD label is medical, so every pair looks alike to it and the scores bunch up near the top",
    "SapBERT was trained on pairs of medical terms that name the same thing, so its job is telling apart two terms that are both medical. That is the thing being asked of it here",
    "",
    "a low median means the pairs that should not match are being pushed down, which is what leaves room between them and the pairs that should",
    "the number to look at is the distance between the correct pairs and everything else, not the median on its own. see the gap sheet",
    "",
    "",
    "the 937 pairs in the validation data. these are the ones a cutoff should let through",
    "the other 720,515 of the 721,452 possible pairs. these are the ones a cutoff should keep out",
    "",
    "median and quartiles for each group, per model",
    "the two medians, the distance between them, and what each model's best cutoff lets through",
    "the same counts at every cutoff tested, so the tradeoff can be read off directly",
    "",
    "this is the cosine similarity step on its own. the co-occurrence step comes after it and removes some of what gets through, which is why the false positive counts in the other workbooks are far smaller than the other pairs counts here",
    "",
    "the label only matrices in data/generated/, and the validation sheet Validation_ICD9_ICD10"),
  check.names = FALSE)

hdr <- createStyle(textDecoration = "bold", valign = "bottom")
wb <- createWorkbook()
add <- function(name, df, widths) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}
add("read me", readme, c(30, 112))
add("the gap", summ, c(14, 20, 24, 26, 12, 21, 21, 19, 19))
add("score distribution", dist, c(14, 38, 11, 16, 16, 9, 16, 10, 8))
add("what each cutoff admits", adm, c(14, 14, 21, 21, 19, 19, 22))
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
