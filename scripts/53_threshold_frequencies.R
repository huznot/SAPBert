source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# she asked for the actual numbers at each cutoff, not the precision and recall
# summaries. everything here is already in absolute_threshold_grid.csv, it has
# just never been laid out as counts. no new runs, nothing recomputed.

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
OUT  <- out_path("absolute_threshold_frequencies.xlsx")

res <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

# the size of the problem. one similarity matrix gives the two counts, the
# validation sheet gives how many pairs are actually correct
sheets <- load_similarity_sheets(
  file.path(GEN, "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx"))
n_icd9 <- length(unique(unlist(lapply(sheets, function(d) setdiff(names(d), names(d)[1])))))
n_t10  <- length(unique(unlist(lapply(sheets, function(d) as.character(d[[1]])))))
n_pairs <- n_icd9 * n_t10

man <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
n_correct <- nrow(man)

cat(sprintf("%d icd-9 codes x %d icd-10-ca codes = %s possible pairs\n",
            n_icd9, n_t10, format(n_pairs, big.mark = ",")))
cat(sprintf("%d of those pairs are marked correct in the validation data\n\n", n_correct))

stopifnot(all(res$sim_pairs_kept <= n_pairs))
stopifnot(all(res$tp + res$fp == res$codes_emitted))

MODEL_ORDER <- c("ClinicalBERT", "SapBERT", "mpnet")
res$model <- factor(res$model, levels = MODEL_ORDER)

pct <- function(x, d) sprintf("%.2f%%", 100 * x / d)

## sheet 2. the similarity step on its own, before co-occurrence is brought in.
## this is the "how many of all possible pairs does the cutoff call similar"
## question she asked directly
step <- unique(res[res$mode == "absolute",
                   c("model", "threshold", "sim_pairs_kept", "icd9_with_any_sim")])
step <- step[order(step$model, step$threshold), ]
step <- data.frame(
  Model = as.character(step$model),
  `Cosine cutoff` = step$threshold,
  `Possible pairs` = n_pairs,
  `Pairs at or above the cutoff` = step$sim_pairs_kept,
  `Share of all pairs` = pct(step$sim_pairs_kept, n_pairs),
  `Pairs dropped` = n_pairs - step$sim_pairs_kept,
  `ICD-9 codes with at least one` = step$icd9_with_any_sim,
  `ICD-9 codes with none` = n_icd9 - step$icd9_with_any_sim,
  check.names = FALSE)

cat("pairs surviving the cosine cutoff, before co-occurrence\n")
print(step, row.names = FALSE)

## sheet 3. the full grid as counts. one row per rule, which is a model, a
## cutoff, how deep the co-occurrence list goes and which of the four
## combination rules is used
counts <- function(d) {
  d <- d[order(d$model, d$threshold, d$top_n, d$flag), ]
  data.frame(
    Model = as.character(d$model),
    `Cosine cutoff` = d$threshold,
    `Top N co-occurring` = d$top_n,
    Rule = d$scenario,
    `Possible pairs` = n_pairs,
    `Pairs above the cutoff` = d$sim_pairs_kept,
    `Mappings produced` = d$codes_emitted,
    `Correct pairs to find` = d$tp + d$fn,
    `True positives` = d$tp,
    `False positives` = d$fp,
    `False negatives` = d$fn,
    Precision = d$precision, Recall = d$recall, F1 = d$f1, Accuracy = d$accuracy,
    check.names = FALSE)
}
abs_counts <- counts(res[res$mode == "absolute", ])
rel_counts <- counts(res[res$mode == "relative", ])

## sheet 4. the single best rule per model, so she does not have to hunt for it
best <- do.call(rbind, lapply(MODEL_ORDER, function(m) {
  d <- res[res$model == m & res$mode == "absolute", ]
  d[which.max(d$f1), ]
}))
best_t <- counts(best)

cat("\nbest rule per model, by f1\n")
print(best_t[, c("Model", "Cosine cutoff", "Top N co-occurring", "Mappings produced",
                 "True positives", "False positives", "False negatives", "F1")],
      row.names = FALSE)

## sheet 5. the plain cutoff against the old rule, in counts. she wanted to know
## whether the simpler cutoff costs anything and this is the answer as numbers
side <- do.call(rbind, lapply(MODEL_ORDER, function(m) {
  a <- res[res$model == m & res$mode == "absolute", ]
  r <- res[res$model == m & res$mode == "relative", ]
  a <- a[which.max(a$f1), ]; r <- r[which.max(r$f1), ]
  data.frame(
    Model = m,
    Rule = c("plain cutoff", "fraction of the column maximum"),
    Cutoff = c(a$threshold, r$threshold),
    `Top N co-occurring` = c(a$top_n, r$top_n),
    `Pairs above the cutoff` = c(a$sim_pairs_kept, r$sim_pairs_kept),
    `Mappings produced` = c(a$codes_emitted, r$codes_emitted),
    `True positives` = c(a$tp, r$tp),
    `False positives` = c(a$fp, r$fp),
    `False negatives` = c(a$fn, r$fn),
    F1 = c(a$f1, r$f1),
    check.names = FALSE)
}))

cat("\nplain cutoff against the old multiplier, both at their own best setting\n")
print(side, row.names = FALSE)

## sheet 1. what the columns mean, so the workbook stands on its own
readme <- data.frame(
  Item = c(
    "What this is",
    "ICD-9-CM codes",
    "ICD-10-CA codes",
    "Possible pairs",
    "Correct pairs to find",
    "",
    "Pairs above the cutoff",
    "Mappings produced",
    "True positives",
    "False positives",
    "False negatives",
    "",
    "Cosine cutoff",
    "Top N co-occurring",
    "Rule",
    "",
    "Source"),
  Detail = c(
    "counts behind the threshold results, for the ICD-9-CM to ICD-10-CA track",
    n_icd9,
    n_t10,
    sprintf("%s, every ICD-9 code against every ICD-10-CA code", format(n_pairs, big.mark = ",")),
    sprintf("%d, the pairs the validation data marks as correct", n_correct),
    "",
    "how many of the possible pairs score at or above the cosine cutoff",
    "how many mappings the pipeline outputs after co-occurrence is brought in",
    "mappings that match the validation data",
    "mappings that do not",
    "correct pairs the pipeline never produced",
    "",
    "a plain cosine score. a pair is kept if it scores at least this",
    "how far down the co-occurrence list the pipeline is allowed to look",
    "which combination of the cosine list and the co-occurrence list is used",
    "",
    "results/review/absolute_threshold_grid.csv, written by scripts/45_absolute_threshold_grid.R"),
  check.names = FALSE)

hdr <- createStyle(textDecoration = "bold", valign = "bottom")
wb <- createWorkbook()

add <- function(name, df, widths) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}

add("read me", readme, c(26, 88))
add("pairs above the cutoff", step, c(14, 14, 15, 26, 18, 14, 27, 22))
add("counts by rule", abs_counts, c(14, 14, 19, 32, 15, 21, 18, 20, 15, 15, 16, 11, 9, 8, 10))
add("best rule per model", best_t, c(14, 14, 19, 32, 15, 21, 18, 20, 15, 15, 16, 11, 9, 8, 10))
add("plain cutoff vs old rule", side, c(14, 30, 9, 19, 21, 18, 15, 15, 16, 8))
add("old rule, all counts", rel_counts, c(14, 14, 19, 32, 15, 21, 18, 20, 15, 15, 16, 11, 9, 8, 10))

saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
