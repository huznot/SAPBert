source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(library(officer))

OUT <- out_path("cutoff_summary.docx")
res <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

models <- c("ClinicalBERT", "SapBERT", "mpnet")
bestf1 <- function(m, md, fl = NULL) {
  d <- res[res$model == m & res$mode == md, ]
  if (!is.null(fl)) d <- d[d$flag == fl, ]
  sprintf("%.3f", max(d$f1))
}

main <- data.frame(
  Model = models,
  `Old rule` = sapply(models, bestf1, md = "relative"),
  `Plain cutoff` = sapply(models, bestf1, md = "absolute"),
  check.names = FALSE, row.names = NULL)

scales <- data.frame(
  Model = models,
  `Median similarity` = c("0.841", "0.220", "0.161"),
  check.names = FALSE)

scen <- data.frame(
  Scenario = c("1  top cosine or top co-occurrence",
               "2  top cosine or in both lists",
               "3  top co-occurrence or in both lists",
               "4  any of the three"),
  `Old rule` = sapply(1:4, function(f) bestf1("SapBERT", "relative", f)),
  `Plain cutoff` = sapply(1:4, function(f) bestf1("SapBERT", "absolute", f)),
  check.names = FALSE)

p   <- function(x, t) body_add_par(x, t, style = "Normal")
h   <- function(x, t) body_add_par(x, t, style = "heading 2")
tbl <- function(x, d) body_add_par(
  body_add_table(x, d, style = "table_template", first_row = TRUE), "", style = "Normal")

doc <- read_docx()
doc <- body_add_par(doc, "Cosine cutoff test", style = "heading 1")
doc <- p(doc, "The old threshold was a fraction of each column's own best score, not a similarity. I swapped it for a plain cutoff and reran everything.")

doc <- h(doc, "Did it help?")
doc <- tbl(doc, main)
doc <- p(doc, "No. All three are slightly worse, by less than a point of F1, which is noise on this validation set. The reason to switch is that 0.995 was never a similarity, so we could not explain it in the paper.")

doc <- h(doc, "What is worth talking about")
doc <- tbl(doc, scales)
doc <- p(doc, "A cutoff of 0.80 keeps 551,755 pairs on ClinicalBERT and 324 on SapBERT. There is no single cutoff that works across models, so any cutoff has to be reported per model.")

doc <- h(doc, "The four scenarios, SapBERT")
doc <- tbl(doc, scen)
doc <- p(doc, "Scenario 4 is still the best. Scenario 2 is the only one the plain cutoff clearly improves.")

doc <- h(doc, "Separate finding")
doc <- p(doc, "Taking the ICD code out of the embedding input helps SapBERT and mpnet too, not just ClinicalBERT. SapBERT goes from 0.524 to 0.552. That is the best result we have.")

doc <- h(doc, "Two things I want to ask you")
doc <- p(doc, "1. Should the cosine side hand over more than one code? Only the single top scorer goes through now, and that is where recall is going.")
doc <- p(doc, "2. The gaps between models are 0.005 to 0.012. Do we need a bootstrap before saying SapBERT beats mpnet?")

doc <- h(doc, "If you want to see it")
doc <- p(doc, "All 1344 runs: results/review/absolute_threshold_grid.csv")
doc <- p(doc, "The nine codes: results/review/9_codes_review_UPDATED.xlsx")
doc <- p(doc, "249 and 239 in detail: results/review/e13_d48_review.xlsx")
doc <- p(doc, "Proof of the code in the label problem: scripts/41_identical_labels.R")
doc <- p(doc, "Longer writeup with the full sweep: results/review/cutoff_report.docx")

print(doc, target = OUT)
cat("wrote", OUT, "\n")
