source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(library(officer))

OUT <- out_path("cutoff_report.docx")
res <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

T <- c(0.50,0.60,0.70,0.75,0.80,0.85,0.90,0.92,0.95,0.97,0.98,0.99)
models <- c("ClinicalBERT","SapBERT","mpnet")

best_at <- function(mdl, thr) {
  d <- res[res$mode == "absolute" & res$model == mdl & res$threshold == thr, ]
  d[which.max(d$f1), ]
}

sweep <- data.frame(Cutoff = sprintf("%.2f", T))
for (m in models) {
  f <- sapply(T, function(t) sprintf("%.3f", best_at(m, t)$f1))
  cv <- sapply(T, function(t) best_at(m, t)$icd9_with_any_sim)
  sweep[[paste(m, "F1")]] <- f
  sweep[[paste(m, "codes")]] <- cv
}

scen_lab <- c("1. Top cosine or top co-occurrence",
              "2. Top cosine or in both lists",
              "3. Top co-occurrence or in both lists",
              "4. Any of the three")
scen <- do.call(rbind, lapply(1:4, function(fl) {
  do.call(rbind, lapply(models, function(m) {
    a <- res[res$mode == "absolute" & res$model == m & res$flag == fl, ]
    b <- res[res$mode == "relative" & res$model == m & res$flag == fl, ]
    a <- a[which.max(a$f1), ]; b <- b[which.max(b$f1), ]
    data.frame(Scenario = scen_lab[fl], Model = m,
               `Absolute F1` = sprintf("%.3f", a$f1),
               `Relative F1` = sprintf("%.3f", b$f1),
               Precision = sprintf("%.3f", a$precision),
               Recall = sprintf("%.3f", a$recall),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}))

hp <- do.call(rbind, lapply(models, function(m) {
  d <- res[res$mode == "absolute" & res$model == m, ]
  d <- d[which.max(d$precision), ]
  data.frame(Model = m, Cutoff = sprintf("%.2f", d$threshold), `Top N` = d$top_n,
             Scenario = d$flag, Precision = sprintf("%.3f", d$precision),
             Recall = sprintf("%.3f", d$recall), F1 = sprintf("%.3f", d$f1),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

bestrow <- do.call(rbind, lapply(models, function(m) {
  a <- res[res$mode == "absolute" & res$model == m, ]; a <- a[which.max(a$f1), ]
  b <- res[res$mode == "relative" & res$model == m, ]; b <- b[which.max(b$f1), ]
  data.frame(Model = m,
             `Best absolute F1` = sprintf("%.3f", a$f1),
             `at cutoff` = sprintf("%.2f", a$threshold),
             `Best relative F1` = sprintf("%.3f", b$f1),
             `at multiplier` = sprintf("%.3f", b$threshold),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

p <- function(x, txt, style = "Normal") body_add_par(x, txt, style = style)
tbl <- function(x, d) {
  x <- body_add_table(x, d, style = "Table Professional", first_row = TRUE)
  body_add_par(x, "", style = "Normal")
}

doc <- read_docx()
doc <- body_add_par(doc, "Cosine cutoff test", style = "heading 1")
doc <- p(doc, "Muhammad Irfan. ICD-9-CM to ICD-10-CA crosswalk.")
doc <- p(doc, "")

doc <- body_add_par(doc, "What I did", style = "heading 2")
doc <- p(doc, "The pipeline never used a plain similarity threshold. For every ICD-9 code it worked out the highest similarity score in that code's column and kept anything above a fraction of it. So a threshold of 0.995 did not mean 99.5 percent similarity. It meant 99.5 percent of whatever that column happened to top out at.")
doc <- p(doc, "I replaced that with a plain cutoff. Keep a pair if its similarity is above a fixed number, the same number for every code. I tested twelve cutoffs from 0.50 to 0.99. I also reran the old relative rule on the same matrices so the two can be compared fairly. Everything else in the pipeline is unchanged.")
doc <- p(doc, "All three models were rebuilt with the ICD code taken out of the text before embedding. That fix had only been applied to ClinicalBERT before.")
doc <- p(doc, "Each cutoff was run against 7 co-occurrence depths and Lynn's 4 scenarios, for 3 models. That is 1344 runs. They are all in results/review/absolute_threshold_grid.csv.")

doc <- body_add_par(doc, "Check that the code is right", style = "heading 2")
doc <- p(doc, "My relative run gives exactly the same numbers as the existing full_grid_best.csv for ClinicalBERT. Threshold 0.995, top 25, scenario 4, precision 0.623, recall 0.366, F1 0.461. Same in both. So the new script is not changing anything on its own.")

doc <- body_add_par(doc, "Main result", style = "heading 2")
doc <- p(doc, "The plain cutoff does not do better than the old relative rule. It is slightly worse for all three models, by less than one point of F1, which is too small to mean anything on this validation set.")
doc <- tbl(doc, bestrow)
doc <- p(doc, "So if the reason to switch is better accuracy, the numbers do not support it. The reason to switch is that the number is honest. With the old rule you cannot tell a reader what 0.995 means, because it is not a similarity. With a plain cutoff you can.")

doc <- body_add_par(doc, "The three models are on different scales", style = "heading 2")
doc <- p(doc, "This is the most useful thing that came out of the test. The models do not put their scores in the same range at all.")
doc <- tbl(doc, data.frame(
  Model = models,
  `Median similarity` = c("0.841","0.220","0.161"),
  `Lowest` = c("0.505","-0.176","-0.179"),
  `Highest` = c("1.000","1.000","1.000"),
  check.names = FALSE))
doc <- p(doc, "ClinicalBERT rates the average unrelated pair at 0.84. SapBERT rates it at 0.22. So a cutoff of 0.80 does almost nothing on ClinicalBERT, it still keeps 551755 of 721452 pairs. The same 0.80 on SapBERT keeps 324 pairs and leaves a third of all ICD-9 codes with no similarity candidate at all.")
doc <- p(doc, "This means there is no single cutoff number that works across models. Any cutoff has to be reported per model with the score range next to it.")
doc <- p(doc, "It also means ClinicalBERT is not separating clinical meaning very well. On ICD-9 249 there are 191 codes packed inside 0.036 of each other, including secondary parkinsonism, ascorbic acid deficiency and bartonellosis. The old relative rule hid this because it graded each column on a curve.")

doc <- body_add_par(doc, "Cutoff sweep", style = "heading 2")
doc <- p(doc, "Best F1 at each cutoff, taking the best co-occurrence depth and scenario at that cutoff. The codes column is how many of the 354 ICD-9 codes still get at least one similarity candidate.")
doc <- tbl(doc, sweep)
doc <- p(doc, "ClinicalBERT is flat from 0.50 to 0.92 because a cutoff in that range barely removes anything from its squashed range. That flatness is not the model being stable, it is the cutoff not doing any work.")
doc <- p(doc, "SapBERT and mpnet peak low, at 0.60 and 0.75, then fall off. Raising the cutoff past that point buys a little precision and loses a lot more recall.")
doc <- p(doc, "The codes column is the real difference between the two rules. Under the old relative rule this number is always 354, because the top scoring code in a column always clears a bar set as a fraction of itself. The similarity branch could never say no. With a plain cutoff it can.")

doc <- body_add_par(doc, "Lynn's four scenarios", style = "heading 2")
doc <- p(doc, "Same four scenarios, unchanged. Each one rerun under both rules. The F1 shown is the best that scenario can reach at any cutoff and any co-occurrence depth. Precision and recall are from the absolute run.")
doc <- tbl(doc, scen)
doc <- p(doc, "Scenario 4 is still the best for SapBERT and mpnet under either rule, and scenario 1 is close behind. Scenario 3 is the weakest every time. Dropping the cosine branch's own top pick costs more than the agreement between branches makes back.")
doc <- p(doc, "Scenario 2 is the one place the plain cutoff clearly wins, for all three models. SapBERT goes from 0.514 to 0.539. That scenario only emits a code the cosine branch ranked first or that both branches picked independently, so it gains the most from a cutoff that can actually throw out a weak column.")

doc <- body_add_par(doc, "A high precision option", style = "heading 2")
doc <- p(doc, "Chasing the best F1 hides a setting that may be more useful in practice. With scenario 2 and a short co-occurrence list, a plain cutoff reaches precision the old rule cannot.")
doc <- tbl(doc, hp)
doc <- p(doc, "SapBERT at a 0.90 cutoff with top 3 co-occurrence gets 0.994 precision at 0.181 recall. It maps about one code in five and is almost never wrong when it does.")
doc <- p(doc, "If the goal is to auto map what can be auto mapped and send the rest to a human coder, that is a better setting than the F1 peak. It is only reachable with a plain cutoff. Worth deciding which of the two goals the paper is aiming at, because they pick different settings.")

doc <- body_add_par(doc, "Separate finding", style = "heading 2")
doc <- p(doc, "Taking the ICD code out of the embedding input helps SapBERT and mpnet too, not just ClinicalBERT. Under the unchanged relative rule, SapBERT goes from 0.524 to 0.552 and mpnet from 0.527 to 0.540. SapBERT at 0.552 is the best F1 on this track so far.")

doc <- body_add_par(doc, "Still open", style = "heading 2")
doc <- p(doc, "Only the single highest scoring surviving code is ever handed on by the cosine branch. Emitting the top few instead has not been tested and is the change most likely to help recall, which is where all three models are weakest.")
doc <- p(doc, "The gaps between models are 0.005 to 0.012 in F1 and I am treating them as an ordering. With 937 validated pairs and no confidence intervals that is not solid yet. A bootstrap over ICD-9 codes would settle whether SapBERT really beats mpnet.")
doc <- p(doc, "The co-occurrence branch was left alone. At high cutoffs it is doing more of the work than before, not less, so its share of the result has grown and has not been checked separately.")

doc <- body_add_par(doc, "Files", style = "heading 2")
doc <- p(doc, "scripts/45_absolute_threshold_grid.R runs the whole thing.")
doc <- p(doc, "results/review/absolute_threshold_grid.csv has all 1344 runs.")
doc <- p(doc, "results/review/9_codes_review_UPDATED.xlsx is the reviewed unmatched codes.")
doc <- p(doc, "results/review/e13_d48_review.xlsx is the per code look at 249 and 239.")

print(doc, target = OUT)
cat("wrote", OUT, "\n")
