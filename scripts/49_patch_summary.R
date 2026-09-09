source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(library(officer))

# rebuilds the summary from the wording exactly as it stands in the hand edited
# doc. word splits text across runs once it saves, so body_replace_all_text
# cannot find anything and in place editing is not an option. every paragraph
# below is transcribed from the current document. if the doc is edited by hand
# again, transcribe the new wording here before running this.

DOC <- out_path("cutoff_summary.docx")

## tables carried over unchanged from the current document ------------------
main <- data.frame(
  Model = rep(c("ClinicalBERT","SapBERT","mpnet"), each = 2),
  Rule = rep(c("old rule","plain cutoff"), 3),
  Precision = c("0.623","0.590","0.732","0.658","0.711","0.674"),
  Recall    = c("0.366","0.368","0.443","0.464","0.435","0.443"),
  F1        = c("0.461","0.453","0.552","0.544","0.540","0.535"),
  Accuracy  = c("0.300","0.293","0.381","0.374","0.370","0.365"),
  check.names = FALSE)

scales <- data.frame(
  Model = c("ClinicalBERT","SapBERT","mpnet"),
  `Median similarity` = c("0.841","0.220","0.161"),
  check.names = FALSE)

scen <- data.frame(
  `Which codes get output` = rep(c(
    "1   the best cosine match, plus the most co-occurring code",
    "2   the best cosine match, plus codes both sides picked",
    "3   the most co-occurring code, plus codes both sides picked",
    "4   all of the above"), each = 2),
  Rule = rep(c("old rule","plain cutoff"), 4),
  Precision = c("0.616","0.590","0.597","0.501","0.659","0.666","0.623","0.511"),
  Recall    = c("0.362","0.368","0.313","0.391","0.322","0.296","0.366","0.393"),
  F1        = c("0.456","0.453","0.411","0.439","0.433","0.410","0.461","0.444"),
  Accuracy  = c("0.295","0.293","0.258","0.281","0.276","0.257","0.300","0.285"),
  check.names = FALSE)

src <- data.frame(
  Model = c(rep("ClinicalBERT",3), rep("SapBERT",4), rep("mpnet",4)),
  Source = c("Top of both branches","Top of co-occurrence only","Top of cosine only",
             "Top of both branches","Top of co-occurrence only","Top of cosine only",
             "In both lists, top of neither",
             "Top of both branches","Top of co-occurrence only","Top of cosine only",
             "In both lists, top of neither"),
  Correct = c(102,124,119, 136,85,175,39, 129,96,158,32),
  Wrong   = c(4,112,124, 4,111,28,83, 1,115,21,64),
  `Hit rate` = c("96%","53%","49%", "97%","43%","86%","32%", "99%","45%","88%","33%"),
  check.names = FALSE)

rank <- data.frame(
  Model = c("ClinicalBERT","mpnet","SapBERT"),
  `1` = c(202,310,311), `2 to 5` = c(45,23,25), `6 to 20` = c(32,9,3),
  `over 20` = c(66,3,6), `not in column` = c(0,0,0),
  `rank 1 share` = c("59%","90%","90%"), check.names = FALSE)

## the two tables being corrected --------------------------------------------
# co-occurrence depth, rebuilt with the scenario held at 4 for all three models.
# the version in the doc let each model sit on its own best scenario, and
# clinicalbert's is scenario 1, which only ever takes the single top
# co-occurring code, so top n could not move it and the row looked flat
depth <- read.csv(out_path("trend_cooccurrence_depth.csv"), check.names = FALSE,
                  stringsAsFactors = FALSE)

# correct codes per icd-9 code, now over all 354 codes rather than the 345 that
# have a target. the 9 excluded codes were reviewed: 4 were validation errors
# and do have a code, 5 genuinely have none
targets <- data.frame(
  `Correct codes for one ICD-9 code` = c("0","1","2","3","4 or more"),
  `ICD-9 codes` = c(5,131,76,52,90),
  Share = sprintf("%.0f%%", 100 * c(5,131,76,52,90) / 354),
  check.names = FALSE)

## ---------------------------------------------------------------------------
p   <- function(x, t) body_add_par(x, t, style = "Normal")
h   <- function(x, t) body_add_fpar(x, fpar(ftext(t, fp_text(bold = TRUE, font.size = 13)),
                                          fp_p = fp_par(padding.top = 10, padding.bottom = 2)))
tbl <- function(x, d) body_add_par(
  body_add_table(x, d, style = "table_template", first_row = TRUE), "", style = "Normal")

doc <- read_docx()
doc <- body_add_fpar(doc, fpar(ftext("Cosine cutoff test", fp_text(bold = TRUE, font.size = 16)),
                               fp_p = fp_par(padding.bottom = 6)))
doc <- p(doc, "The pipeline picks ICD-10-CA codes for an ICD-9 code two ways. Cosine similarity between the two labels, and how often the two codes appear together in the data. The threshold sits on the cosine side and filters which codes are allowed through to the next step.")
doc <- p(doc, "For each ICD-9 code, the old threshold took the highest similarity in that code's column and kept anything above a fraction of it. At 0.995 that means 99.5 percent of whatever that column happened to top out at. So the bar moved for every code. A code whose best match scored 0.99 was judged against 0.985. A code whose best match scored 0.62 was judged against 0.617.")
doc <- p(doc, "I replaced it with a plain cutoff. Keep a pair if the similarity is above a fixed number, the same number for every ICD-9 code. I tested twelve values from 0.50 to 0.99, and reran the old rule on the same data to compare.")

doc <- h(doc, "Overall performance")
doc <- p(doc, "Looking at best F1 for each row")
doc <- tbl(doc, main)
doc <- p(doc, "The plain cutoff is slightly worse on F1 for all three models, by less than a point. The plain cutoff gives up precision and buys back a little recall, and it does this the same way on all three models. SapBERT goes from 0.732 precision and 0.443 recall to 0.658 and 0.464. ")
doc <- p(doc, "The plain cutoff values behind the table are 0.90 for ClinicalBERT, 0.60 for SapBERT and 0.75 for mpnet.")

doc <- h(doc, "Median Similarities")
doc <- tbl(doc, scales)
doc <- p(doc, "A cutoff of 0.80 keeps 551,755 pairs on ClinicalBERT and 324 on SapBERT. The best cutoff is 0.90 for ClinicalBERT, 0.60 for SapBERT and 0.75 for mpnet. ")

doc <- h(doc, "The four scenarios, ClinicalBERT")
doc <- p(doc, "The scenarios are the four rules for deciding which codes come out at the end. ")
doc <- tbl(doc, scen)

doc <- h(doc, "Co-occurrence depth")
doc <- p(doc, "Top N is how many co-occurrence candidates the pipeline is allowed to consider. For each ICD-9 code it takes the N ICD-10-CA codes that most often appear alongside it, and only those can be picked. Only Top N changes down the rows. The cutoff stays at each model's best value and the scenario is fixed at 4 for all three, so anything that moves was caused by Top N.")
doc <- tbl(doc, depth)
doc <- p(doc, "Every step deeper costs precision and buys recall. SapBERT goes from 0.724 precision at top 3 to 0.586 at top 30, while recall climbs 0.426 to 0.491. F1 peaks in the middle, at 10 for SapBERT and 20 for mpnet.")
doc <- p(doc, "ClinicalBERT falls the fastest, from 0.531 precision down to 0.299, so more than half of what a wider net adds is wrong. Its F1 is best at the narrowest setting tested. A wider net is affordable for SapBERT and it is not for ClinicalBERT.")

doc <- h(doc, "Where the correct answers come from")
doc <- p(doc, "Every code the pipeline emitted, split by which part of it put the code there, and whether the code was right.")
doc <- tbl(doc, src)
doc <- p(doc, "When both branches independently pick the same code it is right 96 to 99 percent of the time, on all three models. ")
doc <- p(doc, "The two branches are not equally good. On SapBERT the cosine branch alone is right 86 percent of the time and co-occurrence alone is right only 43 percent. On ClinicalBERT the cosine branch drops to 49 percent, while co-occurrence is slightly higher at 53%.")
doc <- p(doc, "The weakest source is the codes that are in both lists without either branch ranking them first, at 32 percent. ")

doc <- h(doc, "Where a correct code ranks")
doc <- tbl(doc, rank)
doc <- p(doc, "This is the cosine ranking on its own, before any cutoff, chapter filter or co-occurrence. Most ICD-9 codes have more than one correct ICD-10-CA code, so this counts the best placed of them. Rank 1 means a correct code is top of all 2038 candidates, not that every correct code is.")

doc <- h(doc, "How many correct codes each ICD-9 code has")
doc <- tbl(doc, targets)
doc <- p(doc, "937 validated pairs over 354 ICD-9 codes, median 2. 63 percent of ICD-9 codes have more than one correct ICD-10-CA code and one has 13.")
doc <- p(doc, "The zero row is the nine excluded codes after review. Four of them turned out to be validation errors and do have a correct code, 339 to G44, 175 to C50, 515 to J84 and 327 to G47. Five genuinely have none, 338, 249, 239, 445 and 209.")

print(doc, target = DOC)
cat("wrote", DOC, "\n")
