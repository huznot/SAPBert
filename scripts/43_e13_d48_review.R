source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(source("scripts/pipeline_lib.R"))
suppressMessages({library(readxl); library(openxlsx)})

LAB  <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"
COOC <- "data/original/Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"
GEN  <- "data/generated"
OUT  <- out_path("e13_d48_review.xlsx")

arms <- list(
  list(model = "ClinicalBERT", file = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx", thr = 0.995, top_n = 25),
  list(model = "SapBERT",      file = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",      thr = 0.95,  top_n = 25),
  list(model = "mpnet",        file = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx",        thr = 0.95,  top_n = 30)
)

cases <- list(c(icd9 = "249", target = "E13"), c(icd9 = "239", target = "D48"))

t9  <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level")
t10 <- read_excel(LAB, sheet = "ICD-10-CA-3Level")
name9  <- setNames(as.character(t9[[2]]),  as.character(t9[[1]]))
name10 <- setNames(as.character(t10[[2]]), as.character(t10[[1]]))
cooc <- load_cooccurrence_df(COOC)

column_for <- function(path, icd9) {
  for (sh in excel_sheets(path)) {
    if (icd9 %in% names(read_excel(path, sheet = sh, n_max = 0))) {
      d <- read_excel(path, sheet = sh)
      return(data.frame(code = as.character(d[[1]]), score = as.numeric(d[[icd9]]),
                        stringsAsFactors = FALSE))
    }
  }
  NULL
}

rows <- list()
for (a in arms) {
  path <- file.path(GEN, a$file)
  co <- get_cooccurrence_codes_from_df(cooc, a$top_n, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
  for (cs in cases) {
    icd9 <- cs[["icd9"]]; target <- cs[["target"]]
    x <- column_for(path, icd9)
    x <- x[order(-x$score), ]
    x$rank <- seq_len(nrow(x))
    ch9 <- find_icd9cm_chapter(icd9)
    x$chapter <- vapply(x$code, function(cc) as.character(find_icd10ca_chapter(cc)), character(1))
    x$dist <- vapply(x$code, function(cc)
      compute_chapter_distance(ch9, find_icd10ca_chapter(cc), chapter_alignment_10), numeric(1))
    x$survives_chapter_filter <- ifelse(!is.na(x$dist) & x$dist < 1, "yes", "no")
    cut <- a$thr * max(x$score)
    x$clears_cutoff <- ifelse(x$score >= cut, "yes", "no")
    x$label <- ifelse(x$code %in% names(name10), name10[x$code], "")
    x$in_top_cooccurrence <- ifelse(
      x$code %in% as.character(co$ICD_10_CA[as.character(co$ICD_9_CM) == icd9]), "yes", "no")

    ts <- x$score[x$code == target]
    beaters <- x[x$score > ts, ]
    keep <- x[x$survives_chapter_filter == "yes" & x$clears_cutoff == "yes", ]

    cat(sprintf("\n%s   %s -> %s\n", a$model, icd9, target))
    cat(sprintf("%s scores %.4f, rank %d of %d\n", target, ts, x$rank[x$code == target], nrow(x)))
    if (nrow(keep)) {
      w <- keep[which.max(keep$score), ]
      cat(sprintf("pipeline gives %s  %.4f  %s\n", w$code, w$score, w$label))
    } else cat("pipeline gives nothing\n")
    b <- head(beaters, 5)
    cat(sprintf("beaten by %d codes, top %d below   (chapter ok / cutoff ok)\n", nrow(beaters), nrow(b)))
    for (i in seq_len(nrow(b)))
      cat(sprintf("  %-4s %.4f  %-3s %-3s  %s\n", b$code[i], b$score[i],
          b$survives_chapter_filter[i], b$clears_cutoff[i], b$label[i]))

    out <- x[, c("rank","code","label","score","chapter","dist",
                 "survives_chapter_filter","clears_cutoff","in_top_cooccurrence")]
    out$model <- a$model; out$icd9 <- icd9; out$target <- target
    out$cutoff <- cut; out$target_score <- ts
    out$beats_target <- ifelse(out$score > ts, "yes", "no")
    rows[[paste(a$model, icd9)]] <- out
  }
}
analysis <- c(

"Both of these codes are excluded in the validation data and that exclusion is correct. There is no three character ICD-10-CA code for either one. So E13 and D48 are not the right answers, and this sheet is not a list of things the pipeline missed. It is a record of what the models do with two codes that have nothing correct to return.",
"",
"",
"249  secondary diabetes mellitus",
"",
"ClinicalBERT ranks E13 191st out of 2038. The codes above it are Secondary parkinsonism, Diabetes mellitus in pregnancy, Ascorbic acid deficiency, Secondary hypertension, Bartonellosis, Hypertrichosis, Late syphilis. That column is not ordered by anything clinical. The model is going on the word secondary and on the general shape of a short ICD label. The whole top 191 spans only 0.036 of similarity, 0.9570 down to 0.9214, so nothing is separated from anything, and the cutoff at 0.9522 lets in 24 chapter surviving codes that have nothing to do with diabetes. E54 Ascorbic acid deficiency is handed over only because it is the highest scoring code that survives the chapter filter.",
"",
"SapBERT and mpnet behave much better. E13 sits at rank 3 and rank 4, and every code above it is another diabetes code, E11, E14, E10. Those columns are ordered sensibly. They still return a false positive, E11 and E14, but at least it is in the right part of the classification.",
"",
"That gap is the useful part. ClinicalBERT returns a vitamin deficiency for a diabetes code. SapBERT returns type 2 diabetes. Both are wrong against the validation data, but they are not equally wrong, and F1 counts them the same.",
"",
"",
"239  neoplasms of unspecified nature",
"",
"All three models rank D48 second. It is beaten by exactly one code every time, and every time that code is another neoplasm code. D15 on ClinicalBERT by 0.0021, C80 on SapBERT by 0.0526, D36 on mpnet by 0.0164.",
"",
"The competing labels are mostly shared boilerplate. Benign neoplasm of other and unspecified intrathoracic organs against Neoplasm of uncertain or unknown behaviour of other and unspecified sites. The only phrase that separates them is uncertain or unknown behaviour, a few words inside a long label that is otherwise near identical. A pooled sentence embedding averages that away.",
"",
"So the models are reading these labels the way you would expect. The problem is that there is no correct code to return and the pipeline has no way to say so.",
"",
"",
"What this says overall",
"",
"1. The pipeline cannot abstain. On both of these codes the right output is nothing, and nothing is not something it can produce. The similarity branch always returns its top scoring survivor, because the cutoff is a fraction of the column maximum and the top scorer always clears it. Co-occurrence always returns its top N. Every code these two produce is a false positive.",
"",
"2. False positives are not all the same size and the metrics cannot see that. ClinicalBERT returning Ascorbic acid deficiency for secondary diabetes and SapBERT returning Type 2 diabetes mellitus both count as one false positive.",
"",
"3. ClinicalBERT scores are badly compressed. 191 codes inside 0.036 on the 249 column, while SapBERT spans 0.85 down to 0.60 on the same column. That is a real argument for SapBERT and mpnet, and it also means a relative cutoff is not comparable across models.",
"",
"4. Neither E13 nor D48 is in the top co-occurrence codes for its ICD-9 code, so the two branches are not backing each other up here.",
"",
"5. All numbers come from the label only matrices in data/generated, that is with the ICD code no longer pasted in front of the label."
)

wb <- createWorkbook()
addWorksheet(wb, "analysis")
writeData(wb, "analysis", data.frame(review = analysis), colNames = FALSE)
setColWidths(wb, "analysis", 1, 150)
for (nm in names(rows)) {
  d <- rows[[nm]]
  sh <- gsub(" ", "_", nm)
  addWorksheet(wb, sh)
  writeData(wb, sh, d[d$beats_target == "yes" | d$code == d$target[1], ])
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, 1:ncol(d), c(6,7,60,9,7,7,10,9,10,14,7,8,9,12,7))
}
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
