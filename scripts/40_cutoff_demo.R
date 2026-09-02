suppressMessages(source("scripts/pipeline_lib.R"))
suppressMessages(library(readxl))

with_code <- "data/original/Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"
no_code   <- "data/generated/cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx"
labels    <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"
threshold <- 0.995

lab10 <- read_excel(labels, sheet = "ICD-10-CA-3Level")
name_of <- setNames(as.character(lab10[[2]]), as.character(lab10[[1]]))

targets <- c("339" = "G44", "338" = "R52", "175" = "C50", "515" = "J84",
             "327" = "G47", "249" = "E13", "239" = "D48", "445" = "I74")

get_column <- function(path, icd9) {
  for (sh in excel_sheets(path)) {
    if (icd9 %in% names(read_excel(path, sheet = sh, n_max = 0))) {
      d <- read_excel(path, sheet = sh)
      return(data.frame(code = as.character(d[[1]]), score = as.numeric(d[[icd9]])))
    }
  }
  NULL
}

chapter_ok <- function(icd9, codes) {
  ch9 <- find_icd9cm_chapter(icd9)
  d <- vapply(codes, function(x)
        compute_chapter_distance(ch9, find_icd10ca_chapter(x), chapter_alignment_10), numeric(1))
  !is.na(d) & d < 1
}

emitted <- function(path, icd9) {
  x <- get_column(path, icd9)
  keep <- chapter_ok(icd9, x$code) & x$score >= threshold * max(x$score)
  if (!any(keep)) return(NA_character_)
  x$code[keep][which.max(x$score[keep])]
}

show_top <- function(path, icd9, tag) {
  x <- get_column(path, icd9)
  x$label <- name_of[x$code]
  top <- head(x[order(-x$score), ], 5)
  correct <- targets[icd9]
  cat("\nICD-9", icd9, "  correct target", correct, "  ", tag, "\n")
  for (i in 1:5) cat(sprintf("   %-4s %.4f  %s\n", top$code[i], top$score[i], top$label[i]))
  cat(sprintf("   %s scores %.4f and ranks %d\n", correct,
      x$score[x$code == correct], which(x$code[order(-x$score)] == correct)))
  e <- emitted(path, icd9)
  cat("   after the chapter filter only the highest scoring code is handed on, so the\n")
  cat("   similarity side gives", ifelse(is.na(e), "nothing", e), "\n")
}

show_top(with_code, "339", "code pasted in front of the label")
show_top(no_code,   "339", "label only")

cat("\n\nwhich code the similarity side hands over, for all 8\n\n")
r <- do.call(rbind, lapply(names(targets), function(k) {
  a <- emitted(with_code, k)
  b <- emitted(no_code, k)
  data.frame(icd9 = k, correct = unname(targets[k]),
             with_code = ifelse(is.na(a), "nothing", a),
             right = ifelse(!is.na(a) && a == targets[k], "yes", "no"),
             label_only = ifelse(is.na(b), "nothing", b),
             right_ = ifelse(!is.na(b) && b == targets[k], "yes", "no"))
}))
names(r)[c(4, 6)] <- c("right", "right ")
print(r, row.names = FALSE)
cat("\ncorrect with the code pasted in:", sum(r$right == "yes"), "of 8\n")
cat("correct with the label only:     ", sum(r$`right ` == "yes"), "of 8\n")
