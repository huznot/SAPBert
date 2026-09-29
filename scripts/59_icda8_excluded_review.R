source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(openxlsx))

# the review of the 52 icd-9 codes with no icda-8 match, laid out the way the
# nine codes on the icd-9 to icd-10-ca track were done: the ones that look
# wrong first, then every code described, then the chapters, then the signals.
#
# her instruction was to assume the validation data is right unless something
# obvious turns up, so only six are marked wrong and the reasoning travels with
# each one. the accuracy effect is the last tab, it is the smallest part.
#
# reads: 52_codes_no_icda8.xlsx      (52_icda8_excluded_list.R, the verdicts)
#        unmatched_*_by_code.csv     (36_unmatched_descriptives.R, the signals)
#        52_codes_accuracy_effect.csv (58_icda8_excluded_effect.R)

ORIG <- "data/original"
GEN  <- "data/generated"
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("52_codes_review.xlsx")

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_8_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_8_9_sapbert_base_nocode.xlsx")

base <- read_excel(out_path("52_codes_no_icda8.xlsx"), sheet = "no ICDA-8 match")
base$`ICD-9-CM` <- as.character(base$`ICD-9-CM`)
stopifnot(nrow(base) == 52)

a8    <- read_excel(LAB, sheet = "ICDA-8-3Level")
a8lab <- setNames(as.character(a8$ICDA_8_LABEL), as.character(a8$ICDA_8))

# "163 other and unspecified respiratory organs" -> "163"
base$proposed <- ifelse(base$`Closest ICDA-8 code` == "none", NA_character_,
                        sub(" .*$", "", base$`Closest ICDA-8 code`))

## what each model actually scores for these codes. for the six with a proposed
## match this says whether the label similarity could ever have found it, which
## is the "why is the correct code missing" column from the nine code review
top_by_model <- list(); rank_by_model <- list()
for (m in names(MODELS)) {
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[m]]))
  cols <- lapply(sheets, function(df) names(df))
  tops <- character(nrow(base)); rks <- character(nrow(base))
  for (i in seq_len(nrow(base))) {
    k <- base$`ICD-9-CM`[i]
    sh <- which(vapply(cols, function(cn) k %in% cn, logical(1)))
    if (!length(sh)) { tops[i] <- "not in the matrix"; rks[i] <- ""; next }
    df <- sheets[[sh[1]]]
    d <- data.frame(t = as.character(df[[names(df)[1]]]), s = df[[k]], stringsAsFactors = FALSE)
    d <- d[!is.na(d$s), ]; d <- d[order(-d$s), ]
    tops[i] <- sprintf("%s %s (%.3f)", d$t[1], a8lab[d$t[1]], d$s[1])
    p <- base$proposed[i]
    rks[i] <- if (is.na(p)) "" else {
      r <- match(p, d$t)
      if (is.na(r)) "not scored" else sprintf("#%d of %d (%.3f)", r, nrow(d), d$s[r])
    }
  }
  top_by_model[[m]] <- tops; rank_by_model[[m]] <- rks
  cat("read", m, "\n")
}

base$`ClinicalBERT top cosine`   <- top_by_model$ClinicalBERT
base$`SapBERT top cosine`        <- top_by_model$SapBERT
base$`Rank of proposed, ClinicalBERT` <- rank_by_model$ClinicalBERT
base$`Rank of proposed, SapBERT`      <- rank_by_model$SapBERT
base$`In the co-occurrence data` <- ifelse(base$`Co-occurrence rows` > 0, "yes", "no")

wrong <- base[base$`Validation correct?` == "validation data wrong", ]
right <- base[base$`Validation correct?` != "validation data wrong", ]
stopifnot(nrow(wrong) == 6, nrow(right) == 46)

## tab 1. the six that look wrong
six <- data.frame(
  `ICD-9-CM` = wrong$`ICD-9-CM`,
  `ICD-9-CM label` = wrong$`ICD-9-CM label`,
  `ICD-9 chapter` = wrong$`ICD-9 chapter`,
  `Proposed ICDA-8` = wrong$proposed,
  `ICDA-8 label` = unname(a8lab[wrong$proposed]),
  `Why it looks wrong` = wrong$Notes,
  `ClinicalBERT top cosine` = wrong$`ClinicalBERT top cosine`,
  `Rank of proposed, ClinicalBERT` = wrong$`Rank of proposed, ClinicalBERT`,
  `SapBERT top cosine` = wrong$`SapBERT top cosine`,
  `Rank of proposed, SapBERT` = wrong$`Rank of proposed, SapBERT`,
  `In the co-occurrence data` = wrong$`In the co-occurrence data`,
  check.names = FALSE)
six <- six[order(as.numeric(six$`ICD-9-CM`)), ]

## tab 2. the other 46, left as the validation data has them
rest <- data.frame(
  `ICD-9-CM` = right$`ICD-9-CM`,
  `ICD-9-CM label` = right$`ICD-9-CM label`,
  `ICD-9 chapter` = right$`ICD-9 chapter`,
  `Why there is no ICDA-8 rubric` = right$Notes,
  `ClinicalBERT top cosine` = right$`ClinicalBERT top cosine`,
  `SapBERT top cosine` = right$`SapBERT top cosine`,
  `In the co-occurrence data` = right$`In the co-occurrence data`,
  `Top co-occurring ICDA-8` = right$`Top co-occurring ICDA-8`,
  check.names = FALSE)
rest <- rest[order(as.numeric(rest$`ICD-9-CM`)), ]

## tab 3. which chapters the 52 come from
ch <- as.data.frame(table(base$`ICD-9 chapter`), stringsAsFactors = FALSE)
names(ch) <- c("ICD-9 chapter", "Codes with no ICDA-8 match")
ch$`Of which look wrong` <- vapply(ch$`ICD-9 chapter`,
  function(x) sum(wrong$`ICD-9 chapter` == x), integer(1))
ch$Share <- sprintf("%.0f%%", 100 * ch$`Codes with no ICDA-8 match` / nrow(base))
ch <- ch[order(-ch$`Codes with no ICDA-8 match`), ]

## tab 4. the two signals per code, with the codes that do have a match as the
## comparison group. straight out of 36_unmatched_descriptives.R
sim <- read.csv(out_path("unmatched_similarity_by_code.csv"), stringsAsFactors = FALSE)
sim <- sim[sim$track == "8_9" & sim$model %in% names(MODELS), ]
sig <- data.frame(
  Model = sim$model, Group = sim$group,
  `ICD-9-CM` = as.character(sim$icd9), `ICD-9-CM label` = sim$icd9_label,
  `CCS category` = sim$ccs_category,
  `ICDA-8 codes scored` = sim$n,
  Mean = sim$mean, Lowest = sim$min, `25th` = sim$q1, Median = sim$median,
  `75th` = sim$q3, Highest = sim$max,
  `Nearest ICDA-8` = sim$nearest_target, `Nearest label` = sim$nearest_target_label,
  check.names = FALSE)
sig <- sig[order(sig$Model, sig$Group != "no reference match",
                 as.numeric(sig$`ICD-9-CM`)), ]

## tab 5. the accuracy effect, last and smallest
eff <- read.csv(out_path("52_codes_accuracy_effect.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)
names(eff) <- gsub("\\.", " ", names(eff))

cat(sprintf("\n%d codes reviewed, %d look wrong, %d left as the validation data has them\n",
            nrow(base), nrow(wrong), nrow(right)))
cat("in the co-occurrence data:", sum(base$`In the co-occurrence data` == "yes"), "of 52\n\n")
print(six[, c("ICD-9-CM", "Proposed ICDA-8", "Rank of proposed, SapBERT")], row.names = FALSE)

hdr <- createStyle(textDecoration = "bold", valign = "bottom", fgFill = "#D9D9D9",
                   border = "bottom", borderStyle = "medium", wrapText = TRUE)
body <- createStyle(valign = "top", wrapText = TRUE)
CB <- createStyle(fgFill = "#EAF1FB"); SB <- createStyle(fgFill = "#FDF1E4")
wb <- createWorkbook()
add <- function(name, df, widths, filter = FALSE, freeze_col = 1, wrap = TRUE) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr, withFilter = filter)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = freeze_col)
  if (wrap) addStyle(wb, name, body, rows = 2:(nrow(df) + 1), cols = seq_along(widths),
                     gridExpand = TRUE, stack = TRUE)
  if ("Model" %in% names(df))
    for (p in list(list("ClinicalBERT", CB), list("SapBERT", SB))) {
      r <- which(df$Model == p[[1]]) + 1
      if (length(r)) addStyle(wb, name, p[[2]], rows = r, cols = seq_along(widths),
                              gridExpand = TRUE, stack = TRUE)
    }
}
add("the six that look wrong", six, c(10, 34, 28, 11, 38, 60, 34, 17, 34, 17, 12), freeze_col = 2)
add("the other 46", rest, c(10, 34, 28, 52, 34, 34, 12, 40), filter = TRUE, freeze_col = 2)
add("by chapter", ch, c(42, 16, 13, 8))
add("similarity and co-occurrence", sig,
    c(14, 22, 10, 34, 28, 12, 9, 9, 9, 9, 9, 9, 13, 34), filter = TRUE, freeze_col = 4)
add("accuracy effect", eff, c(14, 8, 30, 12, 12, 16, 14, 15, 15, 14, 16),
    freeze_col = 2, wrap = FALSE)
saveWorkbook(wb, OUT, overwrite = TRUE)
cat("\nwrote", OUT, "\n")
