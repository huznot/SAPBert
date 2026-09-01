# formats the reviewed list of the 9 icd-9-cm codes with no icd-10-ca match as
# an xlsx to send out.
#
# it reads the xlsx it wrote last time, because the verdicts and notes are
# edited by hand in excel and those edits are the point. the csv is only the
# starting point, used when no xlsx exists yet. columns are selected by name at
# the end so re-running is safe and does not duplicate anything.
#
# the two label columns are looked up from ICD_Codes_Labels.xlsx rather than
# typed, so they carry the official ICD-10-CA wording. a cell holding more than
# one code gets one label per code in the same order.
#
# plain formatting on purpose: bold header, widths that fit the text, wrapped
# text, frozen top row, gridlines left on.
#
# usage: Rscript 39_unmatched_xlsx.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages({library(openxlsx); library(readxl)})

OUT <- "results"
SH  <- "Codes with no match"
XLS <- file.path(OUT, "unmatched_9_correct_targets.xlsx")
LAB <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"

if (file.exists(XLS)) {
  d <- read.xlsx(XLS)
  names(d) <- trimws(gsub(".", " ", names(d), fixed = TRUE))
} else {
  d <- read.csv(file.path(OUT, "unmatched_9_correct_targets.csv"), stringsAsFactors = FALSE)
  names(d) <- c("ICD-9-CM", "ICD-9-CM label", "Correct ICD-10-CA", "ICD-10-CA label",
                "in code list", "chapter aligned", "similarity", "sim rank",
                "Pipeline output", "found by", "Verdict", "Notes")
}

tgt <- read_excel(LAB, sheet = "ICD-10-CA-3Level")
lookup <- setNames(as.character(tgt[[2]]), as.character(tgt[[1]]))

CODE <- "[A-Z][0-9][0-9]"
codes_in <- function(cell) regmatches(cell, gregexpr(CODE, cell))[[1]]

# a cell can hold one code, several separated by slashes, or free text such as
# "none". label anything shaped like a 3 character icd-10 code; where a cell has
# no code, return NA so the hand written text can be kept instead
label_of <- function(x) {
  vapply(x, function(cell) {
    cs <- codes_in(cell)
    if (!length(cs)) return(NA_character_)
    paste(ifelse(cs %in% names(lookup), lookup[cs], "not in the ICD-10-CA code list"),
          collapse = " / ")
  }, character(1), USE.NAMES = FALSE)
}

new_correct <- label_of(d$`Correct ICD-10-CA`)
d$`ICD-10-CA label` <- ifelse(is.na(new_correct), d$`ICD-10-CA label`, new_correct)
d$`Pipeline output label` <- label_of(d$`Pipeline output`)

# the pipeline column had the label typed into some cells and not others, so
# strip it back to codes now that the label has a column of its own
d$`Pipeline output` <- vapply(d$`Pipeline output`, function(cell) {
  cs <- codes_in(cell)
  if (!length(cs)) cell else paste(cs, collapse = " / ")
}, character(1), USE.NAMES = FALSE)

# one plain sentence per code on why the pipeline landed where it did. the
# numbers come from the ClinicalBERT matrix at threshold 0.995, where the cutoff
# is 0.995 x the highest score in that code's column. kept in a named vector and
# re-applied every run, with any previous copy stripped off first, so re-running
# never stacks the same sentence twice
why <- c(
  "339" = "G44 was not the highest score in the column, an unrelated infection code A69 was, at 0.9155. that set the cutoff at 0.9109 and G44 at 0.9050 fell just under it. G43 then came in through co-occurrence, not similarity",
  "338" = "the chapter filter threw out R52 before the pipeline could pick it, so it had to choose from what was left",
  "175" = "this one worked, C50 scored 0.9532 and cleared the 0.9504 cutoff. C30 and C79 are extras that came in through co-occurrence",
  "249" = "the top score in the column was O24 diabetes in pregnancy at 0.9507, which set the cutoff at 0.9459. E13 at 0.9252 was cut, and E11 came in through co-occurrence",
  "239" = "the top score was D15 at 0.9416, setting the cutoff at 0.9369, and D48 at 0.9319 missed it by 0.005. C71 came in through co-occurrence",
  "445" = "the top score was M66 a tendon rupture code at 0.8527, setting the cutoff at 0.8485, so I74 at 0.8010 was cut on similarity and only came in through co-occurrence",
  "209" = "the pipeline has no way to answer no match. co-occurrence always hands back its top N codes, so something is always returned even when nothing fits"
)

MARK <- ". Why the pipeline missed it: "
d$Notes <- sub(paste0("[.]? ?", substring(MARK, 3), ".*$"), "", d$Notes)
w <- why[as.character(d$`ICD-9-CM`)]
d$Notes <- ifelse(is.na(w), d$Notes, paste0(trimws(d$Notes), MARK, w))

d <- d[, c("ICD-9-CM", "ICD-9-CM label", "Correct ICD-10-CA", "ICD-10-CA label",
           "Pipeline output", "Pipeline output label", "Verdict", "Notes")]

wb <- createWorkbook()
addWorksheet(wb, SH)
writeData(wb, SH, d, headerStyle = createStyle(textDecoration = "bold", valign = "bottom"))
addStyle(wb, SH, createStyle(valign = "top", wrapText = TRUE),
         rows = 2:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE, stack = TRUE)
setColWidths(wb, SH, cols = 1:8, widths = c(10, 28, 16, 32, 16, 32, 26, 80))
freezePane(wb, SH, firstActiveRow = 2)
saveWorkbook(wb, XLS, overwrite = TRUE)
cat("wrote", XLS, "\n")
