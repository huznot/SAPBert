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

d <- d[, c("ICD-9-CM", "ICD-9-CM label", "Correct ICD-10-CA", "ICD-10-CA label",
           "Pipeline output", "Pipeline output label", "Verdict", "Notes")]

wb <- createWorkbook()
addWorksheet(wb, SH)
writeData(wb, SH, d, headerStyle = createStyle(textDecoration = "bold", valign = "bottom"))
addStyle(wb, SH, createStyle(valign = "top", wrapText = TRUE),
         rows = 2:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE, stack = TRUE)
setColWidths(wb, SH, cols = 1:8, widths = c(10, 28, 16, 32, 16, 32, 26, 55))
freezePane(wb, SH, firstActiveRow = 2)
saveWorkbook(wb, XLS, overwrite = TRUE)
cat("wrote", XLS, "\n")
