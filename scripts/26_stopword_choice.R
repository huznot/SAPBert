suppressMessages({library(readxl); library(dplyr); library(stopwords)})

LABELS <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"

SHEETS <- list(
  `ICD-9-CM`  = c("CCS ICD-9-CM-3Level", "ICD_9_CM",  "ICD_9_CM_LABEL"),
  `ICD-10-CA` = c("ICD-10-CA-3Level",    "ICD_10_CA", "ICD_10_CA_LABEL"),
  `ICDA-8`    = c("ICDA-8-3Level",       "ICDA_8",    "ICDA_8_LABEL"))

base_clean <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9 ]", " ", x)
  trimws(gsub("[ ]+", " ", x))
}
strip_words <- function(x, w) {
  vapply(strsplit(x, " "),
         function(t) paste(t[!(t %in% w)], collapse = " "), character(1))
}

# read the labels once, everything below reuses this
cat("reading labels...\n")
DATA <- lapply(SHEETS, function(p) {
  d <- read_excel(LABELS, sheet = p[1])
  list(code = as.character(d[[p[2]]]), lab = base_clean(d[[p[3]]]))
})

LISTS <- list(
  snowball        = stopwords("en", source = "snowball"),
  nltk            = stopwords("en", source = "nltk"),
  smart           = stopwords("en", source = "smart"),
  `stopwords-iso` = stopwords("en", source = "stopwords-iso"))
# snowball, but keeping single characters and negations
LISTS[["snowball, keep letters + not"]] <-
  setdiff(LISTS$snowball, c(letters, "no", "not", "nor"))

cat("\n============ 1. what each list would delete ============\n")
for (nm in names(LISTS)) {
  w <- LISTS[[nm]]
  cat(sprintf("\n%-30s %4d words\n", nm, length(w)))
  cat(sprintf("   deletes single letters : %s\n",
              ifelse(any(nchar(w) == 1), paste(w[nchar(w) == 1], collapse = " "),
                     "NONE, letters kept")))
  cat(sprintf("   deletes negations      : %s\n",
              ifelse(any(c("no","not","nor") %in% w),
                     paste(intersect(c("no","not","nor"), w), collapse = " "),
                     "NONE, negations kept")))
}

cat("\n\n============ 2. what breaks ============\n")
show_damage <- function(nm, codes_of_interest) {
  w <- LISTS[[nm]]
  cat(sprintf("\n-- %s --\n", nm))
  for (sn in names(DATA)) {
    d <- DATA[[sn]]
    idx <- which(d$code %in% codes_of_interest)
    if (!length(idx)) next
    sp <- strip_words(d$lab[idx], w)
    for (i in seq_along(idx))
      cat(sprintf("   %-9s %-5s %-38s -> %s\n", sn, d$code[idx[i]],
                  substr(d$lab[idx[i]], 1, 38), sp[i]))
  }
}
INTEREST <- c("B15", "E50", "E55", "260", "265", "268", "C84", "X00", "X01")
for (nm in c("nltk", "snowball", "snowball, keep letters + not")) show_damage(nm, INTEREST)
