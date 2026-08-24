# compares the standard english stopword dictionaries on the icd label text
#
# the question is which standardized list to use for the clinicalbert re-run.
# every standard list treats single letters as stopwords, which breaks labels
# like "vitamin a deficiency" and "acute hepatitis a". this shows how much each
# list damages and which codes it merges.
#
# run from the repo root, takes a few seconds:
#   source("scripts/26_stopword_choice.R")

suppressMessages({library(readxl); library(dplyr); library(stopwords)})

LABELS <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"
VALID  <- "data/original/ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"

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
n_dup <- function(v) { tb <- table(v); sum(tb[tb > 1]) }

# read the labels once, everything below reuses this
cat("reading labels...\n")
DATA <- lapply(SHEETS, function(p) {
  d <- read_excel(LABELS, sheet = p[1])
  list(code = as.character(d[[p[2]]]), lab = base_clean(d[[p[3]]]))
})

# codes that are actually used as targets in the manual crosswalk. a collision
# only matters if the merged codes are ones we have to map to
tgt10 <- unique(as.character(read_excel(VALID, sheet = "Validation_ICD9_ICD10")$`ICD-10-CA`))
tgt8  <- unique(as.character(read_excel(VALID, sheet = "Validaion_ICD9_ICD8")$`ICDA-8`))
REAL <- list(`ICD-9-CM` = character(0), `ICD-10-CA` = tgt10, `ICDA-8` = tgt8)

LISTS <- list(
  snowball        = stopwords("en", source = "snowball"),
  nltk            = stopwords("en", source = "nltk"),
  smart           = stopwords("en", source = "smart"),
  `stopwords-iso` = stopwords("en", source = "stopwords-iso"))
# the proposed one: snowball, but keep single characters and negations
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

cat("\n\n============ 2. collisions, two codes becoming the same text ============\n")
summ <- list()
for (nm in names(LISTS)) {
  w <- LISTS[[nm]]
  tot <- 0; real <- 0; changed <- 0
  for (sn in names(DATA)) {
    d <- DATA[[sn]]
    sp <- strip_words(d$lab, w)
    tot <- tot + (n_dup(sp) - n_dup(d$lab))
    changed <- changed + sum(sp != d$lab)
    # collisions where at least one side is a code the crosswalk maps to
    tb <- tibble(code = d$code, lab = d$lab, sp = sp) %>%
      group_by(sp) %>% filter(n() > 1, n_distinct(lab) > 1) %>% ungroup()
    if (nrow(tb)) real <- real + sum(tb$code %in% REAL[[sn]])
  }
  summ[[nm]] <- tibble(list = nm, words = length(w), merged = tot,
                       merged_real = real, stripped = changed)
}
cat("merged      = pairs of codes that become identical text\n")
cat("merged_real = of those, ones the manual crosswalk actually maps to\n")
cat("stripped    = labels the list changes at all, i.e. how much work it does\n\n")
print(as.data.frame(bind_rows(summ)), row.names = FALSE)

cat("\n\n============ 3. what actually breaks ============\n")
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

cat("\n\n============ 4. recommendation ============\n")
cat("snowball is the safest standard list, but it still deletes \"a\", so\n")
cat("\"vitamin a deficiency\" and \"acute hepatitis a\" lose their meaning, and it\n")
cat("deletes \"not\", which merges the X00/X01 fire codes.\n\n")
cat("keeping single characters and no/not/nor fixes every collision and barely\n")
cat("reduces how much cleaning happens (see stripped above).\n")
