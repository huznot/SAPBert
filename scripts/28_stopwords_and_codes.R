# what changes when stop words are removed and when the code number is taken
# out of the text, all with clinicalbert, everything else held fixed.
# reads the cached grids, nothing is recomputed, takes about a second.

R <- "results/full_grid"
get <- function(tag) {
  p <- file.path(R, paste0(tag, ".csv"))
  if (!file.exists(p)) return(NULL)
  read.csv(p, stringsAsFactors = FALSE)
}
best <- function(tag, tr) {
  d <- get(tag)
  if (is.null(d)) return(NA_real_)
  x <- d[d$track == tr, ]
  if (!nrow(x)) return(NA_real_)
  max(x$f1)
}
tname <- function(x) ifelse(x == "10_9", "ICD-9 to ICD-10-CA", "ICD-9 to ICDA-8")
fmt <- function(v) if (is.na(v)) "  --  " else sprintf("%.3f", v)

cat("\nall four cells use clinicalbert. only the text fed to it changes.\n")
cat("score is the best f1 over the same parameter grid in every cell.\n")

for (tr in c("10_9", "8_9")) {
  a <- best("clinicalbert_base", tr)
  b <- best("clinicalbert_base_nocode", tr)
  c1 <- best("clinicalbert_stopwords", tr)
  d1 <- best("clinicalbert_stopwords_nocode", tr)

  cat(sprintf("\n\n===== %s =====\n\n", tname(tr)))
  cat(sprintf("  %-22s %14s %14s\n", "", "code in text", "code removed"))
  cat(sprintf("  %-22s %14s %14s\n", "stop words kept", fmt(a), fmt(b)))
  cat(sprintf("  %-22s %14s %14s\n", "stop words removed", fmt(c1), fmt(d1)))

  cat("\n  effect of removing stop words:\n")
  cat(sprintf("    with the code in text   %+.3f\n", c1 - a))
  cat(sprintf("    with the code removed   %+.3f\n", d1 - b))
  cat("  effect of removing the code number:\n")
  cat(sprintf("    stop words kept         %+.3f\n", b - a))
  cat(sprintf("    stop words removed      %+.3f\n", d1 - c1))
}

cat("\n\n===== which stop word list =====\n\n")
cat("  two versions of the same standard list were run:\n")
cat("    snowball as published        175 words, deletes single letters a and i\n")
cat("    snowball keeping letters     170 words, keeps a, i, no, not, nor\n\n")
for (tr in c("10_9", "8_9")) {
  raw <- best("clinicalbert_stopwords_raw", tr)
  safe <- best("clinicalbert_stopwords", tr)
  cat(sprintf("  %-22s as published %s   keeping letters %s\n",
              tname(tr), fmt(raw), fmt(safe)))
}
cat("\n  keeping the letters matters because 'vitamin a deficiency' and\n")
cat("  'acute hepatitis a' lose the letter that identifies them otherwise.\n")
cat("  see 26_stopword_choice.R for the labels this affects.\n")
