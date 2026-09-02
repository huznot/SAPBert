# what changes when stop words are removed and when the code number is taken
# out of the text, all with clinicalbert, everything else held fixed.
# reads the cached grids, nothing is recomputed, takes about a second.

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

R <- "results/grid/conditions"
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
# a cell that has not been generated yet should say so, not print NA
dlt <- function(a, b) if (is.na(a) || is.na(b)) "not run yet" else sprintf("%+.3f", a - b)

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
  cat(sprintf("    with the code in text   %s\n", dlt(c1, a)))
  cat(sprintf("    with the code removed   %s\n", dlt(d1, b)))
  cat("  effect of removing the code number:\n")
  cat(sprintf("    stop words kept         %s\n", dlt(b, a)))
  cat(sprintf("    stop words removed      %s\n", dlt(d1, c1)))
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

# bar chart of the four cells, one panel per track
suppressMessages({library(ggplot2)})
cells <- expand.grid(track = c("10_9", "8_9"),
                     stop = c("stop words kept", "stop words removed"),
                     code = c("code in text", "code removed"),
                     stringsAsFactors = FALSE)
tag_of <- function(stop, code) {
  base <- if (stop == "stop words kept") "clinicalbert_base" else "clinicalbert_stopwords"
  if (code == "code removed") paste0(base, "_nocode") else base
}
cells$f1 <- mapply(function(tr, s, cd) best(tag_of(s, cd), tr),
                   cells$track, cells$stop, cells$code)
cells <- cells[!is.na(cells$f1), ]
if (nrow(cells)) {
  cells$track_label <- tname(cells$track)
  g <- ggplot(cells, aes(code, f1, fill = stop)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", f1)),
              position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
    facet_wrap(~track_label, scales = "free_y") +
    scale_fill_manual(values = c("stop words kept" = "#6a9fd4",
                                 "stop words removed" = "#e08a4b")) +
    labs(x = NULL, y = "best F1", fill = NULL,
         title = "ClinicalBERT: stop words and code numbers",
         subtitle = "same parameter grid in every cell, only the input text changes") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(out_path("plot_stopwords_codes.png"), g, width = 8, height = 4.5, dpi = 150)
  cat("\nwrote results/plot_stopwords_codes.png\n")
}
