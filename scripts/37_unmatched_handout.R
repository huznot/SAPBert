# builds a standalone handout of the unmatched code descriptives, to send to
# someone who does not have the repo open. reads the csvs that
# 36_unmatched_descriptives.R writes, so run that first.
#
# one page. the 9 icd-10-ca codes are listed in full, the 52 icda-8 codes are
# summarised with the three highest and three lowest, since a 52 row table is
# not readable on a handout. the full tables are the csvs.
#
# writes results/unmatched_codes_handout.html. to get the pdf:
#   chrome --headless --disable-gpu --no-pdf-header-footer \
#     --print-to-pdf=results/unmatched_codes.pdf results/unmatched_codes_handout.html
#
# usage: Rscript 37_unmatched_handout.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages({library(dplyr)})

OUT <- "results"
MODEL <- "ClinicalBERT"

sim  <- read.csv(file.path(OUT, "unmatched_similarity_by_code.csv"), stringsAsFactors = FALSE)
cooc <- read.csv(file.path(OUT, "unmatched_cooccurrence_by_code.csv"), stringsAsFactors = FALSE)
pairs <- read.csv(file.path(OUT, "unmatched_top_cooccurrence_pairs.csv"), stringsAsFactors = FALSE)

f3 <- function(x) formatC(x, format = "f", digits = 3)
fn <- function(x) formatC(x, format = "d", big.mark = ",")
# quartiles land on halves when the count is even, so only show a decimal then
fq <- function(x) ifelse(x == round(x), formatC(x, format = "d", big.mark = ","),
                         formatC(x, format = "f", digits = 1, big.mark = ","))

# minimal html table, no library, the numbers are already rounded upstream
tbl <- function(df, headers, align_right) {
  th <- paste0("<th>", headers, "</th>", collapse = "")
  rows <- apply(df, 1, function(r) {
    td <- vapply(seq_along(r), function(i)
      sprintf("<td class=\"%s\">%s</td>", if (align_right[i]) "num" else "txt", r[i]), character(1))
    paste0("<tr>", paste(td, collapse = ""), "</tr>")
  })
  paste0("<table><thead><tr>", th, "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

# ends = NULL lists every code, ends = n keeps only the n highest and n lowest
sim_tbl <- function(tk, ends = NULL) {
  d <- sim %>% filter(model == MODEL, track == tk, group == "no reference match") %>%
    arrange(desc(max))
  if (!is.null(ends)) d <- bind_rows(head(d, ends), tail(d, ends))
  d <- d %>% transmute(icd9, icd9_label,
                       mean = f3(mean), min = f3(min), q1 = f3(q1), median = f3(median),
                       q3 = f3(q3), max = f3(max),
                       nearest = paste0(nearest_target, ", ", nearest_target_label))
  tbl(d, c("Code", "Label", "Average", "Min", "25th", "Median", "75th", "Max", "Closest target code"),
      c(FALSE, FALSE, rep(TRUE, 6), FALSE))
}

cooc_tbl <- function(tk) {
  top <- pairs %>% filter(track == tk) %>% group_by(icd9) %>%
    slice_max(freq, n = 1, with_ties = FALSE) %>% ungroup()
  d <- cooc %>% filter(track == tk, group == "no reference match", !is.na(partners)) %>%
    arrange(desc(max)) %>%
    left_join(top %>% transmute(icd9, tcode = target, tlab = target_label), by = "icd9") %>%
    transmute(icd9, icd9_label, partners = fn(partners),
              min = fn(min), q1 = fq(q1), median = fq(median), q3 = fq(q3),
              max = fn(max), top = paste0(tcode, ", ", tlab))
  tbl(d, c("Code", "Label", "Partners", "Min", "25th", "Median", "75th", "Max", "Most frequent partner"),
      c(FALSE, FALSE, rep(TRUE, 6), FALSE))
}

# group ranges, used in the sentences that say whether a score is high or low
rng <- function(tk, grp) {
  d <- sim %>% filter(model == MODEL, track == tk, group == grp)
  sprintf("%s to %s", f3(min(d$max)), f3(max(d$max)))
}
n_of <- function(track, grp) sum(sim$model == MODEL & sim$track == track & sim$group == grp)

# middle value of the highest co-occurrence count, for the codes that do match
cooc_mid <- function(tk) {
  d <- cooc %>% filter(track == tk, group == "has a reference match", !is.na(max))
  formatC(median(d$max), format = "f", digits = 1, big.mark = ",")
}

html <- sprintf('
<meta charset="utf-8">
<title>Codes with no match</title>
<style>
  @page { size: letter landscape; margin: 9mm 12mm; }
  body { font-family: Arial, Helvetica, sans-serif; font-size: 9pt; color: #000; line-height: 1.28; }
  h1 { font-size: 12.5pt; margin: 0 0 8px 0; }
  h2 { font-size: 10pt; margin: 9px 0 3px 0; }
  p { margin: 3px 0; max-width: 250mm; }
  table { border-collapse: collapse; width: 100%%; font-size: 7.2pt; line-height: 1.1; margin: 4px 0; }
  th { text-align: left; border-bottom: 1px solid #000; padding: 1px 5px 2px 5px; }
  td { border-bottom: 0.5px solid #ccc; padding: 1px 5px; vertical-align: top; }
  td.num { text-align: right; white-space: nowrap; }
</style>

<h1>ICD-9-CM codes with no match in the target system &mdash; 27 August 2026</h1>

<p>There are 9 ICD-9-CM codes with no ICD-10-CA match and 52 with no ICDA-8 match. The algorithm
uses two numbers to decide: how often a pair of codes appears together in patient records, and the
cosine similarity between the two code labels. Similarity here is from ClinicalBERT.</p>

<h2>ICD-9-CM to ICD-10-CA, the 9 codes</h2>

<p>Co-occurrence, all 9. Partners is the number of different ICD-10-CA codes the code shares a
record with, and the columns after it describe the counts across those partners. The source data
does not report a count below 7.</p>

%s

<p>These counts are not low. Code 338 shares a record with 827 different ICD-10-CA codes, up to
13,431 times. Across the %d codes that do have a match the middle value of the highest count is
%s, so this group sits above it rather than below. Only 445 and 339 are low.</p>

<p>Cosine similarity, each code scored against all 2038 ICD-10-CA codes.</p>

%s

<p>The highest score for these 9 runs from %s, against %s for the %d codes that do have a match.
The closest ICD-10-CA code is clinically unrelated in 8 of the 9: code 175, malignant neoplasm of
male breast, is closest to C30, malignant neoplasm of nasal cavity and middle ear. The exception is
code 338, whose closest target R52 carries the same label, yet the crosswalk records no match.</p>

<h2>ICD-9-CM to ICDA-8, the 52 codes</h2>

<p>None of the 52 appear in the co-occurrence data at all, so there is nothing to describe for that
signal. That is why the algorithm returns no match for them, and also why it returns no match for
39 codes that do have one, since 32 of the %d matched codes are missing from the same data.</p>

<p>Cosine similarity against all 858 ICDA-8 codes. The highest score for these 52 runs from %s,
against %s for the %d codes that do have a match. The three highest and the three lowest are below
as examples, and I have the other 46 if they are useful.</p>

%s

<p>On both crosswalks the codes with no match sit inside the range of the codes that do, so no
cut-off on either signal separates them without also losing a large share of the correct
mappings.</p>
',
  cooc_tbl("10_9"), n_of("10_9", "has a reference match"), cooc_mid("10_9"),
  sim_tbl("10_9"), rng("10_9", "no reference match"),
  rng("10_9", "has a reference match"), n_of("10_9", "has a reference match"),
  n_of("8_9", "has a reference match"),
  rng("8_9", "no reference match"), rng("8_9", "has a reference match"),
  n_of("8_9", "has a reference match"), sim_tbl("8_9", ends = 3))

writeLines(html, file.path(OUT, "unmatched_codes_handout.html"))
cat("wrote results/unmatched_codes_handout.html\n")
