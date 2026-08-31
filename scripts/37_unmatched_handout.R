# builds a standalone handout of the unmatched code descriptives, to send to
# someone who does not have the repo open. reads the csvs that
# 36_unmatched_descriptives.R writes, so run that first.
#
# two or three pages. the 9 icd-10-ca codes are listed in full, the 52 icda-8
# codes are shown as the four highest and four lowest, since a 52 row table is
# not readable on a handout. the full tables are the csvs. no em dashes and no
# date in the heading, both asked for.
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
summ <- read.csv(file.path(OUT, "unmatched_descriptives_summary.csv"), stringsAsFactors = FALSE)

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

# middle value of the highest co-occurrence count, for the codes that do match,
# and how many of the unmatched codes beat it
cooc_mid_raw <- function(tk) {
  median(cooc$max[cooc$track == tk & cooc$group == "has a reference match" & !is.na(cooc$max)])
}
cooc_mid <- function(tk) formatC(cooc_mid_raw(tk), format = "f", digits = 1, big.mark = ",")
n_above_mid <- function(tk) {
  sum(cooc$track == tk & cooc$group == "no reference match" &
        !is.na(cooc$max) & cooc$max > cooc_mid_raw(tk))
}

# the cutoff that would clear every unmatched code, and what it costs
cut_of <- function(tk, field) {
  d <- summ %>% filter(model == MODEL, track == tk, signal == "cosine similarity",
                       group == "no reference match")
  if (field == "cut") f3(d$cutoff_clearing_all[1]) else d$cutoff_cost[1]
}

html <- sprintf('
<meta charset="utf-8">
<title>Codes with no match</title>
<style>
  @page { size: letter landscape; margin: 10mm 14mm; }
  body { font-family: Arial, Helvetica, sans-serif; font-size: 9.5pt; color: #000; line-height: 1.28; }
  h1 { font-size: 13pt; margin: 0 0 8px 0; }
  h2 { font-size: 11pt; margin: 10px 0 4px 0; }
  .break + h2 { margin-top: 0; }
  p { margin: 4px 0; max-width: 250mm; }
  table { border-collapse: collapse; width: 100%%; font-size: 7.3pt; line-height: 1.12; margin: 5px 0; }
  th { text-align: left; border-bottom: 1px solid #000; padding: 2px 5px 3px 5px; }
  td { border-bottom: 0.5px solid #ccc; padding: 1px 5px; vertical-align: top; }
  td.num { text-align: right; white-space: nowrap; }
  .break { page-break-before: always; }
</style>

<h1>ICD-9-CM codes with no match in the target system</h1>

<p>There are 9 ICD-9-CM codes with no ICD-10-CA match and 52 with no ICDA-8 match. The reference
crosswalk lists them with no partner, so for these codes the correct output is no mapping at all.</p>

<p>The algorithm decides using two numbers. One is how often a pair of codes shows up on the same
patient record. The other is the cosine similarity between the two code labels, which here comes
from ClinicalBERT. Both are below, for each code on its own.</p>

<h2>ICD-9-CM to ICD-10-CA, the 9 codes</h2>

<p>Co-occurrence first, since that is the one we expected to be low. Partners is how many different
ICD-10-CA codes each one turned up with, and the five columns after it describe the counts across
those partners. Nothing below 7 is reported in the source data, which is why every minimum is 7.</p>

%s

<p>They are not low. Code 338, pain not elsewhere classified, turned up with 827 different
ICD-10-CA codes, and with E11, type 2 diabetes mellitus, 13,431 times. Code 515 reaches 5,194 and
code 239 reaches 4,862. For comparison I took the %d ICD-9-CM codes that do have a match and found
the middle value of their highest count, which is %s. %d of these 9 are above it.</p>

<p>Two do behave the way we expected. Code 445, atheroembolism, has 2 partners and never gets past
16, and code 339, other headache syndromes, has 58 partners and a top count of 54.</p>

<p>The count seems to track how often a code gets written down, not whether it has a valid target.
Pain, sleep disorders and unspecified neoplasms get recorded on all sorts of patients, so they sit
next to whatever else the patient has.</p>

<p>Cosine similarity next. Each of the 9 is scored against all 2038 ICD-10-CA codes, so each row
below is a summary of 2038 numbers.</p>

%s

<p>The best score for these 9 runs from %s. For the %d codes that do have a match it runs from %s.
The two ranges sit almost on top of each other, so a high score on its own does not tell us whether
a code has a real match.</p>

<p>The closest code is also wrong in 8 of the 9, and not narrowly wrong. Code 175 is malignant
neoplasm of male breast and its closest ICD-10-CA code is C30, malignant neoplasm of nasal cavity
and middle ear, at 0.955. Code 515 is postinflammatory pulmonary fibrosis and its closest is N72,
inflammatory disease of cervix uteri. The model is picking up shared words rather than shared
meaning.</p>

<p>The one worth a second look is code 338. Its closest target is R52, and both are labelled pain
not elsewhere classified, but the crosswalk still records no match for it. I wanted to flag it in
case that one is an error in the reference rather than in the algorithm.</p>

<div class="break"></div>
<h2>ICD-9-CM to ICDA-8, the 52 codes</h2>

<p>These behave completely differently. None of the 52 appear in the co-occurrence data at all. Not
a low count, no row. A code is only listed there if it shares a record with at least one ICDA-8
code, and none of these ever did, so there is nothing to describe for that signal.</p>

<p>That works out in our favour, because it means the algorithm returns no match for all 52, which
is the right answer. But it gets there for a reason that costs us elsewhere. It will not accept a
pair without co-occurrence data, and 32 of the %d ICD-9-CM codes that do have a valid ICDA-8 match
are missing from that file too. It stays silent on those as well, and loses 39 mappings that are
correct.</p>

<p>Similarity looks much the same as it did on the other crosswalk. Each code is scored against all
858 ICDA-8 codes. The best score for the 52 runs from %s, against %s for the %d that do have a
match, so again there is no gap between the two groups. Rather than list all 52 I have put the four
highest and the four lowest below, and I have the rest if you want them.</p>

%s

<h2>Where this leaves the two signals</h2>

<p>Neither one separates these codes from the ones that do have a match. On ICD-10-CA, a similarity
cut-off high enough to clear all 9 would have to sit above %s, and that also removes %d of the %d
codes that are correct. On ICDA-8 the cut-off would be %s and it removes %d of %d. Co-occurrence
does not rescue it either, since the 9 sit on the high side of it and the 52 are not in the data at
all.</p>

<p>So I do not think a threshold on these two numbers is the way to handle these codes, which is
why I have not built one yet.</p>
',
  cooc_tbl("10_9"), n_of("10_9", "has a reference match"), cooc_mid("10_9"),
  n_above_mid("10_9"),
  sim_tbl("10_9"), rng("10_9", "no reference match"),
  n_of("10_9", "has a reference match"), rng("10_9", "has a reference match"),
  n_of("8_9", "has a reference match"),
  rng("8_9", "no reference match"), rng("8_9", "has a reference match"),
  n_of("8_9", "has a reference match"), sim_tbl("8_9", ends = 4),
  cut_of("10_9", "cut"), cut_of("10_9", "cost"), n_of("10_9", "has a reference match"),
  cut_of("8_9", "cut"), cut_of("8_9", "cost"), n_of("8_9", "has a reference match"))

writeLines(html, file.path(OUT, "unmatched_codes_handout.html"))
cat("wrote results/unmatched_codes_handout.html\n")
