# builds a standalone handout of the unmatched code descriptives, to send to
# someone who does not have the repo open. reads the csvs that
# 36_unmatched_descriptives.R writes, so run that first.
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

sim_tbl <- function(tk) {
  d <- sim %>% filter(model == MODEL, track == tk, group == "no reference match") %>%
    arrange(desc(max)) %>%
    transmute(icd9, icd9_label,
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
  @page { size: letter landscape; margin: 14mm 12mm; }
  body { font-family: Georgia, "Times New Roman", serif; font-size: 10.5pt; color: #111; line-height: 1.45; }
  h1 { font-size: 17pt; margin: 0 0 2px 0; }
  h2 { font-size: 12.5pt; margin: 22px 0 6px 0; border-bottom: 1px solid #999; padding-bottom: 3px; }
  h3 { font-size: 11pt; margin: 16px 0 4px 0; }
  .date { color: #555; font-size: 9.5pt; margin-bottom: 14px; }
  p { margin: 6px 0; max-width: 175mm; }
  table { border-collapse: collapse; width: 100%%; font-family: Helvetica, Arial, sans-serif;
          font-size: 7.6pt; margin: 8px 0 4px 0; }
  th { text-align: left; border-bottom: 1.2px solid #333; padding: 3px 5px; font-weight: 600; }
  td { border-bottom: 0.5px solid #ddd; padding: 2.5px 5px; vertical-align: top; }
  td.num { text-align: right; white-space: nowrap; }
  tbody tr:nth-child(even) { background: #f7f7f7; }
  img { width: 100%%; max-width: 250mm; margin: 8px 0; }
  .note { font-size: 9pt; color: #444; }
  .break { page-break-before: always; }
</style>

<h1>ICD-9-CM codes with no match in the target system</h1>
<div class="date">Descriptive analysis, 27 August 2026</div>

<p>There are 9 ICD-9-CM codes with no ICD-10-CA match and 52 with no ICDA-8 match.
The reference crosswalk lists them with no partner, so the right answer for them
is no mapping at all.</p>

<p>The algorithm uses two numbers to decide. How often the two codes appear
together in patient records, and the cosine similarity between the two code
labels. This describes both of those for each of these codes on its own, rather
than as an average over the crosswalk. The similarity scores are ClinicalBERT.</p>

<h2>ICD-9-CM to ICD-10-CA, the 9 codes</h2>

<h3>Co-occurrence</h3>

<p>All 9 codes appear in the co-occurrence data. Partners is the number of
different ICD-10-CA codes each one shares a record with, and the columns after it
describe the counts across those partners. The source data does not report a count below 7, which is why every code has a minimum of 7.</p>

%s

<p>These counts are not low. Code 338 shares a record with 827 different
ICD-10-CA codes, and its highest count is 13,431. Codes 515, 239 and 327 all
reach the thousands. For comparison, across the %d codes that do have a match the
middle value of the highest count is %s, so this group sits slightly above
that rather than below it.</p>

<p>The two exceptions are 445 and 339 at the bottom of the table. Those are the
only two that look the way a code with no valid target was expected to look.</p>

<p>What the counts seem to reflect is how often a code is recorded, not whether
it has a valid target. Pain, sleep disorders and unspecified neoplasms are
recorded constantly, so they turn up next to whatever else the patient has.</p>

<h3>Cosine similarity</h3>

<p>Each code is scored against all 2038 ICD-10-CA codes, so each row below
describes 2038 numbers.</p>

%s

<p>These scores are not low either. The highest score for these 9 runs from %s.
Across the %d codes that do have a match it runs from %s, so the two groups sit
on top of each other.</p>

<p>The closest ICD-10-CA code is clinically unrelated in 8 of the 9 rows. Code
175, malignant neoplasm of male breast, comes out closest to C30, malignant
neoplasm of nasal cavity and middle ear, at 0.955. Code 515, postinflammatory
pulmonary fibrosis, comes out closest to N72, inflammatory disease of cervix
uteri. The scores are as high as the ones the algorithm gets right, and the
answers are not close.</p>

<p>The one exception is code 338, whose closest target is R52, and both are
labelled pain not elsewhere classified. The reference crosswalk still records no
match for it, so that one may be worth checking rather than assuming the
algorithm is wrong.</p>

<div class="break"></div>
<h2>ICD-9-CM to ICDA-8, the 52 codes</h2>

<h3>Co-occurrence</h3>

<p>There is nothing to tabulate here. None of the 52 codes appear in the
co-occurrence data at all. A code is only listed if it shares a record with at
least one ICDA-8 code, so for these 52 there is no data rather than a low count.</p>

<p>This is why the algorithm returns no match for them, which is the answer we
want. It gets there by requiring co-occurrence data before it will accept a pair.
The same gap covers 32 of the %d codes that do have a valid ICDA-8 match, so it
stays silent on those as well.</p>

<h3>Cosine similarity</h3>

<p>Each code is scored against all 858 ICDA-8 codes.</p>

%s

<p>The highest score for these 52 runs from %s. Across the %d codes that do have
a match it runs from %s. As with ICD-10-CA, the two groups overlap.</p>

<div class="break"></div>
<h2>The same thing as a picture</h2>

<p>Each line is one code with no match. The grey line runs from its lowest to its
highest similarity score against any target code, the red bar covers the middle
half of its scores, the white dot is the median and the dark dot is the highest.</p>

<img src="plot_unmatched_similarity_spread.png">

<p>The next one puts the codes with no match against the codes that do have one.
Each point is a single ICD-9-CM code, placed at its highest similarity score.</p>

<img src="plot_unmatched_similarity_max.png">

<p class="note">The scores for the codes with no match sit inside the range of
the codes that do have one, on both crosswalks, so there is no cut-off on either
signal that would separate them without also removing a large share of the
mappings that are correct.</p>
',
  cooc_tbl("10_9"), n_of("10_9", "has a reference match"), cooc_mid("10_9"),
  sim_tbl("10_9"), rng("10_9", "no reference match"),
  n_of("10_9", "has a reference match"), rng("10_9", "has a reference match"),
  n_of("8_9", "has a reference match"),
  sim_tbl("8_9"), rng("8_9", "no reference match"),
  n_of("8_9", "has a reference match"), rng("8_9", "has a reference match"))

writeLines(html, file.path(OUT, "unmatched_codes_handout.html"))
cat("wrote results/unmatched_codes_handout.html\n")
