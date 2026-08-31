# builds a standalone handout of the unmatched code descriptives, to send to
# someone who does not have the repo open. reads the csvs written by
# 36_unmatched_descriptives.R and 35_unmatched_codes.R, so run those first.
#
# built against the 27 aug meeting, and it covers, in her order:
#   co-occurrence for the 9, are any of them really high
#   cosine similarity against every target code, average min max 25th 75th
#   the same for the codes that do have a match, as the comparison group
#   the groups described together, and then each code on its own
#   what the algorithm actually emits for them, which is the false positives
#   whether a threshold on either signal would catch them
#   what the two thresholds are really doing, pipeline_lib.R lines 153 to 190
#   where the work lives in the repo, which she asked for by name
#
# no em dashes and no date in the heading, both asked for. the html is one
# single quoted R string, so apostrophes in the prose need a backslash.
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

sim   <- read.csv(file.path(OUT, "unmatched_similarity_by_code.csv"), stringsAsFactors = FALSE)
cooc  <- read.csv(file.path(OUT, "unmatched_cooccurrence_by_code.csv"), stringsAsFactors = FALSE)
pairs <- read.csv(file.path(OUT, "unmatched_top_cooccurrence_pairs.csv"), stringsAsFactors = FALSE)
summ  <- read.csv(file.path(OUT, "unmatched_descriptives_summary.csv"), stringsAsFactors = FALSE)
hand  <- read.csv(file.path(OUT, "unmatched_code_handling.csv"), stringsAsFactors = FALSE)

f3 <- function(x) ifelse(is.na(x), "no data", formatC(x, format = "f", digits = 3))
fn <- function(x) ifelse(is.na(x), "no data", formatC(x, format = "d", big.mark = ","))
# quartiles land on halves when the count is even, so only show a decimal then
fq <- function(x) ifelse(is.na(x), "no data",
                  ifelse(x == round(x), formatC(x, format = "d", big.mark = ","),
                         formatC(x, format = "f", digits = 1, big.mark = ",")))

# minimal html table, no library, the numbers are already rounded upstream
tbl <- function(df, headers, align_right, cls = "") {
  th <- paste0("<th>", headers, "</th>", collapse = "")
  rows <- apply(df, 1, function(r) {
    td <- vapply(seq_along(r), function(i)
      sprintf("<td class=\"%s\">%s</td>", if (align_right[i]) "num" else "txt", r[i]), character(1))
    paste0("<tr>", paste(td, collapse = ""), "</tr>")
  })
  paste0("<table class=\"", cls, "\"><thead><tr>", th, "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

# ---- the two groups described together, which is the first thing she asked
# for, and it needs the codes that do have a match as the comparison ----
GRP <- c("no reference match", "has a reference match")
lab <- c("no reference match" = "no match", "has a reference match" = "has a match")
tgt <- c("10_9" = "ICD-10-CA", "8_9" = "ICDA-8")

grp_rows <- function(sig) {
  summ %>% filter(signal == sig, is.na(model) | model == MODEL) %>%
    mutate(o = match(group, GRP)) %>% arrange(track, o)
}

grp_sim_tbl <- function() {
  d <- grp_rows("cosine similarity") %>%
    transmute(target = tgt[track], group = lab[group], codes = fn(codes),
              mean = f3(mean_similarity), min = f3(min_of_max), q1 = f3(q1_of_max),
              median = f3(median_of_max), q3 = f3(q3_of_max), max = f3(max_of_max))
  tbl(d, c("Target system", "Group", "Codes", "Average of every score",
           "Best score, lowest", "25th", "Median", "75th", "Best score, highest"),
      c(FALSE, FALSE, rep(TRUE, 7)), "grp")
}

grp_cooc_tbl <- function() {
  d <- grp_rows("co-occurrence frequency") %>%
    transmute(target = tgt[track], group = lab[group], codes = fn(codes),
              absent = fn(no_cooccurrence_data), partners = fq(median_partners),
              min = fn(min_of_max), q1 = fq(q1_of_max), median = fq(median_of_max),
              q3 = fq(q3_of_max), max = fn(max_of_max))
  tbl(d, c("Target system", "Group", "Codes", "Not in the file", "Median partners",
           "Highest count, lowest", "25th", "Median", "75th", "Highest count, highest"),
      c(FALSE, FALSE, rep(TRUE, 8)), "grp")
}

# ---- per code ----
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

# ---- numbers used in the prose ----
rng <- function(tk, grp) {
  d <- sim %>% filter(model == MODEL, track == tk, group == grp)
  sprintf("%s to %s", f3(min(d$max)), f3(max(d$max)))
}
n_of <- function(track, grp) sum(sim$model == MODEL & sim$track == track & sim$group == grp)

# the weakest best score in the whole crosswalk. it still produces a candidate,
# which is what shows the similarity rule has no floor
weakest_best <- function(tk) f3(min(sim$max[sim$model == MODEL & sim$track == tk]))
floor_at <- function(tk, thr) f3(thr * min(sim$max[sim$model == MODEL & sim$track == tk]))

cooc_mid_raw <- function(tk) {
  median(cooc$max[cooc$track == tk & cooc$group == "has a reference match" & !is.na(cooc$max)])
}
cooc_mid <- function(tk) formatC(cooc_mid_raw(tk), format = "f", digits = 1, big.mark = ",")
n_above_mid <- function(tk) {
  sum(cooc$track == tk & cooc$group == "no reference match" &
        !is.na(cooc$max) & cooc$max > cooc_mid_raw(tk))
}
# matched codes actually in the co-occurrence file, the denominator the
# co-occurrence cutoff cost is counted against
n_cooc_present <- function(tk) {
  sum(cooc$track == tk & cooc$group == "has a reference match" & !is.na(cooc$max))
}
cooc_span <- function(tk, which) {
  v <- cooc$max[cooc$track == tk & !is.na(cooc$max)]
  fn(if (which == "low") min(v) else max(v))
}

# what the pipeline actually emits for these codes, from 35_
h <- function(tk, col) hand[[col]][hand$track == tk & hand$model == MODEL]

# the cutoff that would clear every unmatched code, and what it costs
cut_of <- function(tk, sig, field) {
  d <- summ %>% filter(track == tk, signal == sig, group == "no reference match",
                       is.na(model) | model == MODEL)
  if (field == "cut") {
    if (sig == "cosine similarity") f3(d$cutoff_clearing_all[1]) else fn(d$cutoff_clearing_all[1])
  } else d$cutoff_cost[1]
}

style <- '
<meta charset="utf-8">
<title>Codes with no match</title>
<style>
  @page { size: letter portrait; margin: 25.4mm; }
  body { font-family: Arial, Helvetica, sans-serif; font-size: 10pt; color: #000; line-height: 1.35; }
  h1 { font-size: 14pt; margin: 0 0 12px 0; }
  h2 { font-size: 11.5pt; margin: 18px 0 6px 0; }
  p { margin: 7px 0; }
  table { border-collapse: collapse; width: 100%; font-size: 7pt; line-height: 1.15; margin: 8px 0; }
  th { text-align: left; border-bottom: 1px solid #000; padding: 2px 5px 3px 5px; }
  td { border-bottom: 0.5px solid #ccc; padding: 1px 5px; vertical-align: top; }
  td.num { text-align: right; white-space: nowrap; }
  h2 { page-break-after: avoid; }
  table, img { page-break-inside: avoid; }
  table.grp td.txt { white-space: nowrap; }
  img { max-width: 100%; height: auto; margin: 8px 0; }
</style>
'

intro <- '
<h1>ICD-9-CM codes with no match in the target system</h1>

<p>There are 9 ICD-9-CM codes with no ICD-10-CA match and 52 with no ICDA-8 match. The reference
crosswalk lists them with no partner, so for these codes the correct output is no mapping at all,
and any mapping the algorithm produces for them is a false positive.</p>

<p>The algorithm decides using two numbers. One is how often a pair of codes shows up on the same
patient record. The other is the cosine similarity between the two code labels, which here comes
from ClinicalBERT. Below are both numbers, first for the two groups together and then for each of
these codes on its own, with the codes that do have a match alongside for comparison. The
analysis is scripts/36_unmatched_descriptives.R and scripts/35_unmatched_codes.R on the main
branch, written up in report.md Section 11.</p>
'

outcome <- sprintf('
<h2>Algorithm output</h2>

<p>The outcome that the rest of this describes. On ICD-9-CM to
ICD-10-CA the algorithm produces a mapping for all %d of the codes that should get none. Every one
of them is a false positive and not one is correctly rejected. On ICD-9-CM to ICDA-8 it is the
other way round. It stays silent on all %d, which is the right answer, but it also stays silent on
%d of the %d codes that do have a valid match, so it is getting the special cases right and paying
for it on the ordinary ones.</p>
',
  h("10_9", "wrongly_mapped"), h("8_9", "correctly_silent"),
  h("8_9", "wrongly_silent"), h("8_9", "n_matched"))

groups <- sprintf('
<h2>Group summary</h2>

<p>Co-occurrence first, this being the signal expected to be low. A code is only listed in the
co-occurrence file if it shares a record with at least one target code, so a code that is not in
the file has no data rather than a low count. Partners is how many different target codes a code
turned up with, and the last five columns describe each group by the highest count any one of its
codes reached. Pairs below a frequency of 6 were removed when the
co-occurrence data was built, and the lowest count present is 7.</p>

%s

<p>The ICD-10-CA group with no match is not low. Its middle code reaches a higher top count than
the middle code of the group that does have a match. The ICDA-8 group with no match is not low
either, it is simply not there at all.</p>

<p>Then cosine similarity. Each ICD-9-CM code is scored against every target code, 2038 of them for
ICD-10-CA and 858 for ICDA-8, so the average column below is over all of those scores and all of
the codes in the group. The last five columns describe the group by each code\'s single best score,
which is the one that decides whether a mapping is made.</p>

%s

<p>On both crosswalks the group with no match sits inside the range of the group that has one. That
is the finding in one line, and the per code tables below are where it comes from.</p>
',
  grp_cooc_tbl(), grp_sim_tbl())

t10 <- sprintf('
<h2>ICD-9-CM to ICD-10-CA</h2>

<p>Co-occurrence.</p>

%s

<p>The counts are not low. Code 338, pain not elsewhere classified, turned up with 827 different
ICD-10-CA codes, and with E11, type 2 diabetes mellitus, 13,431 times. Code 515 reaches 5,194 and
code 239 reaches 4,862. Across the %d ICD-9-CM codes that do have a match the middle value of the
highest count is %s, and %d of these 9 are above it.</p>

<p>Two behave as expected. Code 445, atheroembolism, has 2 partners and never gets past
16, and code 339, other headache syndromes, has 58 partners and a top count of 54.</p>

<p>The count seems to track how often a code gets written down, not whether it has a valid target.
Pain, sleep disorders and unspecified neoplasms get recorded on all sorts of patients, so they sit
next to whatever else the patient has.</p>

<p>Cosine similarity. Each row below is a summary of that code\'s 2038 scores.</p>

%s

<p>The best score for these 9 runs from %s. For the %d codes that do have a match it runs from %s.
The two ranges sit almost on top of each other, so a high score on its own does not tell us whether
a code has a real match.</p>

<p>The closest code is also wrong in 8 of the 9, and not narrowly wrong. Code 175 is malignant
neoplasm of male breast and its closest ICD-10-CA code is C30, malignant neoplasm of nasal cavity
and middle ear, at 0.955. Code 515 is postinflammatory pulmonary fibrosis and its closest is N72,
inflammatory disease of cervix uteri. The model is picking up shared words rather than shared
meaning, and with both signals pointing the same way all 9 come out as mappings.</p>

<p>The one worth a second look is code 338. Its closest target is R52, and both are labelled pain
not elsewhere classified, but the crosswalk still records no match for it. That one may be an
error in the reference rather than in the algorithm.</p>
',
  cooc_tbl("10_9"), n_of("10_9", "has a reference match"), cooc_mid("10_9"),
  n_above_mid("10_9"), sim_tbl("10_9"), rng("10_9", "no reference match"),
  n_of("10_9", "has a reference match"), rng("10_9", "has a reference match"))

t8 <- sprintf('
<h2>ICD-9-CM to ICDA-8</h2>

<p>These behave completely differently. None of the 52 appear in the co-occurrence data at all. Not
a low count, no row, so there is nothing to tabulate for that signal.</p>

<p>This produces the right answer, since the algorithm returns no match for all 52, but it gets
there for a reason that costs the crosswalk elsewhere. It will not accept a pair without co-occurrence
data, and 32 of the %d ICD-9-CM codes that do have a valid ICDA-8 match are missing from that file
too, which is where the %d lost mappings come from. The 52 are being got right by an accident of
missing data rather than by either threshold doing its job.</p>

<p>Similarity looks much the same as it did on the other crosswalk. The best score for the 52 runs
from %s, against %s for the %d that do have a match, so again there is no gap between the two
groups. The four highest and the four lowest are below, and the full 52 are in the summary
files.</p>

%s
',
  n_of("8_9", "has a reference match"), h("8_9", "wrongly_silent"),
  rng("8_9", "no reference match"), rng("8_9", "has a reference match"),
  n_of("8_9", "has a reference match"), sim_tbl("8_9", ends = 4))

figs <- '
<h2>Distributions</h2>

<p>Each point is one ICD-9-CM code, placed at its highest cosine similarity against any target
code. In both panels the upper row is the codes that do have a match and the lower row is the codes
that do not, and the box covers the middle half of each group.</p>

<img src="plot_unmatched_similarity_max.png">

<p>The same layout for co-occurrence, each code placed at its highest count against any target
code, on a log scale. The lower row of the ICDA-8 panel is empty because none of those 52 codes are
in the file.</p>

<img src="plot_unmatched_cooccurrence.png">
'

thresh <- sprintf('
<h2>Threshold scenarios</h2>

<p>Not on its own. On ICD-10-CA a similarity cut-off high enough to clear all 9 would have to sit
above %s, and that also removes %d of the %d codes that are correct. On ICDA-8 the cut-off would be
%s and it removes %d of %d. Co-occurrence is no better. A count cut-off high enough to clear the 9
would have to sit above %s, which removes %d of the %d matched codes that are in the file. The two
together do not help either, because on ICD-10-CA the 9 pass both tests comfortably, and on ICDA-8
the 52 are already being caught by the missing data rather than by any threshold.</p>
',
  cut_of("10_9", "cosine similarity", "cut"), cut_of("10_9", "cosine similarity", "cost"),
  n_of("10_9", "has a reference match"),
  cut_of("8_9", "cosine similarity", "cut"), cut_of("8_9", "cosine similarity", "cost"),
  n_of("8_9", "has a reference match"),
  cut_of("10_9", "co-occurrence frequency", "cut"),
  cut_of("10_9", "co-occurrence frequency", "cost"), n_cooc_present("10_9"))

mech <- sprintf('
<h2>Threshold implementation</h2>

<p>Step 1 of the Methods defines the similarity cutoff as the maximum similarity score for a code
multiplied by the threshold value, which makes it relative to the code being mapped rather than a
fixed standard. Table 2 of the same document then describes the thresholds as the minimum cosine
similarity required for a match. Only the first of those is what the code does, and the difference
matters for these cases.</p>

<p>Because the cutoff is a fraction of each code\'s own best score, the top target of every code
clears it however low that score is. The weakest code on this crosswalk has a best score of %s,
which puts its cutoff at %s at a threshold of 0.990 and %s at 1.0, so its best target passes at
either end of the range in Table 2. There is no value a cosine similarity can fall below and be
rejected.</p>

<p>Step 2 behaves the same way. It keeps each code\'s top N partners ranked by co-occurrence
frequency, so the count is never compared against a value. The smallest highest count in the file
is %s and it is kept on the same terms as the code with %s. Pairs below a frequency of 6 were
removed when the co-occurrence data was built, but that is a property of the data rather than a
decision the algorithm makes.</p>

<p>So neither step rejects a pair for scoring too low. The only two things that make the algorithm
stay quiet are the Step 3 chapter distance filter removing every candidate, or mapping algorithm 3,
the one selected in the Methods, finding no co-occurrence data for the code. The second is what
happens to all 52 on the ICDA-8 side, and between them they are why all 9 on the ICD-10-CA side
come out as false positives.</p>

<p>This is worth settling before the abstain scenarios are tested. A scenario that rejects a pair
for a low cosine similarity or a low co-occurrence count would be a change to Step 1 or Step 2
rather than a change to a parameter, and on the numbers above it would cost a large share of the
correct mappings.</p>
',
  weakest_best("10_9"), floor_at("10_9", 0.990), floor_at("10_9", 1.0),
  cooc_span("10_9", "low"), cooc_span("10_9", "high"))

writeLines(paste0(style, intro, outcome, groups, t10, t8, figs, thresh, mech),
           file.path(OUT, "unmatched_codes_handout.html"))
cat("wrote results/unmatched_codes_handout.html\n")
