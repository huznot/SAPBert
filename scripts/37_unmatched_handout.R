# builds a standalone handout of the unmatched code descriptives, to send to
# someone who does not have the repo open. reads the csvs written by
# 36_unmatched_descriptives.R and 35_unmatched_codes.R, so run those first.
#
# covers what was asked for in the 27 aug meeting, in her order: co-occurrence
# for the codes with no match, cosine similarity against every target code with
# the average, minimum, maximum, 25th and 75th, the codes that do have a match
# as the comparison group, the groups together and then each code on its own,
# what the algorithm emits for them, and whether a threshold on either signal
# would separate them.
#
# house style, all asked for by name: no em dashes, no date in the heading,
# short noun phrase headings, portrait letter with one inch margins, plain
# words, and no recommendations. the methodology is the pi's, so this describes
# and does not advise. no "as expected" framing either, nothing was predicted.
#
# the html is one single quoted R string, so apostrophes need a backslash.
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

sim   <- read.csv(out_path("unmatched_similarity_by_code.csv"), stringsAsFactors = FALSE)
cooc  <- read.csv(out_path("unmatched_cooccurrence_by_code.csv"), stringsAsFactors = FALSE)
pairs <- read.csv(out_path("unmatched_top_cooccurrence_pairs.csv"), stringsAsFactors = FALSE)
summ  <- read.csv(out_path("unmatched_descriptives_summary.csv"), stringsAsFactors = FALSE)
hand  <- read.csv(out_path("unmatched_code_handling.csv"), stringsAsFactors = FALSE)

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

# ---- the two groups together, with the codes that do have a match alongside ----
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
           "Highest score, lowest", "25th", "Median", "75th", "Highest score, highest"),
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
n_all <- function(track) n_of(track, "no reference match") + n_of(track, "has a reference match")

# the lowest maximum in the crosswalk, and the cutoff it faces at a given
# threshold. together these show the cutoff is relative
weakest_best <- function(tk) f3(min(sim$max[sim$model == MODEL & sim$track == tk]))
floor_at <- function(tk, thr) f3(thr * min(sim$max[sim$model == MODEL & sim$track == tk]))

# middle value of the highest co-occurrence count, per group
cooc_med_raw <- function(tk, grp) {
  median(cooc$max[cooc$track == tk & cooc$group == grp & !is.na(cooc$max)])
}
cooc_med <- function(tk, grp) fq(cooc_med_raw(tk, grp))
n_above_mid <- function(tk) {
  sum(cooc$track == tk & cooc$group == "no reference match" &
        !is.na(cooc$max) & cooc$max > cooc_med_raw(tk, "has a reference match"))
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

# what the pipeline emits for these codes, from 35_
h <- function(tk, col) hand[[col]][hand$track == tk & hand$model == MODEL]

# the cutoff that would exclude every unmatched code, and what it costs
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
  h2 { font-size: 11.5pt; margin: 18px 0 6px 0; page-break-after: avoid; }
  p { margin: 7px 0; }
  table { border-collapse: collapse; width: 100%; font-size: 7pt; line-height: 1.15; margin: 8px 0;
          page-break-inside: avoid; }
  th { text-align: left; border-bottom: 1px solid #000; padding: 2px 5px 3px 5px; }
  td { border-bottom: 0.5px solid #ccc; padding: 1px 5px; vertical-align: top; }
  td.num { text-align: right; white-space: nowrap; }
  table.grp td.txt { white-space: nowrap; }
</style>
'

intro <- sprintf('
<h1>ICD-9-CM codes with no match in the target system</h1>

<p>Of the %d CCS-related ICD-9-CM codes, 9 have no ICD-10-CA match and 52 have no ICDA-8 match in
the reference crosswalk. For these codes the correct output is no mapping, and any mapping the
algorithm produces for them is a false positive.</p>

<p>The algorithm uses two numbers. One is how often a pair of codes appears on the same patient
record, the other is the cosine similarity between the two code labels, taken here from
ClinicalBERT. Both are described below for the two groups together and then for each code
individually, with the codes that do have a match included for comparison.</p>
', n_all("10_9"))

outcome <- sprintf('
<h2>Algorithm output</h2>

<p>On ICD-9-CM to ICD-10-CA the algorithm produces a mapping for all %d codes that should receive
none, so all %d are false positives and none are correctly rejected. On ICD-9-CM to ICDA-8 it
returns no mapping for all %d, and it also returns no mapping for %d of the %d codes that do have a
valid match.</p>
',
  h("10_9", "wrongly_mapped"), h("10_9", "wrongly_mapped"), h("8_9", "correctly_silent"),
  h("8_9", "wrongly_silent"), h("8_9", "n_matched"))

groups <- sprintf('
<h2>Group summary</h2>

<p>Co-occurrence. A code appears in the co-occurrence file only if it shares a record with at least
one target code, so a code absent from the file has no data rather than a low count. Partners is
the number of different target codes a code appears with, and the last five columns describe each
group by the highest count reached by any one of its codes. Pairs below a frequency of 6 were
removed when the co-occurrence data was built, and the lowest count present is 7.</p>

%s

<p>For ICD-10-CA the middle value of the highest count is %s in the group with no match and %s in
the group that has one. For ICDA-8 none of the 52 codes with no match appear in the file.</p>

<p>Cosine similarity. Each ICD-9-CM code is scored against every target code, 2038 for ICD-10-CA
and 858 for ICDA-8, so the average column covers all of those scores across all codes in the group.
The last five columns describe each group by each code\'s single highest score, which is the score
that determines whether a mapping is made.</p>

%s

<p>On both crosswalks the group with no match sits inside the range of the group that has one.</p>
',
  grp_cooc_tbl(), cooc_med("10_9", "no reference match"),
  cooc_med("10_9", "has a reference match"), grp_sim_tbl())

t10 <- sprintf('
<h2>ICD-9-CM to ICD-10-CA</h2>

<p>Co-occurrence for each of the 9 codes.</p>

%s

<p>Code 338, pain not elsewhere classified, appears with 827 different ICD-10-CA codes, and with
E11, type 2 diabetes mellitus, 13,431 times. Code 515 reaches 5,194 and code 239 reaches 4,862. %d
of the 9 exceed %s, the middle value of the highest count across the %d codes that do have a match.
Two are low: code 445, atheroembolism, has 2 partners and a highest count of 16, and code 339,
other headache syndromes, has 58 partners and a highest count of 54.</p>

<p>The codes with the highest counts are ones recorded across a wide range of patients, such as
pain, sleep disorders and neoplasms of unspecified nature.</p>

<p>Cosine similarity for each of the 9 codes. Each row summarises that code\'s 2038 scores.</p>

%s

<p>The highest score for these 9 runs from %s, against %s for the %d codes that do have a match. In
8 of the 9 the closest ICD-10-CA code describes a different condition. Code 175, malignant neoplasm
of male breast, is closest to C30, malignant neoplasm of nasal cavity and middle ear, at 0.955, and
code 515, postinflammatory pulmonary fibrosis, is closest to N72, inflammatory disease of cervix
uteri. These pairs share wording rather than clinical meaning. The exception is code 338, whose
closest target R52 carries the same label, while the crosswalk records no match for it.</p>
',
  cooc_tbl("10_9"), n_above_mid("10_9"), cooc_med("10_9", "has a reference match"),
  n_of("10_9", "has a reference match"), sim_tbl("10_9"),
  rng("10_9", "no reference match"), rng("10_9", "has a reference match"),
  n_of("10_9", "has a reference match"))

t8 <- sprintf('
<h2>ICD-9-CM to ICDA-8</h2>

<p>None of the 52 codes appear in the co-occurrence file, so there is no co-occurrence data to
describe for them. The algorithm requires co-occurrence data before accepting a pair, which is why
it returns no mapping for all 52. The same gap covers 32 of the %d ICD-9-CM codes that do have a
valid ICDA-8 match, and accounts for the %d mappings it does not produce.</p>

<p>Cosine similarity. The highest score for the 52 runs from %s, against %s for the %d codes that
do have a match. The four highest and the four lowest are below.</p>

%s
',
  n_of("8_9", "has a reference match"), h("8_9", "wrongly_silent"),
  rng("8_9", "no reference match"), rng("8_9", "has a reference match"),
  n_of("8_9", "has a reference match"), sim_tbl("8_9", ends = 4))

thresh <- sprintf('
<h2>Threshold scenarios</h2>

<p>A similarity cut-off high enough to exclude all 9 ICD-10-CA codes would have to sit above %s,
which also excludes %d of the %d codes that are correct. On ICDA-8 the equivalent cut-off is %s,
which excludes %d of %d. A co-occurrence cut-off high enough to exclude the 9 would have to sit
above %s, which excludes %d of the %d matched codes present in the file. Requiring both conditions
does not change this, since the 9 ICD-10-CA codes pass both, and the 52 ICDA-8 codes are excluded
by the absence of co-occurrence data rather than by a threshold.</p>
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
multiplied by the threshold value, which makes it relative to the code being mapped. Table 2 of the
same document describes the thresholds as the minimum cosine similarity required for a match. The
code implements the Step 1 definition.</p>

<p>Because the cutoff is a fraction of each code\'s own highest score, the highest scoring target of
every code clears it. The lowest maximum on this crosswalk is %s, which gives that code a cutoff of
%s at a threshold of 0.990 and %s at 1.0, and its highest scoring target clears both.</p>

<p>Step 2 selects each code\'s top N partners ranked by co-occurrence frequency, so the count is not
compared against a value. The lowest maximum count in the file is %s, and it is selected on the
same terms as the code with %s.</p>

<p>Under these definitions a pair is not excluded for a low score on either signal. A code produces
no mapping when the Step 3 chapter distance filter removes all of its candidates, or when mapping
algorithm 3, the combination selected in the Methods, finds no co-occurrence data for it. The
second applies to all 52 ICDA-8 codes.</p>
',
  weakest_best("10_9"), floor_at("10_9", 0.990), floor_at("10_9", 1.0),
  cooc_span("10_9", "low"), cooc_span("10_9", "high"))

writeLines(paste0(style, intro, outcome, groups, t10, t8, thresh, mech),
           out_path("unmatched_codes_handout.html"))
cat("wrote results/unmatched_codes_handout.html\n")
