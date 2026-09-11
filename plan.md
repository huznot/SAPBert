# Plan after the Sep 9 meeting

Ground rules she set:
- No more Word documents. Spreadsheets.
- She wants actual counts, not just summary statistics.
- Overall numbers first, category breakdowns only after she decides where to drill.
- Weekly Zoom from next week. Send her your schedule.

## 1. The 52 codes with no ICDA-8 match  (her priority for this week)

DONE. All three parts are in `results/review/52_codes_no_icda8.xlsx`, written by
`scripts/52_icda8_excluded_list.R`.

1. Validation data checked. 6 marked wrong, 46 left correct.
2. Chapters they come from. On the main sheet and totalled on the chapters sheet.
3. Top 3 ICDA-8 codes each one co-occurs with.

Point 3 came back empty and that is the finding. None of the 52 appear anywhere
in `icd_8_9_co_occurrence_3d.xlsx`, zero of 52. All 9 of the ICD-10-CA excluded
codes do appear in theirs. So the ICDA-8 pipeline staying silent on all 52 is
missing data, not the method knowing when to stop. Sheet 3, "co-occurrence
coverage", puts the two tracks side by side.

## 2. Frequencies for the thresholds, not just summary stats

DONE. `results/review/absolute_threshold_frequencies.xlsx`, written by
`scripts/53_threshold_frequencies.R`. No new runs, it lays out the counts that
were already in `absolute_threshold_grid.csv`.

354 ICD-9-CM codes x 2038 ICD-10-CA codes = 721,452 possible pairs, and 937 of
those pairs are marked correct in the validation data.

Sheets: read me, pairs above the cutoff, counts by rule (all 1008 rows), best
rule per model, plain cutoff vs old rule, old rule all counts.

The one she will look at first is "pairs above the cutoff". It answers the
question she asked out loud in the meeting:

| Model | Cutoff | Pairs above it | Share |
|---|---|---|---|
| ClinicalBERT | 0.80 | 551,755 | 76.48% |
| ClinicalBERT | 0.90 | 78,791 | 10.92% |
| SapBERT | 0.60 | 2,204 | 0.31% |
| mpnet | 0.75 | 764 | 0.11% |

She remembered the 0.80 number as "550,000", so that checks out.

Counts at each model's best rule:

| Model | Cutoff | Top N | Mappings | TP | FP | FN | F1 |
|---|---|---|---|---|---|---|---|
| ClinicalBERT | 0.90 | 25 | 585 | 345 | 240 | 592 | 0.453 |
| SapBERT | 0.60 | 10 | 661 | 435 | 226 | 502 | 0.544 |
| mpnet | 0.75 | 20 | 616 | 415 | 201 | 522 | 0.535 |

## 3. ICD-9 to ICDA-8 without the code in the label

DONE. `results/review/absolute_threshold_icda8_codeinlabel.xlsx`, written by
`scripts/54_icda8_code_in_label.R`.

Her expectation holds. Paired across all 112 grid points, mean change in F1 from
taking the code out of the label:

| Track | Model | Mean change | Label only better | Code in label better |
|---|---|---|---|---|
| ICDA-8 | ClinicalBERT | -0.005 | 38 of 112 | 73 |
| ICDA-8 | SapBERT | -0.004 | 16 of 112 | 88 |
| ICDA-8 | mpnet | -0.002 | 53 of 112 | 43 |
| ICD-10-CA | ClinicalBERT | +0.045 | 111 of 112 | 0 |
| ICD-10-CA | SapBERT | +0.031 | 110 of 112 | 1 |
| ICD-10-CA | mpnet | +0.021 | 105 of 112 | 5 |

So on ICDA-8 it makes no real difference, and if anything the code-included
version is a hair better. On ICD-10-CA taking the code out is a clear gain. The
reason is the one she gave: ICD-9-CM and ICDA-8 are both numeric and the numbers
often line up, so the code text is not misleading. ICD-10-CA is alphanumeric, so
the ICD-9 number never matches and the model was scoring that mismatch.

Best settings with counts, ICDA-8 track:

| Model | Label | Mappings | TP | FP | FN | F1 |
|---|---|---|---|---|---|---|
| ClinicalBERT | code then label | 295 | 224 | 71 | 107 | 0.716 |
| ClinicalBERT | label only | 301 | 227 | 74 | 104 | 0.718 |
| SapBERT | code then label | 352 | 260 | 92 | 71 | 0.761 |
| SapBERT | label only | 317 | 244 | 73 | 87 | 0.753 |
| mpnet | code then label | 309 | 237 | 72 | 94 | 0.741 |
| mpnet | label only | 311 | 238 | 73 | 93 | 0.741 |

## 4. Why SapBERT median similarity is so much lower

DONE. `results/review/absolute_threshold_median_gap.xlsx`, written by
`scripts/55_median_similarity_gap.R`. Splits all 721,452 pairs into the 937 the
validation data marks correct and the other 720,515, and scores the two groups
separately.

| Model | Median, correct | Median, everything else | Gap |
|---|---|---|---|
| ClinicalBERT | 0.916 | 0.841 | 0.075 |
| SapBERT | 0.654 | 0.220 | 0.434 |
| mpnet | 0.704 | 0.161 | 0.543 |

The median on its own is not the point. The distance between the two groups is,
and that is where ClinicalBERT is weakest. At each model's own best cutoff:

| Model | Cutoff | Correct pairs kept | Other pairs kept | Share of other pairs |
|---|---|---|---|---|
| ClinicalBERT | 0.90 | 555 of 937 | 78,236 | 10.86% |
| SapBERT | 0.60 | 529 of 937 | 1,675 | 0.23% |
| mpnet | 0.75 | 403 of 937 | 361 | 0.05% |

ClinicalBERT and SapBERT keep about the same number of correct pairs, 555
against 529, but ClinicalBERT drags 47 times as many wrong ones along with them.
A low median is the good sign, it means the pairs that should not match are
being pushed down, which is what leaves room between the two groups.

The reason: ClinicalBERT was trained on clinical notes, so it mostly learns
whether text is medical. Every ICD label is medical, so every pair looks alike
to it and the scores bunch up near the top. SapBERT was trained on pairs of
medical terms that name the same thing, so telling apart two terms that are both
medical is exactly its job, which is the thing being asked here.

## 5. Category breakdown

Waiting on her. She wants it eventually but said "we'll get to that" and that
130 CCS categories is too much data to do all of. Ask which ones before
starting.

## Not asked for, do not spend time on

- Confidence intervals. She raised it, then said the differences are in the
  same ballpark from a statistician's view. Not a task.
- More threshold variants. She said there is nothing else to do on threshold
  evaluation until she decides.
- mpnet and SapBERT going forward. She said the comparison was worth doing but
  she does not want to pursue the non-ClinicalBERT models further. Keep the
  comparison, do not extend it.
