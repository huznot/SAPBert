# The four things asked for on 9 September

Spreadsheets, not Word documents, and counts rather than only precision and
recall. Every workbook has a "read me" tab first that says what each column is
and where the numbers come from.

All of this is the ICD-9-CM to ICD-10-CA track unless a sheet says otherwise.
There are **354 ICD-9-CM codes and 2038 ICD-10-CA codes, so 721,452 possible
pairs**, and the validation data marks **937 of those pairs correct**.

### 1. [Frequencies at each cutoff](1_threshold_frequencies.xlsx)

How many of the 721,452 pairs each cutoff calls similar, and the true positives,
false positives and false negatives as raw counts for all 1008 rules tested.
Start on the "pairs above the cutoff" tab:

| Model | Cutoff | Pairs above it | Share of all pairs |
|---|---|---|---|
| ClinicalBERT | 0.80 | 551,755 | 76.48% |
| ClinicalBERT | 0.90 | 78,791 | 10.92% |
| SapBERT | 0.60 | 2,204 | 0.31% |
| mpnet | 0.75 | 764 | 0.11% |

Counts at each model's best rule:

| Model | Cutoff | Top N | Mappings | TP | FP | FN | F1 |
|---|---|---|---|---|---|---|---|
| ClinicalBERT | 0.90 | 25 | 585 | 345 | 240 | 592 | 0.453 |
| SapBERT | 0.60 | 10 | 661 | 435 | 226 | 502 | 0.544 |
| mpnet | 0.75 | 20 | 616 | 415 | 201 | 522 | 0.535 |

The plain cutoff is also put next to the old fraction-of-the-maximum rule, in
counts, on its own tab.

### 2. [The code number in the label, ICDA-8](2_icda8_code_in_label.xlsx)

Taking the code out of the label makes no real difference on the ICDA-8 track,
which is what was expected. Paired across all 112 grid points per model:

| Track | Model | Mean change in F1 | Label only better |
|---|---|---|---|
| ICDA-8 | ClinicalBERT | −0.005 | 38 of 112 |
| ICDA-8 | SapBERT | −0.004 | 16 of 112 |
| ICDA-8 | mpnet | −0.002 | 53 of 112 |
| ICD-10-CA | ClinicalBERT | +0.045 | 111 of 112 |
| ICD-10-CA | SapBERT | +0.031 | 110 of 112 |
| ICD-10-CA | mpnet | +0.021 | 105 of 112 |

ICD-9-CM and ICDA-8 are both numeric and the numbers often line up, so the code
text is not misleading. ICD-10-CA is alphanumeric, so the ICD-9 number never
matches and the models were scoring that mismatch.

### 3. [Why SapBERT's median sits so much lower](3_similarity_gap_by_model.xlsx)

Splits all 721,452 pairs into the 937 marked correct and the other 720,515, and
scores the two groups separately.

| Model | Median, correct | Median, everything else | Gap |
|---|---|---|---|
| ClinicalBERT | 0.916 | 0.841 | 0.075 |
| SapBERT | 0.654 | 0.220 | 0.434 |
| mpnet | 0.704 | 0.161 | 0.543 |

The median on its own is not the useful number, the distance between the two
groups is. At each model's own best cutoff:

| Model | Cutoff | Correct pairs kept | Other pairs kept | Share of other pairs |
|---|---|---|---|---|
| ClinicalBERT | 0.90 | 555 of 937 | 78,236 | 10.86% |
| SapBERT | 0.60 | 529 of 937 | 1,675 | 0.23% |
| mpnet | 0.75 | 403 of 937 | 361 | 0.05% |

ClinicalBERT and SapBERT keep about the same number of correct pairs, 555
against 529, but ClinicalBERT brings 47 times as many wrong ones with them. So a
low median is the good sign, it means the pairs that should not match are being
pushed down.

ClinicalBERT was trained on clinical notes, so it mostly learns whether text is
medical, and every ICD label is medical. SapBERT was trained on pairs of terms
that name the same thing, which is the job being asked here.

### 4. [Counts by CCS category](4_counts_by_ccs_category.xlsx)

All 130 categories with TP, FP and FN, at each model's own plain cutoff, plus a
ClinicalBERT-only tab. This is the raw material for deciding where to drill, not
a drill-down itself.

Worth reading before picking categories. The categories are very uneven:

| Category size | Categories | Share of the 937 correct pairs |
|---|---|---|
| 1 code | 59 (45%) | 15% |
| 2 codes | 23 (18%) | 13% |
| 3 to 4 codes | 27 (21%) | 28% |
| 5 to 9 codes | 18 (14%) | 31% |
| 10 or more codes | 3 (2%) | 12% |

Nearly half the categories hold a single ICD-9 code, so their F1 is computed off
one or two pairs and swings between 0 and 1 for reasons that have nothing to do
with the model. The 21 categories holding five or more codes carry 43% of the
correct pairs and are the only ones where a per-category number means much.

## Also relevant, from earlier

- [The 52 ICD-9 codes with no ICDA-8 match](../review/52_codes_no_icda8.xlsx) —
  the validation check, which chapters they come from, and the top three ICDA-8
  codes each one co-occurs with. That last part came back empty, and that is the
  finding: none of the 52 appear anywhere in the ICD-9 to ICDA-8 co-occurrence
  file, while all 9 of the ICD-10-CA excluded codes do appear in theirs. So the
  pipeline staying silent on all 52 is missing data, not the method knowing when
  to stop.
- [The reviewed nine codes](../review/9_codes_review_UPDATED.xlsx)
- [249 and 239 in detail](../review/e13_d48_review.xlsx)

## Rebuilding these

```bash
Rscript scripts/53_threshold_frequencies.R    # 1
Rscript scripts/54_icda8_code_in_label.R      # 2
Rscript scripts/55_median_similarity_gap.R    # 3
Rscript scripts/56_category_counts.R          # 4
```

`53` reads `results/review/absolute_threshold_grid.csv` and recomputes nothing.
The other three reload the similarity matrices, so they take a few minutes each.
