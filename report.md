# Report

## 1. Rebuilding the Pipeline

I rebuilt the original analysis from the source data, keeping the same steps in
the same order, so that later changes could be compared against a working
reproduction rather than against reported figures.

Two differences had to be resolved before the original results could be
reproduced. First, the original validated against the CIHI crosswalk table,
which is not in the shared project folder, so I repointed the validation step
at the validation sheets that are shared. Second, the range of co-occurrence Top
N values differed between the two crosswalks in the original code, running from
1 to 10 for ICD-10-CA and from 5 to 30 for ICDA-8. I swept 3 to 30 for both.
The best ICD-10-CA result occurs at a Top N of 30, which is beyond the point
where the original search stopped, and this is why the F1 score of 0.427 did not
reproduce at first. I also added 0.95 to the list of similarity thresholds,
since a threshold of 1.0 keeps only candidates tied with the highest score.

With both corrections the rebuild reproduces 0.427 and 0.716 exactly.

I also rewrote the chapter alignment lookup and the merging step to work on
whole columns at once rather than row by row. This made the pipeline
approximately 130 times faster with identical output, which is what made the
larger parameter searches described below practical.

## 2. Comparing Language Models

ClinicalBERT was published in 2019. I tested two more recent models in its
place. SapBERT is trained specifically on medical vocabulary and the
relationships between medical terms. The second, all-mpnet-base-v2, is a
general-purpose model with no medical training at all, included as a control.

So that the comparison would be fair, I regenerated the embeddings for all
three models with a single script, using identical text cleaning and identical
settings, and evaluated all of them over the same 112 parameter combinations per
crosswalk.

**Table 1. Best F1 score by embedding model.**

| Model | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| ClinicalBERT (original results) | 0.427 | 0.716 |
| ClinicalBERT (regenerated) | 0.430 | 0.716 |
| SapBERT | 0.530 | 0.821 |
| all-mpnet-base-v2 | 0.533 | 0.769 |

![model comparison](results/plot_f1_accuracy_comparison.png)

**Figure 1. F1 score and accuracy for each embedding model.**

SapBERT is the stronger model, improving F1 by approximately 0.10 on both
crosswalks. The result worth noting is that all-mpnet-base-v2, which has no
medical training, performs as well as SapBERT on ICD-10-CA. Medical
pre-training is therefore not the deciding factor in this task, which was the
first indication that the model was not the main constraint.

---

## 3. Finding Where Correct Mappings Were Lost

This section contains the central finding of the work.

Rather than only measuring the final score, I followed every manually mapped
pair through the pipeline and recorded the stage at which it was lost. A pair
removed at an early stage cannot be recovered by a later one, so the earliest
loss sets an upper limit on everything that follows.

**Table 2. Percentage of manually mapped pairs still present at each stage of
the original pipeline.**

| Stage | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Passed the Step 1 similarity cutoff | 37.1% | 79.5% |
| Present in the Step 2 co-occurrence list | 55.4% | 71.6% |
| Present in the candidate set after Step 3 | 62.6% | 85.2% |
| Reported by the Step 4 mapping algorithms | 42.9% | 78.5% |
| Best F1 obtainable from the candidate set | 0.770 | 0.920 |

The similarity cutoff in Step 1 discards 63% of the correct ICD-10-CA pairs
before any mapping algorithm is applied. Because those pairs are gone, the best
F1 score that Step 4 could possibly have achieved was 0.770, no matter which
mapping algorithm was chosen or how the parameters were tuned.

I then tested whether those discarded pairs were recoverable. Replacing the
relative similarity cutoff with a fixed number of candidates per ICD-9-CM code,
removing the chapter filter, and lengthening the co-occurrence list raises the
best obtainable F1 score from 0.770 to 0.927 for ICD-10-CA and from 0.920 to
0.976 for ICDA-8. The correct answers were present and ranked highly enough to
be retrieved. They were being discarded by the selection rule rather than missed
by the model, and this is the evidence that the model was not the limiting
factor.

I also examined the chapter distance filter in Step 3 on its own. Across 168
matched designs, on both crosswalks, there is not one in which the filter
improves the best obtainable F1 score. It costs 0.016 on ICD-10-CA and 0.005 on
ICDA-8, and it destroys 48 correct ICD-10-CA pairs and 3 correct ICDA-8 pairs
outright. It does remove far more incorrect pairs than correct ones, roughly
5,000 against 21 in a typical design, but the incorrect ones would have been
rejected by Step 4 in any case, whereas a correct pair it deletes is
unrecoverable.

---

## 4. Including the Code Number in the Text

The original analysis combined each ICD code with its label into a single text
string before embedding, so that the model received text such as "250 diabetes
mellitus" rather than "diabetes mellitus". I tested the effect of embedding the
label alone.

The tables below report retrieval results only, measured on the similarity
scores before Step 4 selects anything.

**Table 3. Effect of removing the code number, ICD-9-CM to ICD-10-CA.**

| | Correct target ranked highest | In top 10 | In top 25 | In top 50 |
|---|---|---|---|---|
| ClinicalBERT, code and label | 34.8% ± 2.6 | 27.7% | 35.0% | 43.4% |
| ClinicalBERT, label only | **58.6% ± 2.7** | 33.9% | 41.2% | 48.5% |
| SapBERT, code and label | 79.7% ± 2.1 | 58.5% | 68.3% | 73.5% |
| SapBERT, label only | **90.1% ± 1.6** | 64.8% | 73.8% | 81.5% |
| all-mpnet-base-v2, code and label | 80.9% ± 2.2 | 65.5% | 77.0% | 85.1% |
| all-mpnet-base-v2, label only | **89.9% ± 1.6** | 67.7% | 77.8% | 85.8% |

**Table 4. Effect of removing the code number, ICD-9-CM to ICDA-8.**

| | Correct target ranked highest | In top 10 | In top 25 | In top 50 |
|---|---|---|---|---|
| ClinicalBERT, code and label | **58.3% ± 2.9** | 69.2% | 76.7% | 80.4% |
| ClinicalBERT, label only | 55.0% ± 2.9 | 62.5% | 69.5% | 73.7% |
| SapBERT, code and label | **84.1% ± 2.1** | 93.1% | 95.5% | 97.6% |
| SapBERT, label only | 81.8% ± 2.2 | 90.6% | 92.5% | 94.0% |
| all-mpnet-base-v2, code and label | **77.8% ± 2.4** | 87.3% | 93.1% | 95.8% |
| all-mpnet-base-v2, label only | 74.5% ± 2.5 | 85.2% | 90.9% | 94.9% |

The plus or minus figures are bootstrap standard deviations, obtained as in
Section 8 by resampling the ICD-9-CM codes 2,000 times.

The three gains on ICD-10-CA are large and all sit well clear of zero:
+23.8 ± 2.8 points for ClinicalBERT, +10.4 ± 1.8 for SapBERT and +9.0 ± 1.8 for
all-mpnet-base-v2. The three declines on ICDA-8 are much smaller, -3.3 ± 1.9,
-2.3 ± 1.5 and -3.3 ± 1.6 points, and only the last has an interval excluding
zero. The right reading is that removing the code number clearly helps on
ICD-10-CA and does not help on ICDA-8, rather than that it clearly hurts there.

Because the direction is the same for every model, the explanation lies in the
code sets rather than in any particular model. ICD-9-CM and
ICD-10-CA numbering are unrelated, so the number carries no useful information
and acts as noise. By contrast, 69.8% of ICD-9-CM to ICDA-8 pairs share the same
code number, so for that crosswalk the number is genuinely informative.

ClinicalBERT benefits most, gaining 23.8 percentage points on ICD-10-CA compared
with 10.4 for SapBERT and 9.0 for all-mpnet-base-v2. ClinicalBERT is also the
weakest of the three at reading the labels themselves, which suggests it was
relying on the code number more heavily than the others.

---

## 5. Choosing a Standard Stop Word Dictionary

Stop words are common words such as "the", "of" and "and" that carry little
meaning on their own. The original analysis did not remove them.

I compared four published English dictionaries. The important consideration is
not the size of the dictionary but whether removing its words causes two
different codes to end up with the same cleaned label, since the pipeline cannot
distinguish codes in that situation.

**Table 5. Published stop word dictionaries compared.**

| Dictionary | Words | Codes made identical |
|---|---|---|
| Snowball | 175 | 4 |
| NLTK | 179 | 8 |
| SMART | 571 | 23 |
| stopwords-iso | 1298 | 26 |

Snowball causes the fewest collisions.

Snowball nonetheless has one weakness for this application. It removes the
single letters "a" and "i", so that "vitamin a deficiency" becomes "vitamin
deficiency" and "acute hepatitis a" becomes "acute hepatitis", which loses
exactly the character that distinguishes those codes from their neighbours. It
also removes "no", "not" and "nor", which reverse the meaning of a label.
Retaining those five words leaves a dictionary of 170 words that still alters
2,347 labels, compared with 2,381 for the full list, and causes no collisions at
all.

---

## 6. Applying Stop Word Removal to ClinicalBERT

I reran ClinicalBERT under all four combinations of removing stop words and
removing the code number, using the same parameter search each time, so that the
only difference between the four results is the text given to the model.

**Table 6. Best F1 score for ClinicalBERT under each text preparation.**

| ICD-9-CM to ICD-10-CA | Code number included | Code number removed |
|---|---|---|
| Stop words retained | 0.430 ± 0.015 | 0.465 ± 0.016 |
| Stop words removed | 0.436 ± 0.016 | **0.482 ± 0.016** |

| ICD-9-CM to ICDA-8 | Code number included | Code number removed |
|---|---|---|
| Stop words retained | 0.716 ± 0.022 | 0.718 ± 0.022 |
| Stop words removed | 0.719 ± 0.022 | 0.713 ± 0.023 |

![stop words and code numbers](results/plot_stopwords_codes.png)

**Figure 2. Best F1 score for ClinicalBERT under each combination of text
preparation.**

Each of the four differences was also measured with the paired bootstrap
described in Section 8, which is what separates a real effect from a redraw of
the code set.

**Table 6a. What each text change is worth on ICD-9-CM to ICD-10-CA.**

| Change | Difference in F1 | 95% interval |
|---|---|---|
| Remove the code number, stop words retained | +0.034 ± 0.010 | +0.016 to +0.054 |
| Remove stop words, code number retained | +0.006 ± 0.008 | -0.009 to +0.022 |
| Remove stop words, code number already removed | +0.018 ± 0.009 | +0.001 to +0.036 |
| Both changes together | +0.052 ± 0.010 | +0.033 to +0.070 |

Removing the code number is the change that carries the result, worth 0.034 on
its own with an interval well clear of zero. Removing stop words on its own is
worth 0.006 with an interval that includes zero, so it cannot be claimed as an
improvement from this evidence. Applied on top of the code number removal it is
worth 0.018 with an interval that only just clears zero. The two together give
0.482 against 0.430, and that combined gain is solid, but it is mostly the code
number.

On ICDA-8 none of the four differences clears zero. The result stays between
0.713 and 0.719 whatever is done to the text, and all four differences are
smaller than 0.006 against a standard deviation of 0.022. Nothing in the
cleaning process affects that crosswalk.

I also ran both versions of the Snowball dictionary so that the choice
described in Section 5 could be made on evidence.

**Table 7. Best F1 score with each version of the Snowball dictionary.**

| | Snowball as published (175 words) | Snowball retaining letters and negations (170 words) |
|---|---|---|
| ICD-9-CM to ICD-10-CA | 0.439 ± 0.016 | 0.436 ± 0.016 |
| ICD-9-CM to ICDA-8 | 0.717 ± 0.022 | 0.719 ± 0.022 |

The two versions differ by 0.003 on ICD-10-CA and 0.002 on ICDA-8, in opposite
directions, against standard deviations of 0.016 and 0.022. That is at the limit
of what this data can resolve and neither version can be called better. The
reason to prefer the one that retains the five words is therefore not the F1
score but the labels themselves, since the published version silently removes
the letter that identifies several vitamin deficiency and hepatitis codes.

---

## 7. What the Similarity and Co-occurrence Numbers Look Like

Steps 1 and 2 of the original pipeline each cut a list short using a number.
Step 1 uses the similarity score and Step 2 uses how often two codes appear
together in the health records.

### The similarity score

For each of the 354 ICD-9-CM codes I took its best similarity score against any
target code.

**Table 8. Best similarity score available for each ICD-9-CM code,
ClinicalBERT.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Lowest | 0.853 | 0.865 |
| Middle value | 0.930 | 0.968 |
| Codes whose best match is a word for word identical label | 0 | 98 |

![max similarity distribution](results/plot_freq_dist_max_similarity.png)

**Figure 3. How the best available similarity score is spread across the 354
ICD-9-CM codes.**

Two things follow from this.

The first is that every code has a best match above 0.85, so there is no single
score that separates good matches from bad ones. A score of 0.90 might be the
best match one code has and only the tenth best another code has. This is why
the original analysis compared each candidate against that code's own best score
instead of using one fixed cutoff, and that was the right decision.

The second is that ICDA-8 scores higher across the board. Of the 354 ICD-9-CM
codes, 98 have an ICDA-8 code whose label is word for word identical, so those
mappings need no interpretation at all. Not one ICD-9-CM code has an identical
label to an ICD-10-CA code, because ICD-10-CA was rewritten while ICDA-8 was
not. This is the simplest explanation for why every ICDA-8 result in this
document is higher than the matching ICD-10-CA result. It is a difference
between the code sets, not between the methods.

### How often two codes appear together

**Table 9. How often the most frequent partner of each ICD-9-CM code appears
alongside it.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Codes with co-occurrence data | 347 of 354 | 270 of 354 |
| Codes with none | 7 | 84 |
| Middle value | 1,296 | 326 |
| Highest | 238,548 | 19,099 |

![co-occurrence distribution](results/plot_freq_dist_top_cooccurrence.png)

**Figure 4. How the co-occurrence count of each code's most frequent partner is
spread, on a scale where each step is ten times the last.**

The counts differ enormously from one code to another, from single figures up to
238,548. A count of 300 is therefore a strong signal for one code and a weak one
for another, and the raw count is not comparable between codes. A count
carries meaning only relative to the other counts for the same ICD-9-CM code.

The other point is that 84 ICD-9-CM codes, close to a quarter, have no ICDA-8
co-occurrence data at all. For those codes Step 2 of the original pipeline
contributes nothing and the mapping rests entirely on the labels.

---

## 8. How Much These Numbers Move

Every score in this document is a single number measured on one set of 354
ICD-9-CM codes. That leaves the question of how much of a difference is a real
difference and how much is the particular codes that happen to be in the set.

I measured this with a bootstrap. The idea is to treat the 354 codes as a sample
and ask what would have happened if the sample had been slightly different. I
draw 354 codes at random from the set, with replacement, so some codes appear
twice and some not at all, recompute the score on that draw, and repeat 2,000
times. The spread of those 2,000 scores is the standard deviation of the
reported score. Codes are drawn whole, with all of their pairs, because two
pairs belonging to the same ICD-9-CM code are not independent of one another.

**Table 10. Best F1 score with its bootstrap standard deviation and 95%
interval.**

| Model | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| ClinicalBERT (original results) | 0.427 ± 0.015 (0.397 to 0.456) | 0.716 ± 0.022 (0.671 to 0.758) |
| ClinicalBERT (regenerated) | 0.430 ± 0.015 (0.401 to 0.460) | 0.716 ± 0.022 (0.671 to 0.758) |
| SapBERT | 0.530 ± 0.015 (0.501 to 0.560) | 0.821 ± 0.021 (0.780 to 0.863) |
| all-mpnet-base-v2 | 0.533 ± 0.015 (0.503 to 0.561) | 0.769 ± 0.023 (0.725 to 0.815) |

First row: the score is 0.427, a redraw of the code set typically moves it by
0.015, and 95% of the redraws landed between 0.397 and 0.456. The standard
deviation is about 0.015 on ICD-10-CA and 0.022 on ICDA-8, larger there because
that crosswalk has fewer codes, 302 against 345 after exclusions.

These intervals are wider than a comparison between two models needs. Both
models are scored on the same codes, so a draw holding hard codes pulls both
down together and the gap between them barely moves. Table 11 measures that gap
directly, subtracting one model's score from the other's inside each of the
2,000 draws.

**Table 11. Difference in F1 against the original ClinicalBERT result, measured
inside each bootstrap draw.**

| Comparison | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| ClinicalBERT regenerated - original | +0.003 ± 0.005 (-0.006 to +0.014) | 0.000 |
| SapBERT - original | +0.103 ± 0.012 (+0.079 to +0.126) | +0.106 ± 0.024 (+0.062 to +0.155) |
| all-mpnet-base-v2 - original | +0.106 ± 0.012 (+0.084 to +0.129) | +0.054 ± 0.025 (+0.009 to +0.103) |
| all-mpnet-base-v2 - SapBERT | +0.003 ± 0.010 (-0.016 to +0.022) | -0.052 ± 0.019 (-0.090 to -0.018) |
| SapBERT, filler stripped - base | +0.004 ± 0.007 (-0.009 to +0.017) | -0.016 ± 0.013 (-0.041 to +0.010) |
| all-mpnet-base-v2, filler stripped - base | -0.010 ± 0.007 (-0.023 to +0.002) | -0.004 ± 0.010 (-0.023 to +0.015) |

Second row: SapBERT scores 0.103 higher than the original result, the gap moves
by 0.012 between redraws, and 95% of the redraws put it between 0.079 and 0.126.
The whole interval sits above zero, so the gain does not depend on which codes
were in the set.

Three things follow.

The gains from SapBERT and all-mpnet-base-v2 are larger than the measurement
noise. SapBERT is ahead of the original result in 2,000 of the 2,000 draws on
both crosswalks. These are not close calls.

The difference between regenerating ClinicalBERT and the original ClinicalBERT
result is not distinguishable from zero, which is what should happen, since it
is the same model on the same task and the regeneration is only a check that the
reproduction is faithful.

The difference between SapBERT and all-mpnet-base-v2 on ICD-10-CA, 0.003, is far
smaller than the noise, so those two should be reported as tied rather than
ranked. On ICDA-8 SapBERT is ahead of all-mpnet-base-v2 by 0.052, which is
larger than the noise. Filler word stripping does not clear zero anywhere.

The same treatment is applied to the text preparation results in Section 6 and
to the retrieval results in Section 4, and it changes how two of those should be
read. Removing stop words on its own is worth 0.006 on ICD-10-CA against a
standard deviation of 0.008, so it cannot be claimed as an improvement, and the
declines from removing the code number on ICDA-8 are smaller than their own
intervals for two of the three models. Every reported difference and its
interval is in `results/bootstrap_deltas.csv` and
`results/bootstrap_top1_deltas.csv`.

The revised pipeline described elsewhere is evaluated by five-fold cross
validation rather than on one fixed split, so its spread can be read directly
from the folds. It scores 0.669 ± 0.040 on ICD-10-CA and 0.841 ± 0.057 on
ICDA-8 across the five folds. Fold-to-fold spread is wider than the bootstrap
spread above because each fold is scored on a fifth of the codes.

---

## 9. Performance Across All 130 CCS Categories

The 354 ICD-9-CM codes are grouped into 130 CCS categories, clinical groupings
such as breast cancer, asthma or other nervous system conditions. I scored each
category separately, from the same run at each model's best setting, to see
whether the overall figure is carried by a few areas.

Category sizes are uneven: 59 of the 130 hold a single ICD-9-CM code and the
largest holds 18. A one-code category is scored on one or two pairs, so its F1
swings between 0 and 1 on a single decision. The code count is shown in every
table and chart below.

**Table 12. Spread of F1 across the 130 categories.**

| | Mean over categories | SD | Median | Q1 to Q3 | Perfect | Zero | Pooled F1 |
|---|---|---|---|---|---|---|---|
| ICD-10-CA, ClinicalBERT | 0.532 | 0.252 | 0.500 | 0.333-0.667 | 17 | 3 | 0.427 |
| ICD-10-CA, SapBERT | 0.632 | 0.248 | 0.600 | 0.482-0.800 | 27 | 2 | 0.530 |
| ICD-10-CA, all-mpnet-base-v2 | 0.635 | 0.235 | 0.620 | 0.500-0.800 | 25 | 1 | 0.533 |
| ICDA-8, ClinicalBERT | 0.775 | 0.295 | 1.000 | 0.632-1.000 | 65 | 8 | 0.716 |
| ICDA-8, SapBERT | 0.877 | 0.233 | 1.000 | 0.800-1.000 | 91 | 4 | 0.821 |
| ICDA-8, all-mpnet-base-v2 | 0.827 | 0.274 | 1.000 | 0.667-1.000 | 79 | 6 | 0.769 |

SapBERT on ICD-10-CA: the average category scores 0.632, half of the categories
fall between 0.482 and 0.800, 27 are mapped perfectly and 2 are missed
completely. SapBERT and all-mpnet-base-v2 are level on ICD-10-CA here as they
are everywhere else, and separate on ICDA-8.

The last column is the figure quoted everywhere else in this document, and it is
lower than the first in every row. Pooled F1 counts every code once; the mean
over categories counts every category once, whatever its size. Two categories,
one holding a single code that is mapped correctly and one holding ten codes of
which five are mapped, give a mean over categories of 0.75 and a pooled F1 of
0.55. The small categories are the easy ones, so weighting them equally flatters
the result. Pooled F1 stays the headline number.

![f1 spread across categories](results/plot_ccs_f1_distribution.png)

**Figure 5. How F1 is spread across the CCS categories.** Each bar counts the
categories scoring in that range; the dashed lines are the averages.

The spread between categories, an SD of about 0.25, is much larger than the
uncertainty on the overall score in Section 8, about 0.015. The differences
between clinical areas are real, not measurement noise, and one number for the
whole crosswalk hides a range running from 0 to 1.

**Table 13. Mean F1 by how many ICD-9-CM codes the category holds.**

| | 1 code | 2 codes | 3 to 4 codes | 5 or more |
|---|---|---|---|---|
| ICD-10-CA, ClinicalBERT | 0.624 | 0.525 | 0.452 | 0.384 |
| ICD-10-CA, SapBERT | 0.719 | 0.600 | 0.566 | 0.508 |
| ICD-10-CA, all-mpnet-base-v2 | 0.732 | 0.602 | 0.552 | 0.508 |
| ICDA-8, ClinicalBERT | 0.856 | 0.799 | 0.700 | 0.621 |
| ICDA-8, SapBERT | 0.920 | 0.903 | 0.831 | 0.790 |
| ICDA-8, all-mpnet-base-v2 | 0.885 | 0.871 | 0.728 | 0.747 |

The columns hold 59, 23, 27 and 21 categories on ICD-10-CA, and 58, 23, 27 and
21 on ICDA-8.

The more codes a category holds, the worse the pipeline does in it. This holds
on both crosswalks and for all three models, with one exception: mpnet on ICDA-8
scores 0.747 in the largest categories against 0.728 in the 3 to 4 group, which
is a reversal of 0.019 between two groups of 27 and 21 categories. With SapBERT
on ICD-10-CA the score falls from 0.719 in the single-code categories to 0.508
in those holding five or more.

A large category is a clinical area ICD-9-CM divided into many closely related
codes, such as the 18 under other nervous system conditions. Their labels are
similar and they appear together in the same records, so both signals the
pipeline uses point at several codes at once. This is the multi-target problem
in a second form, and it is a property of the clinical area rather than of any
one category, which is why the pattern is so regular.

### The charts

One bar per category, ranked best to worst and dealt into three columns so all
130 fit on a page. The bracketed number after each name is how many ICD-9-CM
codes that category holds.

![clinicalbert icd-10-ca](results/plot_ccs_f1_all_10_9_clinicalbert.png)

**Figure 6. F1 for every CCS category, ICD-9-CM to ICD-10-CA, ClinicalBERT.**

![sapbert icd-10-ca](results/plot_ccs_f1_all_10_9_sapbert.png)

**Figure 7. F1 for every CCS category, ICD-9-CM to ICD-10-CA, SapBERT.**

![mpnet icd-10-ca](results/plot_ccs_f1_all_10_9_mpnet.png)

**Figure 8. F1 for every CCS category, ICD-9-CM to ICD-10-CA,
all-mpnet-base-v2.**

On all three the top is mostly single-code cancers and other conditions with one
distinctive label, and the bottom is where ICD-10-CA split or regrouped an
ICD-9-CM code, such as heart valve disorders and other and ill-defined heart
disease. The hard categories are the same categories for every model.

![clinicalbert icda-8](results/plot_ccs_f1_all_8_9_clinicalbert.png)

**Figure 9. F1 for every CCS category, ICD-9-CM to ICDA-8, ClinicalBERT.**

![sapbert icda-8](results/plot_ccs_f1_all_8_9_sapbert.png)

**Figure 10. F1 for every CCS category, ICD-9-CM to ICDA-8, SapBERT.**

![mpnet icda-8](results/plot_ccs_f1_all_8_9_mpnet.png)

**Figure 11. F1 for every CCS category, ICD-9-CM to ICDA-8,
all-mpnet-base-v2.**

The ICDA-8 charts have a long flat top of categories at 1.000, which is the
identity mapping described in Section 7 showing up category by category.

The next two charts subtract one condition from another, so a bar to the right
means the first named did better in that category.

![sapbert minus clinicalbert, icd-10-ca](results/plot_ccs_delta_sapbert_vs_clinicalbert_10_9.png)

**Figure 12. Change in F1 by category, SapBERT minus ClinicalBERT, ICD-9-CM to
ICD-10-CA.**

![sapbert minus clinicalbert, icda-8](results/plot_ccs_delta_sapbert_vs_clinicalbert_8_9.png)

**Figure 13. Change in F1 by category, SapBERT minus ClinicalBERT, ICD-9-CM to
ICDA-8.**

SapBERT is better in 72 categories, level in 42 and worse in 16 on ICD-10-CA,
and better in 44, level in 71 and worse in 14 on ICDA-8. The improvement is
spread across the crosswalk rather than coming from a few categories, so no
clinical area is carrying it. The 16 losses on ICD-10-CA are small, 14 of them
0.167 or less, and the largest, 0.333, is a one-code category where a single
pair decides the score.

### The text preparation arms

The same breakdown was run on the four ClinicalBERT text preparations from
Section 6, to see whether the text cleaning helps in particular clinical areas
or only in aggregate.

**Table 14. Spread of F1 across categories for each ClinicalBERT text
preparation.**

| | Mean over categories | SD | Perfect | Zero | Pooled F1 |
|---|---|---|---|---|---|
| ICD-10-CA, stop words and code number retained | 0.532 | 0.246 | 16 | 3 | 0.430 |
| ICD-10-CA, code number removed | 0.571 | 0.273 | 25 | 4 | 0.464 |
| ICD-10-CA, stop words removed | 0.545 | 0.263 | 20 | 4 | 0.436 |
| ICD-10-CA, both removed | 0.580 | 0.268 | 24 | 4 | 0.482 |
| ICDA-8, stop words and code number retained | 0.775 | 0.295 | 65 | 8 | 0.716 |
| ICDA-8, code number removed | 0.774 | 0.293 | 64 | 8 | 0.718 |
| ICDA-8, stop words removed | 0.774 | 0.298 | 64 | 9 | 0.719 |
| ICDA-8, both removed | 0.769 | 0.304 | 66 | 9 | 0.713 |

![text cleaning, icd-10-ca](results/plot_ccs_delta_textclean_vs_base_10_9.png)

**Figure 14. Change in F1 by category from removing both the stop words and the
code number, ICD-9-CM to ICD-10-CA.**

![text cleaning, icda-8](results/plot_ccs_delta_textclean_vs_base_8_9.png)

**Figure 15. Change in F1 by category from removing both the stop words and the
code number, ICD-9-CM to ICDA-8.**

On ICD-10-CA the cleaning helps 57 categories, changes nothing in 52 and hurts
21, for a mean gain of 0.048 per category. It is a broad shift rather than a few
categories moving a long way, and the number of categories mapped perfectly rises
from 16 to 24.

On ICDA-8 it changes nothing in 115 of the 129 categories, helps 8 and hurts 6,
for a mean of -0.006. The category breakdown gives the same answer as the
aggregate in Section 6: there is nothing there. This is worth stating because a
null result at the aggregate level can always be a few large gains cancelling a
few large losses, and here it is not, the cleaning simply does not reach that
crosswalk.

Figures 6 to 15 are in `results/` under the names shown. Per-category numbers
for all seven conditions, with the counts each score is built from, are in
`results/ccs_all_categories_10_9.csv` and `results/ccs_all_categories_8_9.csv`,
and the category counts behind Figures 12 to 15 are in
`results/ccs_delta_summary.csv`.
