# Report

---

## Summary

The original pipeline converted each ICD code label into a numeric
representation using ClinicalBERT, selected candidate target codes whose
similarity scores were close to the highest score for that code, added the most
frequently co-occurring target codes, removed pairs whose chapters did not
align, and reported whatever passed one of four mapping algorithms. The best F1
score was 0.427 for ICD-9-CM to ICD-10-CA and 0.716 for ICD-9-CM to ICDA-8.

We rebuilt the pipeline in R and reproduced both of those figures exactly. We
then replaced ClinicalBERT with a newer model, SapBERT, which raised the
ICD-9-CM to ICD-10-CA F1 score from 0.427 to 0.530. Further investigation showed
that the language model was not what limited performance. We traced each
manually mapped pair through the pipeline and found that the similarity cutoff
in Step 1 removed 63% of the correct ICD-10-CA pairs before any selection took
place. Once a pair is removed at that step it cannot be recovered later, so no
version of Step 4 could have produced an F1 score above 0.770.

We therefore changed the pipeline to keep a fixed number of candidate codes for
each ICD-9-CM code rather than applying a relative cutoff, and added a second
model that scores the retained candidates before the final selection is made.
Measured on ICD-9-CM codes the pipeline had not seen during training, the F1
score is now 0.668 for ICD-10-CA and 0.840 for ICDA-8. Each mapping is also
given a confidence score, which allows 86% of ICD-9-CM codes to be accepted
without review at a precision of 95%.

---

## 1. Rebuilding the Pipeline in R

The original analysis was written as two Python notebooks. We ported both to R,
keeping the same steps in the same order, so that later changes could be
compared against a working reproduction rather than against reported figures.

Two differences had to be resolved before the original results could be
reproduced. First, the original validated against the CIHI crosswalk table,
which is not in the shared project folder, so we repointed the validation step
at the validation sheets that are shared. Second, the range of co-occurrence Top
N values differed between the two crosswalks in the original code, running from
1 to 10 for ICD-10-CA and from 5 to 30 for ICDA-8. We swept 3 to 30 for both.
The best ICD-10-CA result occurs at a Top N of 30, which is beyond the point
where the original search stopped, and this is why the F1 score of 0.427 did not
reproduce at first. We also added 0.95 to the list of similarity thresholds,
since a threshold of 1.0 keeps only candidates tied with the highest score.

With both corrections the R version reproduces 0.427 and 0.716 exactly.

We also rewrote the chapter alignment lookup and the merging step to work on
whole columns at once rather than row by row. This made the pipeline
approximately 130 times faster with identical output, which is what made the
larger parameter searches described below practical.

---

## 2. Comparing Language Models

ClinicalBERT was published in 2019. We tested two more recent models in its
place. SapBERT is trained specifically on medical vocabulary and the
relationships between medical terms. The second, all-mpnet-base-v2, is a
general-purpose model with no medical training at all, included as a control.

So that the comparison would be fair, we regenerated the embeddings for all
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

**Figure 1. F1 score and accuracy for ClinicalBERT and SapBERT at the best
parameter combination for each model and crosswalk.**

SapBERT is the stronger model, improving F1 by approximately 0.10 on both
crosswalks. The result worth noting is that all-mpnet-base-v2, which has no
medical training, performs as well as SapBERT on ICD-10-CA. Medical
pre-training is therefore not the deciding factor in this task, which was the
first indication that the model was not the main constraint.

During this work we found and corrected a fault in our own code, in which two
conditions whose names differed only in capitalisation were written to the same
file on Windows, so that parallel jobs overwrote one another's results. All
eight conditions were rerun after the fix.

---

## 3. Reviewing the Filler Word List

The cleaning step used a hand-written list of words to remove from labels before
embedding. We compared this list against four published stop word lists, which
are standard word lists distributed with the NLTK and stopwords packages.

Long lists proved unsafe for this task. Two different codes can end up with
identical text once enough words are removed, and when that happens the pipeline
cannot tell them apart at all. The SMART list produces 23 such collisions and
the stopwords-iso list produces 26, largely because they remove single letters
and the words "with" and "without", so that "hepatitis a" and "hepatitis b"
become the same string.

We removed the words "other", "others", "unspecified", "with" and "without"
from the local list, which between them merged 10 codes, and added an automatic
check that stops the pipeline if any word list causes two codes to share the
same cleaned label.

---

## 4. Finding Where Correct Mappings Were Lost

This section contains the central finding of the work.

Rather than only measuring the final score, we followed every manually mapped
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

We then tested whether those discarded pairs were recoverable. Replacing the
relative similarity cutoff with a fixed number of candidates per ICD-9-CM code,
removing the chapter filter, and lengthening the co-occurrence list raises the
best obtainable F1 score from 0.770 to 0.927 for ICD-10-CA and from 0.920 to
0.976 for ICDA-8. The correct answers were present and ranked highly enough to
be retrieved. They were being discarded by the selection rule rather than missed
by the model, and this is the evidence that the model was not the limiting
factor.

We also examined the chapter distance filter in Step 3 specifically, and found
that it removes more correct pairs than incorrect ones. Chapter alignment is
still used, but as one input to the scoring step rather than as a filter that
deletes candidates outright.

---

## 5. The Revised Pipeline

The original pipeline used its similarity scores to do two jobs at once. It
used them to decide which target codes to consider, and it used them to decide
which target codes to report. Section 4 showed that this is where the problem
was, because a code rejected by the first job could never be recovered by the
second.

The revised pipeline splits those two jobs apart.

**Stage 1. Collect the possibilities.** For each ICD-9-CM code we gather a wide
list of target codes that might be correct. The list contains the 10 closest
target codes according to SapBERT, the 10 closest according to
all-mpnet-base-v2, the 10 closest according to ClinicalBERT, and the 50 target
codes that most often appear alongside that ICD-9-CM code in the health records.
No chapter filter is applied. This gives roughly 60 to 80 candidates per
ICD-9-CM code. The list is intentionally longer than it needs to be, because
including a wrong code at this point costs nothing, whereas leaving out a right
one cannot be undone.

**Stage 2. Judge the possibilities.** Each candidate pair is then given a score
by a second model. The model is a gradient boosted tree, which is a standard
method that learns a sequence of yes or no rules from examples. We train it on
the pairs that were mapped manually, so it learns what a correct pair tends to
look like.

To judge a pair, the model is given 52 pieces of information about that pair.
They answer six questions.

1. *How similar are the two labels?* The similarity score from each of the three
   models, where that score ranks among the candidates for this ICD-9-CM code,
   and how it compares with the best score that ICD-9-CM code achieved.

2. *Does the target code agree?* The same three measurements taken the other way
   round. Being close to an ICD-9-CM code is not the same as being its match.
   Take ICD-9-CM 141, malignant neoplasm of tongue. Its second closest ICD-10-CA
   code is C01, malignant neoplasm of base of tongue, and looking back the other
   way, 141 is the closest ICD-9-CM code to C01. That agreement in both
   directions is what a correct pair looks like, and C01 is indeed a correct
   answer. Its sixth closest is C61, malignant neoplasm of prostate, but looking
   back the other way 141 is only the ninth closest ICD-9-CM code to C61, because
   the genuine prostate codes are nearer. The disagreement is the signal that
   C61 is a false match. This turned out to be the most useful single addition,
   as shown in Section 7.

3. *Are the two codes each other's first choice?* Whether each is the other's
   closest match, and by how much the two directions disagree.

4. *Is this target code popular with everything?* How many other ICD-9-CM codes
   are also competing for it. A target code that looks like a good match for
   dozens of different ICD-9-CM codes is usually a vague label rather than a real
   match for any of them.

5. *Do the two codes actually occur together?* The co-occurrence frequency, and
   where that frequency ranks among the candidates for this ICD-9-CM code.

6. *Do the labels and chapters line up?* How many words the two labels share, how
   close the spellings are, and whether the two chapters align, together with a
   figure learned from the training data for how often that particular pair of
   chapters produces a real mapping.

**Deciding what to report.** A pair is reported if its score is above a fixed
cutoff, or if it is close to the best score for that ICD-9-CM code, or if it is
the best score for that code. Both cutoffs are chosen using only the training
codes, never the codes being tested.

**Table 3. F1 score at each stage of development.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Original pipeline, ClinicalBERT | 0.427 | 0.716 |
| Original pipeline, SapBERT | 0.530 | 0.821 |
| Original pipeline, SapBERT, unseen codes | 0.546 | 0.824 |
| Revised pipeline, SapBERT with all-mpnet-base-v2 and ClinicalBERT, unseen codes | **0.668** | **0.840** |

The revised pipeline uses all three embedding models together rather than
choosing between them. Section 7 shows that ClinicalBERT can be removed from
this combination without any loss.

### Is it fair to compare a trained model against the original rules?

The revised pipeline learns from the manually mapped pairs, which the original
pipeline did not do in the same way. This raises a reasonable objection, so it
is worth setting out plainly.

The original pipeline also used the manually mapped pairs. It tested 96
combinations of similarity threshold, Top N and mapping algorithm, and reported
the combination that scored highest against those pairs. That is a smaller
amount of fitting, three settings rather than a trained model, but it is fitting
to the same data. The difference is that the original then reported that best
score on the very pairs used to choose it, whereas every figure for the revised
pipeline is measured on ICD-9-CM codes that were not used in training.

The comparison that settles the question is the middle two rows of Table 3,
where both the original rules and the revised pipeline are given the same
embeddings, the same data and the same held-out test codes.

**Table 4. Original rules and revised pipeline, both measured on unseen codes.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Original rules, unseen codes | 0.546 | 0.824 |
| Revised pipeline, unseen codes | 0.668 | 0.840 |

Two further points support this being a real improvement rather than an artefact
of training. First, the original rules score 0.547 in sample and 0.546 on unseen
codes, so they were barely overfitting to begin with. The gap between 0.427 and
0.546 is the change of embedding model, not the change of evaluation method.
Second, when the cross-validation folds are divided by CCS category instead of
by code, so that whole clinical areas are absent from training, the revised
pipeline scores 0.657 rather than 0.661. It is not simply memorising the
categories it was trained on.

Two honest qualifications belong with this.

The first is that the improvement comes from two changes made at the same time,
a wider candidate list and a trained scoring model, and the experiments in this
document do not separate how much each contributes. The wider candidate list is
what raises the reachable ceiling from 0.770 to 0.927, so a meaningful part of
the gain is likely to come from that rather than from the model.

The second is that the revised pipeline genuinely requires manually mapped
examples in order to run at all, which the original did not. Section 7 shows
that between 100 and 180 mapped codes are enough. The intended use is therefore
to map a portion of the codes by hand and use the pipeline to extend that work
to the remainder, rather than to build a crosswalk from nothing.

It is also worth noting that on ICD-9-CM to ICDA-8 the trained model adds only
0.016. Almost all of the performance on that crosswalk comes from the code
labels themselves, for the reason given in Section 12.

---

## 6. How Performance Is Measured

The figures in the original analysis were calculated on the same pairs that were
used to choose the parameters. This tends to overstate performance, because the
parameter search can settle on values that happen to suit the particular codes
being scored. All results in this document, apart from those reproducing the
original analysis, are now calculated using five-fold cross-validation in which
the folds are divided by ICD-9-CM code, so that no code is ever used for both
training and testing.

Every reported mapping now also carries a confidence score. This allows the
output to be sorted by confidence rather than accepted or rejected as a whole,
so that reviewer time can be directed at the cases that need it.

![precision and coverage](results/plot_precision_coverage.png)

**Figure 2. Precision against recall as the confidence threshold is varied.**

**Table 5. Proportion of ICD-9-CM codes in each category when the confidence
threshold is set to give 95% precision.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Accepted automatically | 86.1% | 84.1% |
| Sent for review, low confidence | 13.6% | 11.9% |
| Sent for review, no candidate found | 0.3% | 4.0% |

We also report how often the correct target appears among the highest scoring
candidates. This separates the question of whether the pipeline can find the
right answer from the question of whether the reporting rule chooses to keep it.

**Table 6. Percentage of ICD-9-CM codes for which a correct target appears among
the highest scoring candidates.**

| | Highest | Top 3 | Top 5 |
|---|---|---|---|
| ICD-9-CM to ICD-10-CA | 93.3% | 97.7% | 98.6% |
| ICD-9-CM to ICDA-8 | 88.4% | 93.7% | 94.7% |

---

## 7. Which Parts of the Pipeline Are Doing the Work

To confirm that each part of the scoring stage is earning its place, we removed
one group of inputs at a time and re-measured performance. A large drop means
the group was important.

**Table 7. F1 score with one group of inputs removed.**

| Group removed | ICD-10-CA | Change | ICDA-8 | Change |
|---|---|---|---|---|
| Nothing removed | 0.669 | | 0.854 | |
| Reverse direction scores | 0.606 | -0.063 | 0.843 | -0.011 |
| Co-occurrence | 0.618 | -0.051 | 0.819 | -0.035 |
| Chapter alignment | 0.618 | -0.051 | 0.851 | -0.003 |
| Word overlap | 0.655 | -0.015 | 0.854 | 0.000 |
| ClinicalBERT | 0.657 | -0.012 | 0.857 | +0.003 |
| all-mpnet-base-v2 | 0.664 | -0.005 | 0.835 | -0.019 |

Scoring each pair from the target code's point of view as well as the ICD-9-CM
code's is the single most useful addition for ICD-10-CA. ClinicalBERT
contributes nothing once SapBERT and all-mpnet-base-v2 are present, and can be
dropped from the scoring stage.

We asked two further questions about the scoring stage.

The first was whether more manually mapped codes would improve results. We
trained on increasing numbers of codes and measured performance each time.

![learning curve](results/plot_learning_curve.png)

**Figure 3. F1 score against the number of ICD-9-CM codes used for training.**

Performance stops improving at between 100 and 180 training codes on both
crosswalks. Additional manually mapped codes of the same kind would not raise
the score.

The second was whether the pipeline is simply memorising clinical areas. We
repeated the cross-validation with the folds divided by CCS category, so that
whole clinical areas are absent from training. This costs 0.004 F1 on ICD-10-CA
and nothing at all on ICDA-8, so the pipeline is not relying on having seen
similar codes.

---

## 8. Including the Code Number in the Text

The original analysis combined each ICD code with its label into a single text
string before embedding, so that the model received text such as "250 diabetes
mellitus" rather than "diabetes mellitus". We tested the effect of embedding the
label alone.

The tables below report retrieval results only, measured on similarity scores
before any scoring or selection takes place.

**Table 8. Effect of removing the code number, ICD-9-CM to ICD-10-CA.**

| | Correct target ranked highest | In top 10 | In top 25 | In top 50 |
|---|---|---|---|---|
| ClinicalBERT, code and label | 34.8% | 27.7% | 35.0% | 43.4% |
| ClinicalBERT, label only | **58.6%** | 33.9% | 41.2% | 48.5% |
| SapBERT, code and label | 79.7% | 58.5% | 68.3% | 73.5% |
| SapBERT, label only | **90.1%** | 64.8% | 73.8% | 81.5% |
| all-mpnet-base-v2, code and label | 80.9% | 65.5% | 77.0% | 85.1% |
| all-mpnet-base-v2, label only | **89.9%** | 67.7% | 77.8% | 85.8% |

**Table 9. Effect of removing the code number, ICD-9-CM to ICDA-8.**

| | Correct target ranked highest | In top 10 | In top 25 | In top 50 |
|---|---|---|---|---|
| ClinicalBERT, code and label | **58.3%** | 69.2% | 76.7% | 80.4% |
| ClinicalBERT, label only | 55.0% | 62.5% | 69.5% | 73.7% |
| SapBERT, code and label | **84.1%** | 93.1% | 95.5% | 97.6% |
| SapBERT, label only | 81.8% | 90.6% | 92.5% | 94.0% |
| all-mpnet-base-v2, code and label | **77.8%** | 87.3% | 93.1% | 95.8% |
| all-mpnet-base-v2, label only | 74.5% | 85.2% | 90.9% | 94.9% |

All three models improve on ICD-10-CA and decline on ICDA-8 when the code number
is removed. Because the direction is the same for every model, the explanation
lies in the code sets rather than in any particular model. ICD-9-CM and
ICD-10-CA numbering are unrelated, so the number carries no useful information
and acts as noise. By contrast, 69.8% of ICD-9-CM to ICDA-8 pairs share the same
code number, so for that crosswalk the number is genuinely informative.

ClinicalBERT benefits most, gaining 23.8 percentage points on ICD-10-CA compared
with 10.4 for SapBERT and 9.0 for all-mpnet-base-v2. ClinicalBERT is also the
weakest of the three at reading the labels themselves, which suggests it was
relying on the code number more heavily than the others.

---

## 9. The Errors That Remain

We examined which ICD-9-CM codes are still mapped incorrectly at the 95%
precision operating point.

**Table 10. Completeness of the automatic mapping for each ICD-9-CM code.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| All manual targets found | 39.5% | 81.6% |
| Some but not all found | 46.2% | 4.4% |
| None found | 14.2% | 14.0% |
| Codes with one target, mapped correctly | 86.7% (n=143) | 88.2% (n=271) |
| Codes with several targets, mapped completely | 6.0% (n=201) | 0.0% (n=22) |

Performance on codes with a single target is essentially the same for both
crosswalks. The entire difference between the two crosswalks comes from how many
codes have several targets. An ICD-9-CM code maps to an average of 2.72
ICD-10-CA codes, and 63.2% of codes have more than one. For ICDA-8 the average
is 1.10 and only 7.9% have more than one. These counts agree with the original
methods document, which reports 127 one-to-one against 218 one-to-many for
ICD-10-CA and 278 against 24 for ICDA-8.

The practical consequence is that the remaining problem is not recognising the
right clinical concept, which the pipeline does well, but deciding how many
target codes a single ICD-9-CM code should expand into. We checked whether the
additional targets of a multi-target code tend to sit in consecutive blocks,
which would suggest a way of predicting them, and found that they partly do.
This has not been implemented.

![per category change](results/plot_ccs_breakdown_10_9.png)

**Figure 4. Change in F1 score by CCS category when ClinicalBERT is replaced by
SapBERT, ICD-9-CM to ICD-10-CA.** Of the 130 categories scored by both models,
114 improved and 16 declined, so the improvement is spread across the
classification rather than driven by a few categories.

---

## 10. Choosing a Standard Stop Word Dictionary

Stop words are common words such as "the", "of" and "and" that carry little
meaning on their own. The original analysis did not remove them. We were asked
to identify a suitable published dictionary and to confirm the choice before
applying it.

We compared four published English dictionaries. The important consideration is
not the size of the dictionary but whether removing its words causes two
different codes to end up with the same cleaned label, since the pipeline cannot
distinguish codes in that situation.

**Table 11. Published stop word dictionaries compared.**

| Dictionary | Words | Codes made identical |
|---|---|---|
| Snowball | 175 | 4 |
| NLTK | 179 | 8 |
| SMART | 571 | 23 |
| stopwords-iso | 1298 | 26 |

Snowball causes the fewest collisions, and the four it does cause involve codes
that are not part of the validation data. The collisions caused by NLTK do
involve validation codes. On that basis we recommend Snowball.

Snowball nonetheless has one weakness for this application. It removes the
single letters "a" and "i", so that "vitamin a deficiency" becomes "vitamin
deficiency" and "acute hepatitis a" becomes "acute hepatitis", which loses
exactly the character that distinguishes those codes from their neighbours. It
also removes "no", "not" and "nor", which reverse the meaning of a label.
Retaining those five words leaves a dictionary of 170 words that still alters
2,347 labels, compared with 2,381 for the full list, and causes no collisions at
all.

---

## 11. Applying Stop Word Removal to ClinicalBERT

We reran ClinicalBERT under all four combinations of removing stop words and
removing the code number, using the same parameter search each time, so that the
only difference between the four results is the text given to the model.

**Table 12. Best F1 score for ClinicalBERT under each text preparation.**

| ICD-9-CM to ICD-10-CA | Code number included | Code number removed |
|---|---|---|
| Stop words retained | 0.430 | 0.465 |
| Stop words removed | 0.436 | **0.482** |

| ICD-9-CM to ICDA-8 | Code number included | Code number removed |
|---|---|---|
| Stop words retained | 0.716 | 0.718 |
| Stop words removed | 0.719 | 0.713 |

![stop words and code numbers](results/plot_stopwords_codes.png)

**Figure 5. Best F1 score for ClinicalBERT under each combination of text
preparation.**

Removing stop words improves the ICD-10-CA result by 0.006 when the code number
is present and by 0.017 when it is not. It has no meaningful effect on ICDA-8.
Removing the code number is the larger of the two changes, worth between 0.035
and 0.046 on ICD-10-CA, which agrees with the retrieval results in Section 8.
The two changes combine, and the best result is obtained with both applied,
giving 0.482 against 0.427 for the original ClinicalBERT configuration.

The ICDA-8 result stays between 0.713 and 0.719 whatever is done to the text.
Nothing in the cleaning process affects that crosswalk.

We also ran both versions of the Snowball dictionary so that the choice
described in Section 10 could be made on evidence.

**Table 13. Best F1 score with each version of the Snowball dictionary.**

| | Snowball as published (175 words) | Snowball retaining letters and negations (170 words) |
|---|---|---|
| ICD-9-CM to ICD-10-CA | 0.439 | 0.436 |
| ICD-9-CM to ICDA-8 | 0.717 | 0.719 |

The two versions score the same to within the noise of the measurement. The
reason to prefer the version that retains the five words is therefore not the F1
score but the labels themselves, since the published version silently removes
the letter that identifies several vitamin deficiency and hepatitis codes.

---

## 12. Distributions of the Two Pipeline Inputs

The original pipeline applies a threshold to two quantities, the cosine
similarity and the co-occurrence frequency. We were asked to describe how those
two quantities are actually distributed rather than assume it.

**Table 14. Highest cosine similarity available for each ICD-9-CM code,
ClinicalBERT, 354 codes on each crosswalk.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Lowest | 0.853 | 0.865 |
| Median | 0.930 | 0.968 |
| Codes whose best match has an identical label | 0 | 98 (27.7%) |

![max similarity distribution](results/plot_freq_dist_max_similarity.png)

**Figure 6. Distribution of the highest cosine similarity available for each
ICD-9-CM code.**

The ICDA-8 distribution is shifted upwards and has a pronounced spike at 1.0.
This is because 98 ICD-9-CM codes have an ICDA-8 code whose label is word for
word identical, so those mappings are found without any inference. No ICD-10-CA
code has an identical label to any ICD-9-CM code, because ICD-10-CA was
rewritten while ICDA-8 was not. This is the clearest single explanation for why
the ICDA-8 crosswalk scores higher throughout this document, and it is a
property of the code sets rather than of the method.

Every code has a best match above 0.85, which explains why a single fixed
similarity cutoff cannot work and why the original analysis used a cutoff
expressed relative to each code's own best score.

**Table 15. Co-occurrence frequency of the most frequent target for each
ICD-9-CM code.**

| | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Codes with co-occurrence data | 347 of 354 | 270 of 354 |
| Codes with none | 7 (2.0%) | 84 (23.7%) |
| Median frequency | 1,296 | 326 |
| Highest frequency | 238,548 | 19,099 |

![co-occurrence distribution](results/plot_freq_dist_top_cooccurrence.png)

**Figure 7. Distribution of the co-occurrence frequency of the most frequent
target for each ICD-9-CM code, on a logarithmic scale.**

The frequencies span five orders of magnitude, so a raw count is not comparable
between one code and another. The scoring stage therefore uses the rank of a
target within its own ICD-9-CM code as well as the raw count. Almost a quarter
of ICD-9-CM codes have no ICDA-8 co-occurrence data at all, and for those codes
Step 2 of the original pipeline contributes nothing.

---

## 13. Applying the Method to Other Code Sets

The aim of this work is a method that another group can apply to their own
crosswalk, rather than a finished mapping or a piece of software. It is
therefore worth stating precisely what a group would need in order to use it.

The method as described above draws on three things: the code labels for both
classifications, linked physician and hospital records from which co-occurrence
frequencies are calculated, and a table aligning the chapters of the two
classifications. The code labels are published and are available to anyone. The
other two are not. Linked administrative health data covering a period when both
classifications were in use is unusual, and a chapter alignment table has to be
built by hand for each pair of classifications.

We therefore measured what the method achieves when those two inputs are
withheld. Each setting uses the same folds, the same scoring model and the same
held-out test codes, so only the available inputs differ.

**Table 16. Performance by what a group has available, measured on unseen
codes.**

| Available | ICD-9-CM to ICD-10-CA | ICD-9-CM to ICDA-8 |
|---|---|---|
| Labels, health records and chapter table | 0.669 | 0.854 |
| Labels and chapter table only | 0.621 | 0.815 |
| Labels only | 0.581 | 0.806 |
| Original pipeline, for reference | 0.427 | 0.716 |

A group with nothing beyond the published code labels of the two
classifications, and enough manually mapped codes to train on, obtains 0.581 and
0.806. Both are well above what the original pipeline achieved with the full
data. The co-occurrence data is worth 0.048 on ICD-10-CA and 0.040 on ICDA-8,
and the chapter table a further 0.040 and 0.009. Neither is essential.

Three requirements remain, and they should be stated plainly.

First, the method needs manually mapped codes to train on. Section 7 shows that
between 100 and 180 are sufficient, and that dividing the folds by CCS category
costs almost nothing, so those codes do not need to cover every clinical area.
The intended use is to map a portion by hand and extend it, not to build a
crosswalk from nothing.

Second, every embedding model used here was trained on English text. A
classification whose labels are in another language would need a model trained
on that language. We have not tested this, and it is the largest untested
assumption in any claim that the method transfers.

Third, the finding in Section 8 that the code number should be excluded from the
text does not generalise automatically. It holds when the two numbering systems
are unrelated, as ICD-9-CM and ICD-10-CA are, and reverses when they are not, as
with ICD-9-CM and ICDA-8. A group applying the method would need to make that
determination for their own pair of classifications rather than adopt our
conclusion.

Run `scripts/29_portability.R`.

---

## Scripts

Each of the following prints its results from saved files and completes in a few
seconds. To run them, open `icd_crosswalk.Rproj` in RStudio, then use the
`source` command shown.

**Table 17. Scripts for reproducing the results in this document.**

| Script | Shows |
|---|---|
| `27_show_results.R` | Every headline figure in this document |
| `23_code_prefix_test.R` | Section 8 |
| `26_stopword_choice.R` | Section 10, including the affected labels |
| `28_stopwords_and_codes.R` | Section 11 |
| `25_frequency_distributions.R` | Section 12 |
| `29_portability.R` | Section 13 |
| `24_show_similarity_matrix.R` | A worked example of a cosine similarity matrix |

For example, `source("scripts/27_show_results.R")`.
