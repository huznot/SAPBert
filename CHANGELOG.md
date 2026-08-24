# Change log

Every change made to the automatic ICD mapping pipeline, from the original
2024 version to now.

## Abstract

The original pipeline mapped ICD-9-CM codes to ICD-10-CA and ICDA-8 at the
three-character level by embedding the code labels with ClinicalBERT, keeping
every target whose cosine similarity was within a fraction of that code's best
score, adding the top-n most frequently co-occurring targets, dropping any pair
whose ICD chapters did not align, and emitting whatever carried at least one of
three flags: best similarity, best co-occurrence, or present in both lists.
Performance was reported over a grid of threshold, top-n and flag settings and
scored against a manual crosswalk. Its best setting reached F1 0.427 on
ICD-9 to ICD-10-CA and 0.716 on ICD-9 to ICDA-8, both measured on the same
data the settings were chosen from.

The work since then reimplemented that pipeline in R, verified it reproduces
the original numbers exactly, and then tested where it was losing accuracy.
Substituting the embedding model helped (SapBERT lifted ICD-10-CA F1 from 0.427
to 0.530) but was not the main problem. An error analysis showed that every
correct target was present somewhere in the similarity matrix, yet only 62.6%
of them survived into the candidate set the flag rules chose from, and once a
pair is dropped it can never be recovered. The best score any downstream method
could reach from that candidate set was F1 0.770. Replacing the threshold and
flag rules with a wider top-k retrieval step and a trained gradient-boosted
ranker over 52 features raised held-out F1 to 0.668 on ICD-10-CA and 0.840 on
ICDA-8, evaluated by cross validation grouped on the source code so no code
appears in both training and test. Every mapping now carries a confidence
score, so 86% of ICD-9 codes can be auto-accepted at 95% precision and the
rest routed to a human.

The remaining error is concentrated in one place. Codes with a single correct
target are fully solved 87-88% of the time on both tracks. Codes with several
correct targets are fully solved 6% of the time on ICD-10-CA and 0% on ICDA-8.
Because 63% of ICD-9 codes have more than one ICD-10-CA target and only 8% have
more than one ICDA-8 target, that single failure mode accounts for essentially
the whole gap between the two tracks.

---

## 1. Reimplementation and verification

- Ported the two Jupyter notebooks (`Pipeline_10_9_2024Nov27.ipynb`,
  `Pipeline_8_9_2024Oct30.ipynb`) to R. Same steps, same order, same metric
  definitions.
- Original results were computed against the CIHI crosswalk table, which is not
  in the shared package. Repointed validation at the `Validation_Data.xlsx`
  sheets that are shared, and at the matching exclusion sheets.
- The original swept different top-n ranges per track: `{1,2,3,4,5,10}` for
  ICD-10-CA and `{5,10,15,20,25,30}` for ICDA-8, 96 settings each. Both tracks
  now sweep `{3,5,10,15,20,25,30}`. The ICD-10-CA optimum turns out to sit at
  top_n = 30, three times past where the original grid stopped, which is why
  its reported F1 of 0.427 looked unreproducible at first.
- Similarity thresholds changed from `{0.990, 0.995, 0.999, 1.0}` to
  `{0.95, 0.99, 0.995, 0.999}`. A threshold of 1.0 keeps only exact ties with
  the best score, and nothing below 0.99 had been tested.
- Confirmed the R port reproduces the original ClinicalBERT maxima exactly:
  F1 0.427 (ICD-10-CA) and 0.716 (ICDA-8).
- Rewrote `find_chapter` and `compute_chapter_distance` to work on vectors
  instead of row by row, and removed `rowwise()` from the merge step. About
  130x faster, bit-identical output, verified by a script that diffs old
  against new.

## 2. Embedding model comparison

- Added a Python script that regenerates embeddings and cosine similarity
  matrices for any model, so every arm is produced the same way instead of
  mixing an external pipeline with a local one.
- Ran a full grid (4 thresholds x 7 top-n x 4 flag rules = 112 settings per
  track) for eight conditions rather than one setting per model.

| model | ICD-10-CA F1 | ICDA-8 F1 |
|---|---|---|
| ClinicalBERT (original) | 0.427 | 0.716 |
| ClinicalBERT (regenerated) | 0.430 | 0.716 |
| ClinicalBERT + filler stripping | 0.448 | 0.722 |
| SapBERT | 0.530 | 0.821 |
| SapBERT + filler stripping | 0.534 | 0.806 |
| all-mpnet-base-v2 | 0.533 | 0.769 |
| all-mpnet-base-v2 + filler stripping | 0.523 | 0.766 |

- SapBERT is the best single model. The gain over ClinicalBERT is large
  (+0.10 and +0.11 F1).
- all-mpnet-base-v2 is a general purpose sentence model with no medical
  training and it matches SapBERT on ICD-10-CA. Domain pretraining is not the
  deciding factor on this task.
- Filler word stripping helps ClinicalBERT and hurts SapBERT, so it is off by
  default.
- Found and fixed a bug where two condition tags differing only in case wrote
  to the same file on Windows, so parallel jobs overwrote each other. All eight
  conditions were rerun after the fix.

## 3. Filler word / stopword testing

- Tested the hand written filler list against NLTK, Snowball, SMART and
  stopwords-iso on all three code label sets.
- Large lists are unsafe here. SMART merges 17 codes and stopwords-iso merges
  20, because they remove single letters and words like "with" and "without"
  that carry clinical meaning, so "hepatitis a" and "hepatitis b" become the
  same string.
- Removed `other`, `others`, `unspecified`, `with`, `without` from the local
  list. Those five were merging 10 codes on their own.
- Added a check that fails if the list in use merges any two codes.

## 4. Error analysis, where the accuracy was actually going

Ran a loss accounting for every correct pair in the manual crosswalk.

| stage | ICD-10-CA | ICDA-8 |
|---|---|---|
| correct target exists in the similarity matrix | 100% | 100% |
| survives the similarity threshold | 37.1% | 79.5% |
| present in the top-n co-occurrence list | 55.4% | 71.6% |
| in the candidate pool (either of the above) | 62.6% | 85.2% |
| actually emitted by the flag rules | 42.9% | 78.5% |
| best F1 reachable from that pool | 0.770 | 0.920 |

- The embedding model is not the binding constraint. Every answer is in the
  matrix. The relative threshold and top-n cutoff throw out 37% of correct
  ICD-10-CA targets before anything is ranked, and there is no way to recover
  them later.
- Ran a separate study varying how candidates are generated. Switching from a
  relative similarity threshold to a fixed top-k per code, dropping the chapter
  filter, and widening the co-occurrence list raises the reachable ceiling from
  0.770 to 0.927 on ICD-10-CA and 0.920 to 0.976 on ICDA-8.
- The chapter filter costs more than it saves. Chapter alignment is now a
  feature the model can weigh rather than a hard filter.

## 5. New architecture

Replaced threshold and flag rules with retrieve then rerank.

- Retrieval: keep the top 10 candidates per code from each embedding model plus
  the top 50 co-occurring targets. No hard chapter filter.
- Reranking: an xgboost model scores every candidate pair. 52 features:
  - forward similarity, rank and relative score, per embedding model
  - reverse direction similarity and rank, so a pair is judged from both sides
  - mutual features: whether each is the other's nearest neighbour, the
    geometric mean of both directions, the sum of both ranks, the gap between
    directions
  - how many other ICD-9 codes are competing for the same target
  - reciprocal rank fusion across the three embedding models
  - co-occurrence frequency and rank
  - lexical overlap between labels: token Jaccard, edit distance, shared prefix
    length
  - chapter alignment, plus a smoothed chapter-pair compatibility rate learned
    from the training folds only
- Emission: a pair is output if its score clears an absolute threshold, or is
  within a fraction of the best score for that code, or is the best score for
  that code. Both cutoffs are tuned on inner folds, never on test data.

Results, grouped 5-fold cross validation with folds split by ICD-9 code:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| original ClinicalBERT pipeline (in-sample) | 0.427 | 0.716 |
| SapBERT, same rules (in-sample) | 0.530 | 0.821 |
| SapBERT, same rules (held out) | 0.546 | 0.824 |
| retrieve and rerank (held out) | **0.668** | **0.840** |

Precision 0.742 / recall 0.607 on ICD-10-CA, precision 0.867 / recall 0.815 on
ICDA-8.

## 6. Evaluation changes

- The original numbers were in-sample: the grid picked the best threshold and
  flag combination on the same pairs it was scored against. Everything now
  reports held-out performance from cross validation with folds grouped on the
  source code, so no ICD-9 code contributes to both training and test.
- Added a confidence score per mapping and a precision-at-coverage curve. At a
  95% precision target the triage is:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| auto-accept | 86.1% | 84.1% |
| review, candidates but low confidence | 13.6% | 11.9% |
| review, no plausible candidate | 0.3% | 4.0% |

- Added top-1 and top-k retrieval accuracy, which measures whether the right
  answer is reachable, separately from whether the emission rule keeps it:

| | top-1 | top-3 | top-5 |
|---|---|---|---|
| ICD-10-CA | 93.3% | 97.7% | 98.6% |
| ICDA-8 | 88.4% | 93.7% | 94.7% |

## 7. Ablation and robustness

Held-out F1 with one feature group removed at a time:

| removed | ICD-10-CA | delta | ICDA-8 | delta |
|---|---|---|---|---|
| nothing (full) | 0.669 | | 0.854 | |
| mutual / reverse direction | 0.606 | -0.063 | 0.843 | -0.011 |
| co-occurrence | 0.618 | -0.051 | 0.819 | -0.035 |
| chapter | 0.618 | -0.051 | 0.851 | -0.003 |
| lexical | 0.655 | -0.015 | 0.854 | -0.000 |
| ClinicalBERT arm | 0.657 | -0.012 | 0.857 | +0.003 |
| ensemble fusion | 0.663 | -0.006 | 0.838 | -0.016 |
| mpnet arm | 0.664 | -0.005 | 0.835 | -0.019 |

- The reverse-direction and mutual-neighbour features are the single largest
  contributor on ICD-10-CA.
- ClinicalBERT contributes nothing once SapBERT and mpnet are present.
- Tested `rank:pairwise` and `rank:ndcg` against `binary:logistic`. Ranking
  objectives are worse on ICD-10-CA and slightly better on ICDA-8
  (0.862 vs 0.854).
- Learning curve: performance plateaus at roughly 100-180 training codes on
  both tracks. More manually verified pairs of the same kind would not help.
- Category holdout: folds split by CCS category instead of by code, so whole
  clinical areas are unseen at test time. Costs -0.004 on ICD-10-CA and -0.000
  on ICDA-8, within noise. The model is not memorising categories.
- Retrieval sensitivity: swept 24 retrieval settings. Results are flat across a
  wide range, so the choice of k is not delicately tuned.

## 8. Removing code numbers from the model input

Tested feeding the model only the text label, with the code number stripped
from the string that gets embedded.

| | top-1 hit | recall @ 10 | recall @ 25 | recall @ 50 |
|---|---|---|---|---|
| SapBERT, code + label, ICD-10-CA | 79.7% | 58.5% | 68.3% | 73.5% |
| SapBERT, label only, ICD-10-CA | **90.1%** | 64.8% | 73.8% | 81.5% |
| mpnet, code + label, ICD-10-CA | 80.9% | 65.5% | 77.0% | 85.1% |
| mpnet, label only, ICD-10-CA | **89.9%** | 67.7% | 77.8% | 85.8% |
| SapBERT, code + label, ICDA-8 | **84.1%** | 93.1% | 95.5% | 97.6% |
| SapBERT, label only, ICDA-8 | 81.8% | 90.6% | 92.5% | 94.0% |

- On ICD-10-CA, dropping the code number moves the correct target into first
  place 10 percentage points more often. The code number was acting as noise
  because ICD-9 and ICD-10 numbering are unrelated.
- On ICDA-8 it is slightly worse, because 69.8% of ICDA-8 pairs have the same
  number as the ICD-9 source, so the number is genuinely informative there.
- These label-only matrices have been generated but not yet run through the
  full rerank pipeline.

## 9. Where the remaining error is

Completeness per source code, at the 95% precision operating point:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| all correct targets found | 39.5% | 81.6% |
| some but not all | 46.2% | 4.4% |
| none | 14.2% | 14.0% |
| single-target codes fully solved | 86.7% (n=143) | 88.2% (n=271) |
| multi-target codes fully solved | 6.0% (n=201) | 0.0% (n=22) |

- ICD-9 codes have on average 2.72 ICD-10-CA targets (max 13, 63.2% have more
  than one) and 1.10 ICDA-8 targets (max 6, 7.9% have more than one). This
  matches the original methods document, which counted 127 one-to-one and 218
  one-to-many ICD-9 codes on the ICD-10-CA side, against 278 one-to-one and 24
  one-to-many on the ICDA-8 side.
- Single-target performance is the same on both tracks. The whole track gap is
  the proportion of multi-target codes, not the difficulty of the code systems.
- Checked whether the extra targets of a multi-target code sit in contiguous
  blocks of the target code space, which would allow a block-expansion feature.
  They partly do. Not yet implemented.

## 10. Housekeeping

- Repository made public, MIT licensed.
- Added a script that produces a worked example of the similarity matrix, with
  the manual crosswalk answers marked, for explaining the method
  (`24_show_similarity_matrix.R`).
- Removed the R Markdown report in favour of this document.

## Corrections made along the way

Recording these because they changed reported numbers.

- An early version of the ablation script computed the truth set after
  retrieval instead of before, which dropped the positives retrieval had missed
  out of the recall denominator and inflated F1. A claimed "+0.034 free accuracy
  from fixing retrieval" was withdrawn. Fixed, and the corrected script now
  agrees with the main cross validation to within 0.001.
- Coverage was first reported as 43% / 79%, which was recall against reachable
  pairs. Against the full crosswalk it is 37.5% / 76.1%.
- Review burden was first reported as 57%. That was the share of pairs, not
  codes. It is 14% of codes.
- An early claim that nine codes had no ground truth and were silently skipped
  was wrong. Those codes are on the crosswalk's own exclusion list.

## What is still open

- Run the label-only embeddings through the full rerank pipeline.
- Add block or sibling features aimed specifically at multi-target codes.
- Larger validation set. The current one is 937 ICD-10-CA and 331 ICDA-8
  manually verified pairs over 354 three-digit ICD-9 categories, not the full
  code set. The CIHI crosswalk table used in the original work would roughly
  match what the published LLM translation study evaluated on.
- Wider co-occurrence coverage. Co-occurrence is the second most important
  feature group and it is missing for many pairs.
