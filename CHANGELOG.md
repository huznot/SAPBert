# Change log

## Abstract

The original pipeline embedded ICD code labels with ClinicalBERT, kept targets
whose cosine similarity was close to each code's best score, added the top-n
most frequently co-occurring targets, dropped pairs whose ICD chapters did not
align, and output whatever passed one of three rules. Best F1 was 0.427 on
ICD-9 to ICD-10-CA and 0.716 on ICD-9 to ICDA-8.

I rebuilt it in R and confirmed those numbers. SapBERT in place of ClinicalBERT
took ICD-10-CA from 0.427 to 0.530. But the model was not the limit. Every
correct target was already in the similarity matrix. The similarity cutoff threw
out 37% of them before ranking, and a dropped pair cannot be recovered. Nothing
downstream could have scored above 0.770.

I widened the candidate step and added a second model that scores what survives.
Held-out F1 is now 0.668 on ICD-10-CA and 0.840 on ICDA-8. Each mapping carries
a confidence score, so 86% of codes auto-accept at 95% precision.

The remaining errors are one problem. Codes with a single correct target are
solved 87-88% of the time on both tracks. Codes with several are solved 6% on
ICD-10-CA and 0% on ICDA-8. 63% of ICD-9 codes have multiple ICD-10-CA targets
against 8% for ICDA-8, so that accounts for the whole gap between tracks.

---

## 1. Reimplementation

- Ported both Python notebooks to R. Same steps, same order.
- The original validated against the CIHI crosswalk table, not in the shared
  folder. Repointed at the validation sheets that are shared.
- Original top-n ranges differed by track: 1-10 for ICD-10-CA, 5-30 for ICDA-8.
  Both now sweep 3-30. The ICD-10-CA optimum is at 30, past where the original
  stopped. That is why 0.427 would not reproduce at first.
- Similarity cutoffs were 0.99, 0.995, 0.999, 1.0. Added 0.95. A cutoff of 1.0
  only keeps exact ties with the best score.
- With both fixes the R port reproduces 0.427 and 0.716 exactly.
- Vectorized the chapter lookup and dropped `rowwise()` from the merge. 130x
  faster, identical output, verified by a diff script.

## 2. Model comparison

Regenerated every model's embeddings with one script so the arms are
comparable. 112 settings per track, eight conditions.

| model | ICD-10-CA F1 | ICDA-8 F1 |
|---|---|---|
| ClinicalBERT (original) | 0.427 | 0.716 |
| ClinicalBERT (regenerated) | 0.430 | 0.716 |
| ClinicalBERT, filler stripped | 0.448 | 0.722 |
| SapBERT | 0.530 | 0.821 |
| SapBERT, filler stripped | 0.534 | 0.806 |
| all-mpnet-base-v2 | 0.533 | 0.769 |
| all-mpnet-base-v2, filler stripped | 0.523 | 0.766 |

- SapBERT is the best single model, +0.10 and +0.11 over ClinicalBERT.
- all-mpnet-base-v2 has no medical training and ties SapBERT on ICD-10-CA.
  Domain pretraining is not the deciding factor.
- Filler stripping helps ClinicalBERT, hurts SapBERT. Off by default.
- Fixed a bug where condition tags differing only by case wrote to the same file
  on Windows and parallel jobs overwrote each other. Reran all eight.

## 3. Filler words

- Compared the hand-written list against NLTK, Snowball, SMART, stopwords-iso.
- Large lists are unsafe. SMART merges 17 codes, stopwords-iso merges 20. They
  strip single letters and words like "with" and "without", so "hepatitis a" and
  "hepatitis b" collapse to the same string.
- Removed `other`, `others`, `unspecified`, `with`, `without` from the local
  list. Those five merged 10 codes.
- Added a check that fails if the list merges any two codes.

## 4. Error analysis

Traced every correct pair through the pipeline.

| stage | ICD-10-CA | ICDA-8 |
|---|---|---|
| correct target in the similarity matrix | 100% | 100% |
| past the similarity cutoff | 37.1% | 79.5% |
| in the top-n co-occurrence list | 55.4% | 71.6% |
| in the candidate pool | 62.6% | 85.2% |
| output by the three rules | 42.9% | 78.5% |
| best F1 reachable from the pool | 0.770 | 0.920 |

- The model is not the bottleneck. Every answer is in the matrix. The cutoff
  discards 37% of correct ICD-10-CA targets before ranking.
- Separate study on candidate generation: fixed top-k per code instead of a
  relative cutoff, no chapter filter, wider co-occurrence list. Ceiling goes
  from 0.770 to 0.927 on ICD-10-CA and 0.920 to 0.976 on ICDA-8.
- The chapter filter discards more correct pairs than incorrect ones. Chapter
  alignment is now a feature, not a filter.

## 5. Retrieve and rerank

Candidates: top 10 per code from each model, plus top 50 co-occurring targets.
No chapter filter.

Reranker: xgboost over 52 features.

- similarity, rank and relative score per model
- the same in reverse, scoring the pair from the target side
- mutual nearest neighbour flags, and the gap between the two directions
- how many ICD-9 codes compete for the same target
- reciprocal rank fusion across the three models
- co-occurrence frequency and rank
- lexical overlap: shared tokens, edit distance, shared prefix
- chapter alignment, plus a smoothed rate for how often that chapter pair
  actually maps, learned from training folds only

Output rule: keep a pair if it clears an absolute threshold, or is within a
fraction of the best score for that code, or is the best score for that code.
Both cutoffs tuned on inner folds.

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| original ClinicalBERT | 0.427 | 0.716 |
| SapBERT, original rules | 0.530 | 0.821 |
| SapBERT, original rules, held out | 0.546 | 0.824 |
| retrieve and rerank, held out | **0.668** | **0.840** |

## 6. Evaluation

- Original numbers were in-sample: the grid picked settings on the same pairs it
  scored against. Everything now uses 5-fold cross validation with folds split
  by ICD-9 code, so no code is in both training and test.
- Added confidence scores. At a 95% precision target:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| auto-accept | 86.1% | 84.1% |
| review, low confidence | 13.6% | 11.9% |
| review, no candidate | 0.3% | 4.0% |

- Added top-k accuracy, which separates whether the answer is reachable from
  whether the output rule keeps it:

| | top-1 | top-3 | top-5 |
|---|---|---|---|
| ICD-10-CA | 93.3% | 97.7% | 98.6% |
| ICDA-8 | 88.4% | 93.7% | 94.7% |

## 7. Ablation

One feature group removed at a time, held out.

| removed | ICD-10-CA | delta | ICDA-8 | delta |
|---|---|---|---|---|
| nothing | 0.669 | | 0.854 | |
| reverse direction | 0.606 | -0.063 | 0.843 | -0.011 |
| co-occurrence | 0.618 | -0.051 | 0.819 | -0.035 |
| chapter | 0.618 | -0.051 | 0.851 | -0.003 |
| lexical | 0.655 | -0.015 | 0.854 | -0.000 |
| ClinicalBERT | 0.657 | -0.012 | 0.857 | +0.003 |
| rank fusion | 0.663 | -0.006 | 0.838 | -0.016 |
| mpnet | 0.664 | -0.005 | 0.835 | -0.019 |

- Reverse-direction features are the largest single contributor on ICD-10-CA.
- ClinicalBERT contributes nothing alongside SapBERT and mpnet. Can be dropped.
- `rank:pairwise` and `rank:ndcg` are worse on ICD-10-CA. `rank:ndcg` is
  slightly better on ICDA-8, 0.862 against 0.854.
- Learning curve: plateaus at 100-180 training codes on both tracks. More
  manually mapped pairs of the same kind would not help.
- Category holdout: folds split by CCS category so whole clinical areas are
  unseen. Costs 0.004 on ICD-10-CA, 0.000 on ICDA-8. Not memorizing categories.
- Retrieval sensitivity: 24 settings, results flat. Nothing is finely tuned.

## 8. Code numbers in the model input

Labels were embedded as code plus label. Tested label only.

| | top-1 | recall@10 | recall@25 | recall@50 |
|---|---|---|---|---|
| SapBERT, code + label, ICD-10-CA | 79.7% | 58.5% | 68.3% | 73.5% |
| SapBERT, label only, ICD-10-CA | **90.1%** | 64.8% | 73.8% | 81.5% |
| mpnet, code + label, ICD-10-CA | 80.9% | 65.5% | 77.0% | 85.1% |
| mpnet, label only, ICD-10-CA | **89.9%** | 67.7% | 77.8% | 85.8% |
| SapBERT, code + label, ICDA-8 | **84.1%** | 93.1% | 95.5% | 97.6% |
| SapBERT, label only, ICDA-8 | 81.8% | 90.6% | 92.5% | 94.0% |

- ICD-10-CA gains 10 points of top-1. ICD-9 and ICD-10 numbering are unrelated,
  so the number was noise.
- ICDA-8 loses slightly. 69.8% of ICDA-8 pairs share the source code's number,
  so there the number is signal.
- These are retrieval numbers, SapBERT alone, before reranking. Not the same
  measurement as the 93.3% in section 6, which is the full pipeline.
- Label-only matrices are generated but not yet run through the pipeline.

## 9. Remaining errors

Completeness per source code at the 95% precision point.

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| all targets found | 39.5% | 81.6% |
| some but not all | 46.2% | 4.4% |
| none | 14.2% | 14.0% |
| single-target codes solved | 86.7% (n=143) | 88.2% (n=271) |
| multi-target codes solved | 6.0% (n=201) | 0.0% (n=22) |

- ICD-9 codes average 2.72 ICD-10-CA targets (max 13, 63.2% multi) and 1.10
  ICDA-8 targets (max 6, 7.9% multi). The original methods document counts the
  same: 127 one-to-one against 218 one-to-many on ICD-10-CA, 278 against 24 on
  ICDA-8.
- Single-target performance is equal across tracks. The track gap is entirely
  the proportion of multi-target codes.
- Checked whether the extra targets of a multi-target code sit in contiguous
  blocks. They partly do. Not implemented.
