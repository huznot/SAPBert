# Automated ICD Crosswalk Construction Using Embedding Similarity, Diagnosis Co-occurrence, and Learned Reranking

**Draft for discussion — not for circulation.** Numbers are current as of the
latest run and every figure below is reproducible from this repository. Open
items are marked `[TODO]`.

Authors: Muhammad Huzaifa Irfan, [others], Lisa M. Lix

---

## Abstract

**Background.** Crosswalk tables mapping diagnosis codes between revisions of
the International Classification of Diseases (ICD) are costly to build by hand
but necessary whenever a health system adopts a new coding standard. Recent
work in this group found that GPT-4 translated ICD-10-CA codes to ICD-9-CM
with 84.7% accuracy at the three-character level but was judged insufficient
for deployment, in part because the model gives no reliable indication of
which of its outputs to trust.

**Objective.** To evaluate whether a task-specific pipeline combining code
label embeddings, historical diagnosis co-occurrence, and a learned reranker
can construct ICD crosswalks accurately enough for practical use, and to
determine what limits its accuracy.

**Methods.** ICD-9-CM codes were mapped to ICD-10-CA and to ICDA-8 at the
three-character level. Candidate target codes were retrieved using cosine
similarity from three sentence embedding models and top-N diagnosis
co-occurrence, then reranked by a gradient-boosted model using similarity,
rank-fusion, mutual (forward and reverse) agreement, co-occurrence, learned
chapter compatibility, and lexical overlap features. Performance was estimated
by five-fold cross-validation grouped by ICD-9-CM code, with the reranker,
retrieval settings, and decision thresholds all fitted inside training folds.
The existing rule-based pipeline was evaluated under the identical protocol.

**Results.** Held-out F1 was 0.668 for ICD-9-CM to ICD-10-CA and 0.840 for
ICD-9-CM to ICDA-8, against 0.546 and 0.824 for the rule-based baseline. A
correct target was ranked first for 93.3% and 88.4% of codes respectively.
Requiring 95% precision, 43.2% and 79.0% of the reference crosswalk was
reproduced with no human review. An error decomposition showed that candidate
retrieval, not the embedding model, was the binding constraint: every correct
target was present in the similarity matrices, but the original relative
similarity threshold admitted only 62.6% of them, capping attainable F1 at
0.770. Performance was unchanged when whole clinical categories were held out
(-0.004, -0.000) and plateaued at approximately 150 training codes.

**Conclusions.** A task-specific pipeline reproduces most of a manually
constructed ICD crosswalk and, unlike a general-purpose language model,
attaches calibrated confidence to each mapping so that high-precision
mappings can be accepted automatically and the remainder routed for review.
Candidate retrieval, not label quantity or embedding choice, is the principal
limitation.

---

## 1. Background

ICD coding standards are revised periodically, and each revision strands the
records coded under the previous one. Bridging them requires a crosswalk
mapping every legacy code to its modern equivalent. These are built manually
by clinical coders, are expensive, and must be rebuilt for each pair of
systems and each national variant.

Two automation strategies have been examined in this group. Monchka et al.
evaluated GPT-4 across nine prompting strategies, reporting 48.3% full-code
and 84.7% three-character accuracy for ICD-10-CA to ICD-9-CM translation, and
concluded that performance was insufficient for deployment. That study also
noted that publicly available crosswalks may appear in the model's training
data, so some of the observed accuracy may reflect memorisation.

The alternative examined here combines two signals available to any health
system undertaking a migration: the text of the code labels, and how often
codes co-occur in records coded under both systems. Neither signal requires an
external model or transmitting data to a third party, and co-occurrence counts
derived from local records cannot be contaminated by published crosswalks.

`[TODO]` Situate against Zhang et al. (2019) and general crosswalk automation
literature.

## 2. Methods

### 2.1 Data

Two mapping tasks were evaluated at the three-character level: ICD-9-CM to
ICD-10-CA (937 manually verified pairs over 345 ICD-9-CM codes) and ICD-9-CM
to ICDA-8 (331 pairs over 302 codes). Codes on the accompanying exclusion
lists were removed. Code labels were available for 354 ICD-9-CM, 2,038
ICD-10-CA, and 858 ICDA-8 codes. Diagnosis co-occurrence frequencies between
ICD-9-CM and each target system were derived from records coded under both.

ICD-9-CM codes map to more than one target frequently: 63% of codes have
multiple correct ICD-10-CA targets, averaging 2.7 per code and reaching 13.
The task is therefore set prediction rather than one-to-one translation.

### 2.2 Stage one: candidate retrieval

For each ICD-9-CM code, candidate targets were drawn from the union of the
top-K most similar target codes under each of three embedding models
(SapBERT, all-mpnet-base-v2, ClinicalBERT) and the top-N most frequently
co-occurring targets.

This replaced the previous rule, which retained candidates whose similarity
was within a fixed fraction of the maximum similarity for that code. Section
3.3 shows why.

### 2.3 Stage two: reranking

Each candidate was scored by a gradient-boosted decision tree model
(xgboost, depth 6, 200 rounds) over: per-model cosine similarity, within-code
rank, and relative similarity; reciprocal rank fusion across models and
inter-model disagreement; **mutual agreement**, computed by normalising
similarity within the target code as well as within the source code, so that a
target that ranks every source code highly is discounted; co-occurrence
frequency, rank and share; chapter-pair compatibility learned from training
folds; and lexical features over the code labels (IDF-weighted token overlap,
Jaccard, normalised edit distance).

Chapter compatibility was learned rather than taken from the existing
hand-written alignment table, which was found to be incomplete: it discards 51
correct pairs, including all mappings from ICD-9-CM chapter 3 ("endocrine,
nutritional, metabolic and immunity") to ICD-10 immunity codes D80-D89.
Correcting the table by hand against the validation data would constitute
fitting to the evaluation set.

Candidates were emitted when the predicted probability exceeded an absolute
threshold or fell within a fixed fraction of the best candidate for that code,
the second criterion accommodating codes with several valid targets. Every
code emits at least its highest-scoring candidate.

### 2.4 Evaluation

Five-fold cross-validation grouped by ICD-9-CM code. All fitting and tuning —
reranker, retrieval settings, emission thresholds, chapter encoding — occurred
inside training folds; test folds were used once, for prediction. The
rule-based baseline was evaluated identically, with its threshold, top-N and
flag rule selected on each fold's training codes.

Grouping matters: a split over pairs rather than codes would place some of a
code's correct targets in training and others in test.

## 3. Results

### 3.1 Held-out performance

| Track | Method | Precision | Recall | F1 |
|---|---|---|---|---|
| ICD-9 → ICD-10-CA | Rule-based, in-sample | 0.692 | 0.452 | 0.547 |
| ICD-9 → ICD-10-CA | Rule-based, held-out | 0.691 | 0.451 | 0.546 |
| ICD-9 → ICD-10-CA | **Two-stage, held-out** | **0.742** | **0.607** | **0.668** |
| ICD-9 → ICDA-8 | Rule-based, in-sample | 0.861 | 0.790 | 0.824 |
| ICD-9 → ICDA-8 | Rule-based, held-out | 0.861 | 0.790 | 0.824 |
| ICD-9 → ICDA-8 | **Two-stage, held-out** | **0.867** | **0.815** | **0.840** |

The gain is driven by recall (0.451 to 0.607 on ICD-10-CA) at comparable
precision. The baseline's in-sample and held-out results are nearly identical,
indicating it was limited by its candidate retrieval rather than overfitted by
parameter selection.

A correct target was ranked first for **93.3%** (ICD-10-CA) and **88.4%**
(ICDA-8) of codes, and appeared within the top three for 97.7% and 93.7%.

### 3.2 Precision at coverage

| Precision target | ICD-10-CA reproduced | ICDA-8 reproduced |
|---|---|---|
| 80% | 60.6% | 85.6% |
| 90% | 49.4% | 82.1% |
| **95%** | **43.2%** | **79.0%** |
| 99% | 21.3% (100% precision achieved) | 63.9% |

At the 95% operating point, 86.1% and 84.1% of codes receive at least one
automatically accepted mapping; 13.6% and 11.9% are routed for review; and
0.3% and 4.0% have no plausible candidate.

### 3.3 What limits accuracy

Decomposing the original pipeline:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| Correct targets present in similarity matrix | 100% | 100% |
| Surviving the similarity threshold | 37.1% | 79.5% |
| Entering the candidate pool | 62.6% | 85.2% |
| **Oracle F1 (perfect reranking of that pool)** | **0.770** | **0.920** |

Every correct answer was retrievable; the retrieval rule discarded most of
them. Because a candidate not retrieved cannot be recovered downstream, no
embedding model or reranker could have exceeded F1 0.770 under the original
design. Replacing the relative threshold with top-K retrieval raised the pool's
recall ceiling to 0.949 and 0.994.

### 3.4 Component contributions

Change in held-out F1 when each feature group is removed:

| Removed | ICD-10-CA | ICDA-8 |
|---|---|---|
| Mutual (forward + reverse) agreement | **-0.063** | -0.011 |
| Chapter compatibility | -0.051 | -0.003 |
| Co-occurrence | -0.051 | **-0.036** |
| Lexical | -0.015 | -0.001 |
| ClinicalBERT | -0.012 | **+0.003** |

Mutual agreement is the largest single contributor on ICD-10-CA. ClinicalBERT
and the lexical features contribute little, and removing ClinicalBERT slightly
improves ICDA-8, suggesting the ensemble could be reduced to two models.

Under a full parameter sweep, the general-purpose all-mpnet-base-v2 matched
the domain-specific SapBERT on ICD-10-CA (F1 0.533 vs 0.534) but not on ICDA-8
(0.769 vs 0.821); both exceeded ClinicalBERT on both tracks. Domain-specific
pretraining alone did not determine performance; SapBERT's training objective,
biomedical entity name alignment, corresponds directly to this task.

### 3.5 Robustness

**Unseen clinical areas.** With folds grouped by CCS category rather than by
code, F1 changed by -0.004 (ICD-10-CA) and -0.000 (ICDA-8), both within
run-to-run variation. Performance does not depend on clinically adjacent codes
appearing in training.

**Training set size.** F1 plateaued at approximately 100 to 180 training
codes on both tracks, with a slightly negative slope over the largest sizes
tested. Additional labelled training pairs would not be expected to raise
accuracy; the limitation is the discriminative signal available, not the
quantity of labels.

## 4. Discussion

The principal finding is diagnostic. Accuracy in this pipeline was limited not
by the embedding model, which is where comparative effort had previously gone,
but by a candidate retrieval rule that discarded 63% of correct answers before
ranking began. This is measurable in advance of any modelling work, via the
oracle ceiling, and we suggest reporting it as standard practice for staged
retrieval-and-ranking systems.

Relative to LLM-based translation, this approach attaches calibrated
confidence to each mapping, allowing a defined fraction of the crosswalk to be
accepted without review and the remainder to be triaged. It also requires no
external service and no transfer of health data. Direct numerical comparison
with Monchka et al. is not appropriate: the translation direction is reversed,
the reference standard differs (CIHI versus the local manual crosswalk), the
code subsets differ, and the CIHI crosswalk is one-to-one where this task
averages more than one correct target per code.

`[TODO]` **Priority: re-run against the CIHI reference crosswalk used by
Monchka et al. (1,272 chronic disease codes).** A shared reference standard
would make the comparison direct and roughly triple the validation set.

### Limitations

Evaluation covers 354 three-character ICD-9-CM codes from a single
institution's data, at category level rather than full code specificity. The
manual crosswalks are the reference standard and are themselves imperfect.
Only two target systems were tested, both against ICD-9-CM. The retrieve-and-
rerank architecture is standard in information retrieval; the contribution
here is its application, the diagnostic analysis, and the evaluation protocol,
not the architecture.

`[TODO]` External validation on an independent crosswalk.
`[TODO]` Coverage is capped by co-occurrence data, which currently spans 347
ICD-9-CM codes. Extension requires co-occurrence counts over a wider code set.

---

## Reproducibility

All results are produced by the scripts in `scripts/`, numbered in order.
`09` and `10` contain the error decomposition; `11` and `12` build and
evaluate the system; `13` produces the coverage analysis; `16` through `20`
contain the ablation, learning curve, category holdout, and top-1 analyses.
