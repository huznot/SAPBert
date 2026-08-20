# Automated ICD Crosswalk Construction Using Embedding Similarity, Diagnosis Co-occurrence, and Learned Reranking

Draft for discussion. All numbers are reproducible from this repository. Open
items are marked TODO.

---

## Abstract

**Background.** Crosswalk tables mapping diagnosis codes between ICD revisions
are expensive to build by hand and are required whenever a health system
adopts a new coding standard. Prior work evaluated GPT-4 on this task and
reported 84.7% accuracy at the three-character level for ICD-10-CA to ICD-9-CM
translation, concluding that performance was insufficient for deployment.

**Objective.** To evaluate whether a task-specific pipeline combining code
label embeddings, diagnosis co-occurrence, and a learned reranker can build
ICD crosswalks accurately enough for practical use, and to identify what
limits its accuracy.

**Methods.** ICD-9-CM codes were mapped to ICD-10-CA and ICDA-8 at the
three-character level. Candidate targets were retrieved using cosine
similarity from three embedding models and top-N diagnosis co-occurrence, then
reranked by a gradient-boosted model. Features included per-model similarity
and rank, rank fusion across models, mutual forward and reverse agreement,
co-occurrence, learned chapter compatibility, and lexical overlap.
Performance was estimated by five-fold cross-validation grouped by ICD-9-CM
code, with all fitting and tuning confined to training folds. The existing
rule-based pipeline was evaluated under the same protocol.

**Results.** Held-out F1 was 0.668 for ICD-9-CM to ICD-10-CA and 0.840 for
ICD-9-CM to ICDA-8, compared with 0.546 and 0.824 for the rule-based baseline.
A correct target ranked first for 93.3% and 88.4% of codes. At 95% precision,
37.5% and 76.1% of the manual crosswalk pairs were reproduced without review,
covering 86.1% and 84.1% of source codes.
Error decomposition showed candidate retrieval was the binding constraint:
all correct targets were present in the similarity matrices, but the original
retrieval rule admitted 62.6% of them, capping attainable F1 at 0.770.
Performance was stable when whole clinical categories were held out (-0.004,
-0.000) and plateaued at roughly 150 training codes.

**Conclusions.** The pipeline reproduces most of a manually built ICD
crosswalk and attaches calibrated confidence to each mapping, so
high-precision mappings can be accepted automatically and the rest routed for
review. Candidate retrieval limits accuracy more than label quantity or
embedding choice.

---

## 1. Background

ICD coding standards are revised periodically. Each revision strands records
coded under the previous standard, and bridging them requires a crosswalk
mapping every legacy code to a modern equivalent. Crosswalks are built
manually by clinical coders, are expensive, and must be rebuilt for each pair
of systems and each national variant.

Prior work in this group evaluated GPT-4 across nine prompting strategies for
ICD-10-CA to ICD-9-CM translation, reporting 48.3% full-code and 84.7%
three-character accuracy, and concluded performance was insufficient for
deployment. That study noted that published crosswalks may appear in the
model's training data, so some accuracy may reflect memorisation rather than
inference.

The approach evaluated here uses two signals any health system undertaking a
migration already has: the text of the code labels, and how often codes
co-occur in records coded under both systems. Neither requires an external
service or transferring data to a third party. Co-occurrence counts derived
from local records cannot be contaminated by published crosswalks.

TODO: position against Zhang et al. (2019) and the wider crosswalk automation
literature.

## 2. Methods

### 2.1 Data

Two tasks were evaluated at the three-character level: ICD-9-CM to ICD-10-CA
(937 manually verified pairs over 345 ICD-9-CM codes) and ICD-9-CM to ICDA-8
(331 pairs over 302 codes). Codes on the accompanying exclusion lists were
removed. Labels were available for 354 ICD-9-CM, 2,038 ICD-10-CA, and 858
ICDA-8 codes. Co-occurrence frequencies between ICD-9-CM and each target
system were derived from records coded under both.

Most ICD-9-CM codes map to several targets. 63% have more than one correct
ICD-10-CA target, averaging 2.7 per code and reaching 13. The task is
therefore set prediction rather than one-to-one translation.

### 2.2 Stage one: candidate retrieval

For each ICD-9-CM code, candidates were taken from the union of the top-K most
similar targets under each of three embedding models (SapBERT,
all-mpnet-base-v2, ClinicalBERT) and the top-N most frequently co-occurring
targets.

This replaced the previous rule, which kept candidates whose similarity fell
within a fixed fraction of the maximum similarity for that code. Section 3.3
gives the reason.

### 2.3 Stage two: reranking

Candidates were scored by a gradient-boosted decision tree model (xgboost,
depth 6, 200 rounds). Features covered per-model cosine similarity,
within-code rank and relative similarity; reciprocal rank fusion across models
and inter-model disagreement; mutual agreement, computed by normalising
similarity within the target code as well as the source code so that targets
ranking highly for many sources are discounted; co-occurrence frequency, rank
and share; chapter-pair compatibility learned from training folds; and lexical
overlap between labels (IDF-weighted token overlap, Jaccard, normalised edit
distance).

Chapter compatibility was learned rather than taken from the existing
hand-written alignment table. That table discards 51 correct pairs, including
all mappings from ICD-9-CM chapter 3 (endocrine, nutritional, metabolic and
immunity) to ICD-10 immunity codes D80-D89. Correcting it by hand against the
validation data would amount to fitting the evaluation set.

Candidates were emitted when the predicted probability exceeded an absolute
threshold, or fell within a fixed fraction of the best candidate for that code.
The second criterion accommodates codes with several valid targets. Every code
emits at least its highest-scoring candidate.

### 2.4 Evaluation

Five-fold cross-validation grouped by ICD-9-CM code. The reranker, retrieval
settings, emission thresholds and chapter encoding were all fitted inside
training folds, and test folds were used once for prediction. The rule-based
baseline was evaluated identically, with its threshold, top-N and flag rule
selected on each fold's training codes.

Grouping by code matters here. A split over pairs would place some of a code's
correct targets in training and the rest in test.

## 3. Results

### 3.1 Held-out performance

| Track | Method | Precision | Recall | F1 |
|---|---|---|---|---|
| ICD-9 to ICD-10-CA | Rule-based, in-sample | 0.692 | 0.452 | 0.547 |
| ICD-9 to ICD-10-CA | Rule-based, held-out | 0.691 | 0.451 | 0.546 |
| ICD-9 to ICD-10-CA | Two-stage, held-out | 0.742 | 0.607 | 0.668 |
| ICD-9 to ICDA-8 | Rule-based, in-sample | 0.861 | 0.790 | 0.824 |
| ICD-9 to ICDA-8 | Rule-based, held-out | 0.861 | 0.790 | 0.824 |
| ICD-9 to ICDA-8 | Two-stage, held-out | 0.867 | 0.815 | 0.840 |

The improvement comes from recall (0.451 to 0.607 on ICD-10-CA) at comparable
precision. The baseline's in-sample and held-out results are nearly identical,
which indicates it was limited by candidate retrieval rather than overfitted
through parameter selection.

A correct target ranked first for 93.3% (ICD-10-CA) and 88.4% (ICDA-8) of
codes, and appeared in the top three for 97.7% and 93.7%.

### 3.2 Precision at coverage

Two denominators matter and they answer different questions. Pair-level
recall asks how much of the crosswalk is complete. Code coverage asks how
many source codes receive a confident mapping at all.

| Precision target | ICD-10-CA pairs | ICD-10-CA codes | ICDA-8 pairs | ICDA-8 codes |
|---|---|---|---|---|
| 80% | 52.6% | 96.8% | 82.5% | 93.7% |
| 90% | 42.9% | 92.2% | 79.2% | 90.1% |
| 95% | 37.5% | 86.1% | 76.1% | 84.1% |
| 99% | 18.5% | 49.6% | 61.6% | 67.9% |

Percentages are of the 937 and 331 manually verified pairs, and of the 345
and 302 source codes.

The gap between the two tracks is pair-level, not code-level. ICD-10-CA codes
have 2.7 correct targets on average against 1.09 for ICDA-8, and at 95%
precision the system emits about 1.07 mappings per code. It recovers the
primary mapping for most codes and misses secondary targets. At the 95%
operating point 13.6% and 11.9% of codes are routed for review and 0.3% and
4.0% have no plausible candidate.

### 3.3 What limits accuracy

Decomposing the original pipeline:

| | ICD-10-CA | ICDA-8 |
|---|---|---|
| Correct targets present in similarity matrix | 100% | 100% |
| Surviving the similarity threshold | 37.1% | 79.5% |
| Entering the candidate pool | 62.6% | 85.2% |
| Oracle F1, perfect reranking of that pool | 0.770 | 0.920 |

Every correct answer was retrievable, and the retrieval rule discarded most of
them. A candidate that is not retrieved cannot be recovered downstream, so no
embedding model or reranker could exceed F1 0.770 under the original design.
Replacing the relative threshold with top-K retrieval raised the pool's recall
ceiling to 0.949 and 0.994.

### 3.4 Component contributions

Change in held-out F1 when each feature group is removed:

| Removed | ICD-10-CA | ICDA-8 |
|---|---|---|
| Mutual forward and reverse agreement | -0.063 | -0.011 |
| Chapter compatibility | -0.051 | -0.003 |
| Co-occurrence | -0.051 | -0.036 |
| Lexical | -0.015 | -0.001 |
| ClinicalBERT | -0.012 | +0.003 |

Mutual agreement contributes most on ICD-10-CA. ClinicalBERT and the lexical
features contribute little, and removing ClinicalBERT slightly improves
ICDA-8, so the ensemble could likely be reduced to two models.

Under a full parameter sweep the general-purpose all-mpnet-base-v2 matched the
domain-specific SapBERT on ICD-10-CA (F1 0.533 against 0.534) but not on
ICDA-8 (0.769 against 0.821). Both exceeded ClinicalBERT on both tracks.
Domain-specific pretraining alone did not determine performance. SapBERT's
training objective, biomedical entity name alignment, corresponds closely to
this task.

### 3.5 Robustness

Holding out whole CCS categories rather than random codes changed F1 by -0.004
(ICD-10-CA) and -0.000 (ICDA-8), both within run-to-run variation. Performance
does not depend on clinically adjacent codes appearing in training.

F1 plateaued at roughly 100 to 180 training codes on both tracks, with a
slightly negative slope over the largest sizes tested. Additional labelled
training pairs would not be expected to raise accuracy. The limitation is the
discriminative signal available rather than the quantity of labels.

## 4. Discussion

The main finding is diagnostic. Accuracy was limited by a candidate retrieval
rule that discarded 63% of correct answers before ranking began, rather than
by the embedding model, which is where comparative effort had previously gone.
The oracle ceiling makes this measurable before any modelling work, and we
suggest reporting it as standard practice for staged retrieval and ranking
systems.

Compared with LLM-based translation, this approach attaches calibrated
confidence to each mapping, so a defined fraction of the crosswalk can be
accepted without review and the rest triaged. It requires no external service
and no transfer of health data.

Direct numerical comparison with the GPT-4 results is not appropriate. The
translation direction is reversed, the reference standard differs (CIHI
against the local manual crosswalk), the code subsets differ, and the CIHI
crosswalk is one-to-one where this task averages more than one correct target
per code.

TODO, priority: re-run against the CIHI reference crosswalk used in the GPT-4
study (1,272 chronic disease codes). A shared reference standard would make
the comparison direct and roughly triple the validation set.

### Limitations

Evaluation covers 354 three-character ICD-9-CM codes from one institution's
data, at category level rather than full code specificity. The manual
crosswalks serve as the reference standard and are themselves imperfect. Only
two target systems were tested, both against ICD-9-CM. Retrieve-and-rerank is
a standard information retrieval architecture; the contribution here is the
diagnostic analysis, the evaluation protocol, and the application, not the
architecture.

TODO: external validation on an independent crosswalk.

TODO: coverage is capped by the co-occurrence data, which spans 347 ICD-9-CM
codes. Extending it requires co-occurrence counts over a wider code set.

---

## Reproducibility

All results come from the numbered scripts in `scripts/`. 09 and 10 contain
the error decomposition, 11 and 12 build and evaluate the system, 13 produces
the coverage analysis, and 16 through 20 cover the ablation, learning curve,
category holdout and top-1 analyses.
