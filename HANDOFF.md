# Handoff — ICD crosswalk automation

Written 2026-08-19. Read this first when resuming.

---

## 1. Where the project actually stands

The original 5-item task list is essentially done (details in §6). The work
has since moved to a much more valuable question: **why is accuracy capped,
and how high can it actually go?** That investigation produced a clear,
publishable answer, and a new two-stage system that is mid-evaluation.

### The headline finding (this is the paper)

The embedding model was **never** the bottleneck. Error analysis
(`scripts/09_error_analysis.R`) established:

| Fact | ICD-9→ICD-10-CA | ICD-9→ICDA-8 |
|---|---|---|
| True target codes present *somewhere* in the similarity matrix | **100%** | **100%** |
| ...that survive into the candidate pool | 62.6% | 85.2% |
| **Oracle F1** (a *perfect* reranker on that pool) | **0.770** | **0.920** |
| Actual F1 | 0.524 | 0.761 |

So the pipeline retrieves every correct answer and then **throws a third of
them away before ranking**. No model swap, ensemble, or reranker could ever
have exceeded 0.770 while candidate generation stayed as it was. That single
number reframes the whole project.

**Why they were thrown away:** the similarity filter is *relative to each
code's own maximum* (`keep sim >= 0.95 * max_sim`), which admits only ~37% of
true pairs. Replacing it with plain **top-K retrieval**
(`scripts/10_candidate_generation_study.R`) lifts the ceiling enormously:

| Track | Old ceiling | Top-K ceiling | Oracle F1 |
|---|---|---|---|
| ICD-10-CA | 0.626 | **0.927** | 0.77 → **0.96** |
| ICDA-8 | 0.852 | **0.976** | 0.92 → **0.99** |

Secondary finding: the hand-written chapter-alignment table is **incomplete**
and destroys 51 true pairs outright. Clearest example — ICD-9 chapter 3 is
"endocrine, nutritional, metabolic **AND IMMUNITY**" but is never permitted to
map to ICD-10 D80–D89 (immunity). It is *not* hand-patched, deliberately:
patching it against the validation answers would be leakage. Instead
chapter-pair compatibility is **learned inside training folds**
(`fit_chapter_encoding()` in `12_cv_rerank.R`).

### The new system

Two-stage retrieve-and-rerank:

- `scripts/11_rerank_features.R` — builds the candidate pool + features.
  Pool recall ceiling **0.949** (ICD-10-CA) / **0.994** (ICDA-8).
  Features: per-model similarity/rank/relative-score for SapBERT + mpnet +
  ClinicalBERT, **reciprocal-rank-fusion ensemble** across them, model
  disagreement (SD), co-occurrence frequency/rank/share, chapter features,
  and deliberately **non-neural lexical** features (IDF-weighted token
  overlap, Jaccard, first-token match, normalized edit distance).
- `scripts/12_cv_rerank.R` — xgboost reranker under **5-fold CV grouped by
  ICD-9 code**. Retrieval config (K, N, chapter on/off), decision rule
  (tau, rho) and chapter encoding are **all fitted inside training folds
  only**. The old pipeline is evaluated under the *identical* protocol so
  the comparison is held-out vs held-out.
- `scripts/13_precision_coverage.R` — deployment view (see §3). **Written
  but never yet run** — it needs `results/cv_rerank_predictions.rds` from a
  completed 12_ run.

### FINAL held-out results (both tracks, complete)

`12_cv_rerank.R` and `13_precision_coverage.R` have both **run to
completion**. These are the real numbers — nothing here was tuned on the
data it is scored on.

| Track | Baseline in-sample | Baseline held-out | **Reranker held-out** |
|---|---|---|---|
| ICD-9→ICD-10-CA | 0.547 | 0.546 | **0.614** |
| ICD-9→ICDA-8 | 0.824 | 0.819 | **0.833** |

Two things worth noting:

- **Baseline in-sample ≈ baseline held-out** (0.547 vs 0.546). So the old
  parameter sweep was *not* badly overfit — it was simply **capped**. The
  selection-on-test worry turned out to be minor; the retrieval ceiling was
  the real problem, exactly as §1 predicted.
- The reranker's gain is almost entirely **recall**: 0.451 → **0.567** on
  ICD-10-CA, at roughly unchanged precision (0.691 → 0.668). That is the
  signature of fixing candidate generation, not of a better classifier.

### The deployment result (this is the headline)

From `13_precision_coverage.R`, on out-of-fold predictions:

| Operating point | ICD-10-CA | ICDA-8 |
|---|---|---|
| **Crosswalk auto-built at ~95% precision** | **38.5%** | **76.6%** |
| Auto-built at ~99% precision | 20.0% | 57.5% |

Triage at the 95%-precision operating point:

| Bucket | ICD-10-CA | ICDA-8 |
|---|---|---|
| Auto-accept (no human needed) | 79.1% of codes | 85.1% |
| Review — candidates, low confidence | 20.3% | 13.9% |
| Review — no plausible candidate | 0.6% | 1.0% |

**"Three-quarters of the ICDA-8 crosswalk requires no human at all, at 95%
precision, on codes the system never saw" is the claim to lead with.**

Artifacts: `results/cv_rerank_results.csv`, `_folds.csv`, `_importance.csv`,
`_predictions.rds`, `precision_coverage_curve.csv`,
`_operating_points.csv`, `_triage.csv`, `plot_precision_coverage.png`.

---

## 2. IMMEDIATE NEXT STEP — write the report section

All the computation is done. **The gap is now writing, not running.**
`report.Rmd` currently documents only the full-grid model comparison and
says nothing about §1 at all.

Add a new section built around the **triage / precision-coverage table**,
with F1 as supporting detail rather than the headline. Suggested spine:

1. The diagnosis — 100% of answers retrievable, only 62.6% pooled, oracle
   F1 capped at 0.770. This is the intellectual core; lead with it.
2. The fix — relative threshold → top-K retrieval, ceiling 0.626 → 0.927.
3. The two-stage system and the grouped-CV protocol (stress that the
   baseline is evaluated under the *identical* protocol).
4. Held-out results table, then the precision-coverage/triage table.
5. Honest limits — §7 below.

Then re-render:

```bash
cd /c/Users/IRFANM/Downloads/icd_crosswalk_r
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e "Sys.setenv(RSTUDIO_PANDOC='C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'); rmarkdown::render('report.Rmd', output_format='html_document')"
```

To re-run the pipeline from scratch (only if inputs change — ~35 min):

```bash
cd scripts
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" 11_rerank_features.R
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" 12_cv_rerank.R > ../logs/12_cv_rerank.log 2>&1
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" 13_precision_coverage.R
```

---

## 3. The framing that makes this deployable

Dr. Lix's assistant made the key domain point: some ICD-9 codes **genuinely
have no counterpart** (terminology added/restructured between revisions), so
the goal is "map as much as we can", not "map everything".

That means **F1 is the wrong headline metric for deployment.** The useful
property is not being right about everything — it is *knowing which mappings
are reliable*, so a coder auto-accepts those and reviews only the rest.
`13_precision_coverage.R` produces exactly this:

- precision as a function of coverage
- the operating point achieving 95% / 99% precision, and how much of the
  crosswalk is automated at that precision
- a triage split: **auto-accept / review (low confidence) / review (no
  plausible candidate)**

**"At 95% precision we build X% of the crosswalk automatically" is a far
stronger symposium claim than any F1 number.** Prioritise this.

---

## 4. Ideas not yet tried (ranked by expected value)

1. **Reverse-direction / mutual-rank features.** Currently every feature is
   "how good is this target for this ICD-9 code?" Missing: "is this ICD-9
   code the best match *for that target*?" Mutual-nearest-neighbour is
   typically a large win in entity matching, and task 3 already showed
   reverse-direction info is informative. Cheap to add to `11_`. **Best
   next lever** — and now clearly the binding one: ICD-10-CA held-out recall
   is 0.567 against a pool ceiling of 0.949, so ~40% of retrievable answers
   are being retrieved and then discarded by the ranker.
2. **The inner CV keeps selecting K=10**, i.e. it retreats to a small pool
   because precision on a large pool is hard. Confirms the reranker, not
   retrieval, is the constraint. Stronger features (idea 1), more xgboost
   capacity/tuning, or a listwise/ranking objective instead of per-pair
   binary classification. Check `results/cv_rerank_importance.csv` first to
   see which features actually carry weight before adding more.
3. **Adaptive set size.** Truth averages 2.7 targets/code (range 1–13);
   predicting *how many* targets a code should get is a separate learnable
   problem from scoring each candidate.
4. **Per-CCS-category breakdown of the reranker**, to check gains are not
   concentrated in a few easy categories.
5. Filler-stripped embedding variants as *extra ensemble features* rather
   than as alternative arms.

---

## 5. Practical gotchas (all learned the hard way)

- **Rscript is not on PATH.** Use `"/c/Program Files/R/R-4.6.1/bin/Rscript.exe"`.
- **Never edit an .R file while an Rscript process is running it.** Rscript
  parses incrementally; a mid-run edit kills it with a bogus syntax error.
  This cost a full run already.
- **Windows filesystem is case-insensitive.** Two condition tags differing
  only by case = the same file; parallel jobs silently clobbered each
  other's results. `07_full_grid_comparison.R` now asserts tags are unique
  case-insensitively. Watch for this anywhere filenames come from identifiers.
- `rmarkdown::render()` needs pandoc:
  `Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")`.
- **`adist(a, b)` builds the full cross matrix** and segfaults on large
  vectors — compute pairwise on unique pairs (already fixed in `11_`).
- Bash heredocs in this environment mangle backslashes and apostrophes —
  use the Write tool for R scripts, not `cat <<EOF`.
- `xgboost` and `ranger` were installed from CRAN this session and work.
- **Never add a `Co-Authored-By: Claude` trailer to commits.** The user is
  emphatic. History was rewritten on 2026-08-19 to strip it from three
  commits; `main` is clean — keep it that way.
- Long jobs: run as background processes, each writing a distinct output
  file. 8 cores available.

---

## 6. Original task list status

1. **Filler-word stripping** — DONE. Full-grid, paired-point result: helps
   ClinicalBERT (112/112 points on ICD-10-CA, +0.017 F1), **hurts** SapBERT
   on ICDA-8 (111/112 points worse, −0.012). Model-dependent.
   **Recommendation: do not adopt globally.**
2. **General-purpose embedding model (mpnet)** — DONE, full grid. Ties
   SapBERT on ICD-10-CA (0.533 vs 0.530–0.534), clearly loses on ICDA-8
   (0.769 vs 0.821). Both beat ClinicalBERT everywhere. Conclusion: it is
   *not* domain-specificity that matters (ClinicalBERT is domain-specific
   and worst) but **training objective** — SapBERT trains on biomedical
   entity-name alignment, which is exactly this task.
3. **Bidirectional / round-trip** — DONE as a *reported metric only*.
   Whether it should actively *filter* mappings is still an open PI
   question (`PI_QUESTIONS.md` §0.5).
4. **Expand validation set** — BLOCKED, no additional data in the repo. PI
   deprioritised it.
5. **"Remove codes because cosine similarity drops"** — NOT IMPLEMENTED,
   underspecified. Questions listed in `PI_QUESTIONS.md` §1.

### Resolved: the baseline "discrepancy" was a false alarm
The reference numbers (0.427/0.271, 0.530/0.360, 0.716/0.557, 0.821/0.697)
**all reproduce exactly**. The problem was the top-N grid stopping at 15
when the true optimum is at **top_n = 30**. An earlier session wrongly
concluded "not reproducible" from a run that had been interrupted after 1 of
4 thresholds. Both tracks now sweep {3,5,10,15,20,25,30}. See `report.Rmd`.

---

## 7. Honest assessment (do not oversell this)

- **Near-perfect accuracy on ICD-10-CA is not achievable and should not be
  promised.** 63% of codes map to multiple targets, truth averages 2.7 per
  code, boundaries are genuine judgement calls, and human coders disagree.
  ICDA-8 is far more tractable (0.994 ceiling, mostly 1:1).
- The reranker recovers roughly a third of available headroom on ICD-10-CA
  (0.614 held-out against a ~0.97 oracle). Real work remains there. ICDA-8
  is in much better shape (0.833).
- **Do not quote 0.614 / 0.833 as if they were comparable to the old
  0.530 / 0.821.** The old figures are in-sample; the correct like-for-like
  comparison is against the held-out baseline (0.546 / 0.819), which is what
  `cv_rerank_results.csv` reports and what the report must use.
- The ICD-10-CA improvement is solid and real (+0.068 F1 held-out, driven by
  +0.116 recall). The ICDA-8 improvement is small (+0.014) — say so plainly
  rather than implying the method helps equally on both tracks.
- **The strongest contribution is the diagnosis, not the F1.** "We decomposed
  a published-style pipeline, showed its similarity filter discards 37% of
  recoverable answers and caps *any* downstream method at F1 0.77, fixed it,
  and evaluated held-out" is more defensible than another tenth of a point —
  and reviewers cannot wave it away.
- Every pre-`12_` number in this project is **selection-on-test** (tuned and
  scored on the same 937/331 pairs). They are fine for *comparing* arms, but
  they are not generalization estimates and must never be presented as such.
