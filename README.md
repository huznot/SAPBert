# ICD Crosswalk Automation

Maps old medical diagnosis codes to newer ones automatically.

When a country switches coding systems, decades of health records are stuck in
the old system. Someone has to build a crosswalk, a mapping from every old
code to its modern equivalent and that is normally done by hand, by clinical
coders, over months.

This work builds on an existing pipeline and reports what changes to it are
worth making. It is tested on two migrations:

- **ICD-9-CM → ICD-10-CA** (Canadian, the modern standard)
- **ICD-9-CM → ICDA-8** (historical, going backwards)

![Logos](results/logos.png)

## Branches

This is the `testing` branch. It keeps the two-stage pipeline: a wide candidate
set followed by a trained scoring model, reaching F1 0.668 and 0.840 on unseen
codes against 0.546 and 0.824 for the original rules measured the same way. The
ICD-10-CA gain is well outside the measurement noise; the ICDA-8 gain is not,
see `report.md` Section 12.

`main` follows the existing four-step methodology and does not include the
scoring model. Work here should not be merged into `main` until the change of
methodology has been agreed.

## Running it

Needs R. Open `icd_crosswalk.Rproj` in RStudio, then:

```r
source("scripts/27_show_results.R")   # every headline number, about a second
```

To reproduce the analysis from the data, from `scripts/`:

```r
source("07_full_grid_comparison.R")  # parameter grid, all conditions
source("08_assemble_full_grid.R")    # combine into summary tables
source("09_error_analysis.R")        # where correct mappings are lost
source("23_code_prefix_test.R")      # code number in the embedded text
```

R packages: `dplyr`, `tidyr`, `readxl`, `writexl`, `stringr`, `purrr`,
`ggplot2`, `xgboost`, `jsonlite`, `reticulate`.

Python is only needed to regenerate embeddings (`torch`, `transformers`,
`sentence-transformers`, `pandas`, `openpyxl`):

```bash
python generate_embeddings.py --model mpnet --clean base
```

### All scripts

Numbering is chronological, the order the work was done in. Everything below
still runs.

**Build the inputs**

| | |
|---|---|
| `generate_embeddings.py` | similarity matrices for any model and text cleaning |
| `01_generate_sapbert_embeddings.R` | how `data/sapbert/` was made, superseded by the Python script |
| `pipeline_lib.R` | shared functions, sourced by everything |

**The original pipeline, reproduced and searched**

| | |
|---|---|
| `02_run_comparison.R` | first ClinicalBERT vs SapBERT run |
| `03_visualize_results.R` | per-CCS-category charts |
| `05_bidirectional_and_roundtrip.R` | mapping in the reverse direction |
| `06_extended_comparison.R` | wider parameter sweep, superseded by `07` |
| `07_full_grid_comparison.R` | full grid, one file per condition |
| `08_assemble_full_grid.R` | combines those into the summary tables and report Figure 1 |

**Diagnosing the ceiling**

| | |
|---|---|
| `09_error_analysis.R` | where correct pairs are lost, report Section 4 |
| `10_candidate_generation_study.R` | what a wider candidate step would reach |

**The revised pipeline**

| | |
|---|---|
| `11_rerank_features.R` | candidates and the 52 features |
| `12_cv_rerank.R`, `12b_merge_cv_results.R` | train and evaluate, held out |
| `13_precision_coverage.R` | confidence thresholds and triage |
| `14_predict_crosswalk.R` | map codes with no known answer |

**Follow-up experiments**

| | |
|---|---|
| `16_ablation.R` | which feature groups matter, report Section 7 |
| `17_retrieval_sensitivity.R` | how much the candidate settings matter |
| `18_learning_curve.R` | would more training data help, report Figure 3 |
| `19_category_holdout.R` | does it work on unseen clinical areas |
| `20_top1_accuracy.R` | is the right answer reachable at all |
| `21_error_by_code_type.R` | single vs multi target codes |
| `22_target_block_structure.R` | do multi-target answers sit in blocks |
| `23_code_prefix_test.R` | code number in the text, report Section 8 |
| `29_portability.R` | performance without health records or a chapter table |
| `30_retrieval_vs_model.R` | separates the wider candidate set from the scoring model |
| `31_variability.R` | bootstrap standard deviation on the held-out scores and on the gap between the two systems, report Section 12 |
| `32_category_breakdown.R` | performance across all 130 CCS categories, both systems, report Section 13 |

**Show the results, no recomputation, a second or two each**

| | |
|---|---|
| `27_show_results.R` | every headline number |
| `24_show_similarity_matrix.R` | a worked similarity matrix |
| `25_frequency_distributions.R` | report Section 11 |
| `26_stopword_choice.R` | which stop word dictionary, report Section 9 |
| `28_stopwords_and_codes.R` | report Section 10 |

`04` and `15` were removed. `15` compared stop word dictionaries using a
hardcoded approximation of the word lists rather than the real ones and gave
wrong collision counts; `26` does the same job correctly.

## Using it on other code systems

To point this at a different migration you need, per code system:

- code labels (code → text description)
- co-occurrence counts between old and new codes, from records coded both ways
- a chapter/category grouping, if one exists
- some manually verified mappings for training, a few hundred is enough

Then swap the paths in the `TRACKS` list at the top of `11_rerank_features.R`.

## Layout

```
data/original/    labels, co-occurrence tables, manual crosswalks
data/sapbert/     SapBERT similarity matrices
data/generated/   regenerated / filler-stripped / mpnet matrices
scripts/          pipeline code
results/          csv output and charts
report.md         every change from the original pipeline to now
```
