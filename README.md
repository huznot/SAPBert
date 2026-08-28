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

`main` follows the existing four-step methodology. The changes on it are a
newer embedding model, corrections to the parameter search, an analysis of
where correct mappings are lost, and the text-cleaning work requested in
review. Best F1 is 0.524 and 0.761, against 0.423 and 0.716 for the original.

Those figures score all 354 ICD-9-CM codes, including the 9 and 52 that have no
match in the target system. The original scored only the codes that had one.
`report.md` Section 2 covers what that changed.

`testing` additionally replaces the selection step with a wide candidate set
and a trained scoring model, reaching 0.668 and 0.840 on unseen codes. It is
kept separate because it departs from the existing methodology and has not been
agreed.

## Running it

Needs R. All the data is in the repo, so nothing has to be downloaded or
requested. Scripts find the repo root themselves, so it does not matter which
directory you start them from.

Every headline number, read from the committed results, about a second:

```bash
Rscript scripts/27_show_results.R
```

To rebuild everything from the data instead, in dependency order:

```bash
Rscript run_all.R           # about 45 minutes
Rscript run_all.R --quick   # ~10 minutes, skips the parameter searches
```

Or run a single stage:

```bash
Rscript scripts/07_full_grid_comparison.R   # parameter grid, all conditions
Rscript scripts/08_assemble_full_grid.R     # combine into summary tables
Rscript scripts/09_error_analysis.R         # where correct mappings are lost
Rscript scripts/23_code_prefix_test.R       # code number in the embedded text
```

R packages: `dplyr`, `tidyr`, `readxl`, `writexl`, `stringr`, `purrr`,
`ggplot2`, `xgboost`, `jsonlite`, `reticulate`.

Python is only needed to regenerate embeddings (`torch`, `transformers`,
`sentence-transformers`, `pandas`, `openpyxl`):

```bash
python scripts/generate_embeddings.py --model mpnet --clean base
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
| `paths.R` | finds the repo root so scripts run from any directory |

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
| `16_ablation.R` | which feature groups matter in the scoring model |
| `17_retrieval_sensitivity.R` | how much the candidate settings matter |
| `18_learning_curve.R` | would more training data help |
| `19_category_holdout.R` | does it work on unseen clinical areas |
| `20_top1_accuracy.R` | is the right answer reachable at all |
| `21_error_by_code_type.R` | single vs multi target codes |
| `22_target_block_structure.R` | do multi-target answers sit in blocks |
| `23_code_prefix_test.R` | code number in the text, report Section 5 |
| `29_portability.R` | performance without health records or a chapter table |
| `30_variability.R` | bootstrap standard deviation on every reported score and difference, report Section 9 |
| `31_category_breakdown.R` | performance across all 130 CCS categories, every condition, report Section 10 |

**Show the results, no recomputation, a second or two each**

| | |
|---|---|
| `27_show_results.R` | every headline number |
| `34_export_rds_as_csv.R` | the held-out predictions as csv, for reading without R |
| `35_unmatched_codes.R` | codes with no correct answer, report Section 2 and Figure 1 |
| `36_unmatched_descriptives.R` | similarity and co-occurrence for those codes on their own, report Section 11 |
| `37_unmatched_handout.R` | the same thing as a standalone html page, for sending to someone without the repo |
| `24_show_similarity_matrix.R` | a worked similarity matrix |
| `25_frequency_distributions.R` | report Section 8 |
| `26_stopword_choice.R` | which stop word dictionary, report Section 6 |
| `28_stopwords_and_codes.R` | report Section 7 |

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
results/csv_export/  the held-out predictions as csv
run_all.R         rebuilds everything in dependency order
report.md         every change from the original pipeline to now
```
