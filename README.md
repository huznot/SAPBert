![Logos](results/figures/logos.png)

# icd crosswalk

maps old medical diagnosis codes to newer ones automatically.

when a country switches coding systems, decades of health records are stuck in
the old system. the crosswalk from every old code to its modern equivalent is
normally built by hand, by clinical coders, over months. this builds on an
existing pipeline and tests what changes to it are worth making.

two migrations:

- icd-9-cm to icd-10-ca, the current canadian standard
- icd-9-cm to icda-8, going backwards

## how it works

for each old code the pipeline builds two candidate lists, then picks from them:

1. cosine similarity between the code labels, from a sentence embedding model
2. the codes it most often co-occurs with in administrative records
3. a rule that combines the two lists
4. a chapter filter, so a mapping cannot cross clinical areas

## results

| | |
|---|---|
| [clinicalbert vs sapbert](results/head_to_head/clinicalbert_vs_sapbert.xlsx) | both models on the same settings, overall and by category |
| [counts at each cutoff](results/workbooks/1_threshold_frequencies.xlsx) | how many pairs each cutoff keeps, with true and false positives |
| [code number in the label](results/workbooks/2_icda8_code_in_label.xlsx) | whether embedding the code number with the label matters |
| [why sapbert scores lower](results/workbooks/3_similarity_gap_by_model.xlsx) | correct pairs against everything else, per model |
| [all 130 ccs categories](results/workbooks/4_counts_by_ccs_category.xlsx) | counts per category |
| [the 52 with no icda-8 match](results/review/52_codes_review.xlsx) | the review of those codes |

best f1 is 0.544 on icd-9 to icd-10-ca and 0.761 on icd-9 to icda-8, against
0.423 and 0.716 for the original pipeline. those score all 354 icd-9-cm codes,
including the ones with no match in the target system.

## running it

needs r. all the data is in the repo, so nothing has to be downloaded. scripts
find the repo root themselves, so it does not matter where you start them from.

```bash
Rscript scripts/27_show_results.R   # every headline number, about a second
Rscript run_all.R                   # rebuild everything, about 45 minutes
Rscript run_all.R --quick           # ~10 minutes, skips the parameter searches
```

r packages: `dplyr`, `tidyr`, `readxl`, `stringr`, `purrr`, `ggplot2`,
`jsonlite`, `xgboost`, `stopwords`, `openxlsx`.

python is only needed to rebuild the embeddings:

```bash
python scripts/generate_embeddings.py --model mpnet --clean base
```

## layout

```
data/        labels, co-occurrence tables, manual crosswalks, matrices
scripts/     pipeline code, numbered in the order it was written
docs/        the paper drafts and the change log
results/     everything the scripts produce
run_all.R    rebuilds it all in dependency order
```

scripts call `out_path("name.csv")` and `scripts/paths.R` works out which
results folder it belongs in from the file name, so the layout and the code
cannot drift apart.

`data/` is committed, which is why the repo is about 350 mb. that keeps the
analysis reproducible from a clone without a gpu.

## branches

`main` follows the existing four step methodology. `testing` replaces the
selection step with a wide candidate set and a trained scoring model, reaching
0.668 and 0.840 on unseen codes. it is kept separate because it departs from the
agreed methodology.
