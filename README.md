# icd crosswalk automation: sapbert vs clinicalbert

r pipeline that maps icd-9-cm codes to icd-10-ca and icda-8, comparing
clinicalbert vs sapbert as the embedding model used for the similarity
signal. combines cosine similarity + co-occurrence frequency + chapter
filtering, validated against a manual crosswalk for precision/recall/f1/accuracy.

## packages needed

r: `dplyr`, `readxl`, `writexl`, `tidyr`, `ggplot2`, `reticulate`, `stringr`, `purrr`, `knitr`, `rmarkdown`

python (only if regenerating sapbert embeddings): `torch`, `transformers`

## how to run

run from inside `scripts/` (relative paths assume this), in order:

```r
source("run_sanity_check.R")   # confirms the r port matches published numbers, run first
source("02_run_comparison.R")  # runs the full parameter grid, writes results/*.csv
source("03_visualize_results.R")  # builds charts from results/*.csv
```

`01_generate_sapbert_embeddings.R` is optional, only run it if you need to
regenerate `data/sapbert/*.xlsx` from scratch. the precomputed matrices are
already there.

then knit the report from the project root (not `scripts/`):

```r
rmarkdown::render("report.Rmd")
```

## folder structure

```
data/original/       clinicalbert similarity matrices, co-occurrence tables, manual crosswalks
data/sapbert/         sapbert similarity matrices + embedding row order
scripts/               pipeline code (see above)
results/               csvs + charts written by 02/03
report.Rmd             the report
```
