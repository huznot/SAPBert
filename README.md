# icd crosswalk automation: sapbert vs clinicalbert

r pipeline that maps icd-9-cm codes to icd-10-ca and icda-8, comparing
clinicalbert vs sapbert as the embedding model used for the similarity
signal. combines cosine similarity + co-occurrence frequency + chapter
filtering, validated against a manual crosswalk for precision/recall/f1/accuracy.

## packages needed

r: `dplyr`, `readxl`, `writexl`, `tidyr`, `ggplot2`, `reticulate`, `stringr`, `purrr`, `knitr`, `rmarkdown`, `jsonlite`

python (only for regenerating embeddings -- sapbert, filler-stripped
variants, or the mpnet arm): `torch`, `transformers`, `sentence-transformers`,
`pandas`, `openpyxl`. Plain Python (no reticulate) -- `01_generate_sapbert_embeddings.R`
uses reticulate for the original sapbert matrices, but `generate_embeddings.py`
(used for everything since) calls the models directly in Python, since
report.Rmd documents that reticulate silently returns wrong embeddings for
this.

## how to run

run from inside `scripts/` (relative paths assume this), in order:

```r
source("run_sanity_check.R")          # confirms the r port matches published numbers, run first
source("02_run_comparison.R")         # runs the full parameter grid, writes results/*.csv
source("03_visualize_results.R")      # builds charts from results/*.csv
source("05_bidirectional_and_roundtrip.R")  # reverse-direction mapping + round-trip consistency
source("06_extended_comparison.R")    # filler-word stripping + general-model (mpnet) arms
```

`01_generate_sapbert_embeddings.R` is optional, only run it if you need to
regenerate `data/sapbert/*.xlsx` from scratch. the precomputed matrices are
already there.

`06_extended_comparison.R` needs `data/generated/*.xlsx`, produced by:

```bash
python generate_embeddings.py --model clinicalbert --clean base
python generate_embeddings.py --model clinicalbert --clean stripped
python generate_embeddings.py --model sapbert       --clean stripped
python generate_embeddings.py --model mpnet          --clean base
python generate_embeddings.py --model mpnet          --clean stripped
```

(SapBERT/base reuses the existing `data/sapbert` matrices directly, so it's
not regenerated.) Filler words stripped before embedding are configured in
`filler_words.json`, shared by the R and Python code.

then knit the report from the project root (not `scripts/`):

```r
rmarkdown::render("report.Rmd")
```

## folder structure

```
data/original/       clinicalbert similarity matrices, co-occurrence tables, manual crosswalks
data/sapbert/         sapbert similarity matrices + embedding row order
data/generated/        filler-stripped / mpnet similarity matrices (scripts/generate_embeddings.py)
scripts/               pipeline code (see above)
results/               csvs + charts written by 02/03/05/06
report.Rmd             the report
PI_QUESTIONS.md         open questions for Dr. Lix (ambiguous PI direction, missing data)
```
