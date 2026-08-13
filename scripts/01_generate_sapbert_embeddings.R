# ============================================================================
# 01_generate_sapbert_embeddings.R
#
# Generates SapBERT embeddings for the ICD-9-CM, ICD-10-CA, and ICDA-8
# code+label text (same universe used for the original ClinicalBERT
# embeddings), using the `reticulate` package to call Python's
# `transformers`/`torch` under the hood -- there is no native R
# implementation of transformer models with pretrained weights, so this is
# the standard R approach for this kind of task. The rest of the pipeline
# is pure R.
#
# NOTE: this step downloads the SapBERT model (~440MB) the first time it
# runs, and takes a few minutes on CPU. If you already have the pre-computed
# embeddings in data/sapbert/, you can skip this script entirely and go
# straight to 02_build_similarity_matrices.R.
#
# Requires: pip install torch transformers  (in whatever Python reticulate
# is pointed at -- see py_config() below).
# ============================================================================

library(reticulate)
library(readxl)
library(dplyr)
library(stringr)

# Point reticulate at a Python that has torch + transformers installed.
# Adjust this path for your machine, or remove the use_python() call to let
# reticulate auto-discover a suitable environment.
use_python("/usr/bin/python3", required = TRUE)
py_config()

# IMPORTANT: the actual model forward pass is done via a plain Python
# function defined here and called from R -- NOT by calling the PyTorch
# model object directly using R's `model(...)` wrapper syntax. Testing
# showed that reticulate's translation of `model(input_ids = ..., ...)`
# into a PyTorch nn.Module `__call__` silently produces WRONG numerical
# output (different values, not just a different data type) compared to
# calling the exact same model the normal Python way -- while a plain
# user-defined Python function called via `py$fn(...)` gives byte-for-byte
# correct results. This was caught by comparing against a known-correct
# Python-generated embedding for the same code+label text; the same
# calling pattern is used throughout this script and its correctness was
# separately re-verified against the Python reference before use.
py_run_string("
import torch
from transformers import AutoTokenizer, AutoModel

_tokenizer = None
_model = None

def _load_sapbert():
    global _tokenizer, _model
    if _model is None:
        name = 'cambridgeltl/SapBERT-from-PubMedBERT-fulltext'
        _tokenizer = AutoTokenizer.from_pretrained(name)
        _model = AutoModel.from_pretrained(name)
        _model.eval()

def embed_batch(texts, max_length=64):
    _load_sapbert()
    toks = _tokenizer(list(texts), padding='max_length', max_length=max_length,
                       truncation=True, return_tensors='pt')
    with torch.no_grad():
        out = _model(**toks)
        cls = out[0][:, 0, :]  # CLS token of the last hidden layer
    return cls.numpy()
")

BASE <- "../data/original"
OUT  <- "../data/sapbert"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

clean_label <- function(s) {
  s <- tolower(as.character(s))
  s <- str_replace_all(s, "[^a-z0-9\\s]", " ")
  s <- str_squish(s)
  s
}

cat("Loading code/label files...\n")
icd8  <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"), sheet = "ICDA-8-3Level")
icd9  <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"), sheet = "CCS ICD-9-CM-3Level")
icd10 <- read_excel(file.path(BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"), sheet = "ICD-10-CA-3Level")

icd8  <- icd8  %>% mutate(text = paste(ICDA_8, clean_label(ICDA_8_LABEL)))
icd9  <- icd9  %>% mutate(text = paste(ICD_9_CM, clean_label(ICD_9_CM_LABEL)))
icd10 <- icd10 %>% mutate(text = paste(ICD_10_CA, clean_label(ICD_10_CA_LABEL)))

cat(sprintf("icd8: %d, icd9: %d, icd10: %d\n", nrow(icd8), nrow(icd9), nrow(icd10)))

cat("Loading SapBERT model (downloads ~440MB the first time)...\n")

embed_texts <- function(texts, batch_size = 32L, max_length = 64L) {
  all_embs <- list()
  n <- length(texts)
  i <- 1L
  while (i <= n) {
    j <- min(i + batch_size - 1L, n)
    batch <- as.list(texts[i:j])
    emb <- py$embed_batch(batch, as.integer(max_length))  # see py_run_string() block above
    all_embs[[length(all_embs) + 1]] <- emb
    cat(sprintf("  %d/%d\n", j, n))
    i <- j + 1L
  }
  do.call(rbind, all_embs)
}

cat("Embedding ICD-9-CM...\n")
emb9 <- embed_texts(icd9$text)
cat("Embedding ICD-10-CA...\n")
emb10 <- embed_texts(icd10$text)
cat("Embedding ICDA-8...\n")
emb8 <- embed_texts(icd8$text)

saveRDS(emb9,  file.path(OUT, "icd9_embeddings_SapBERT.rds"))
saveRDS(emb10, file.path(OUT, "icd10_embeddings_SapBERT.rds"))
saveRDS(emb8,  file.path(OUT, "icd8_embeddings_SapBERT.rds"))

write.csv(icd9  %>% select(ICD_9_CM),  file.path(OUT, "icd9_order.csv"),  row.names = FALSE)
write.csv(icd10 %>% select(ICD_10_CA), file.path(OUT, "icd10_order.csv"), row.names = FALSE)
write.csv(icd8  %>% select(ICDA_8),    file.path(OUT, "icd8_order.csv"),  row.names = FALSE)

cat(sprintf("Done. Shapes: icd9 %s, icd10 %s, icd8 %s\n",
            paste(dim(emb9), collapse = "x"),
            paste(dim(emb10), collapse = "x"),
            paste(dim(emb8), collapse = "x")))
