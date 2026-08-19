# Builds the candidate pool and feature table for the reranker.
#
# Takes top-50 similarity per model plus top-50 co-occurrence, no chapter
# filter, and attaches all features. Narrower retrieval settings are just row
# filters on this (sim_rank <= K, cooc_rank <= N, chapter_ok), so 12_ can
# sweep them without rebuilding anything.
#
# The manual crosswalk is only joined at the end as the label y. Nothing here
# is fitted.
#
# Writes results/rerank_features_<track>.rds

source("pipeline_lib.R")

ORIG_BASE    <- "../data/original"
SAPBERT_BASE <- "../data/sapbert"
GEN_BASE     <- "../data/generated"
OUT_DIR      <- "../results"

MAX_SIM_K  <- 50   # top-K similarity candidates per ICD-9 code, per model
MAX_COOC_N <- 50   # top-N co-occurrence candidates per ICD-9 code

LABELS_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
VAL_XLSX    <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

icd9_labels <- read_excel(LABELS_XLSX, sheet = "CCS ICD-9-CM-3Level") %>%
  transmute(ICD_9_CM = as.character(ICD_9_CM), icd9_label = tolower(ICD_9_CM_LABEL),
            CCS_ID = as.character(CCS_ID))
icd10_labels <- read_excel(LABELS_XLSX, sheet = "ICD-10-CA-3Level") %>%
  transmute(target = as.character(ICD_10_CA), target_label = tolower(ICD_10_CA_LABEL))
icda8_labels <- read_excel(LABELS_XLSX, sheet = "ICDA-8-3Level") %>%
  transmute(target = as.character(ICDA_8), target_label = tolower(ICDA_8_LABEL))

TRACKS <- list(
  `10_9` = list(tcn = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                find_fn = find_icd10ca_chapter, align = chapter_alignment_10,
                manual_target_col = "ICD-10-CA", manual_sheet = "Validation_ICD9_ICD10",
                excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                target_labels = icd10_labels,
                sim = list(
                  sapbert      = file.path(SAPBERT_BASE, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                  mpnet        = file.path(GEN_BASE, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"),
                  clinicalbert = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"))),
  `8_9`  = list(tcn = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                find_fn = find_icda8_chapter, align = chapter_alignment_8,
                manual_target_col = "ICDA-8", manual_sheet = "Validaion_ICD9_ICD8",
                excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                target_labels = icda8_labels,
                sim = list(
                  sapbert      = file.path(SAPBERT_BASE, "cosine_similarity_matrices_8_9_SapBERT.xlsx"),
                  mpnet        = file.path(GEN_BASE, "cosine_similarity_matrices_8_9_mpnet_base.xlsx"),
                  clinicalbert = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")))
)

long_similarity <- function(path) {
  sheets <- load_similarity_sheets(path)
  purrr::map_dfr(sheets, function(df) {
    id_col <- names(df)[1]
    tgt <- as.character(df[[id_col]])
    purrr::map_dfr(setdiff(names(df), id_col), function(col) {
      tibble(ICD_9_CM = as.character(col), target = tgt, sim = df[[col]])
    })
  }) %>% filter(!is.na(sim))
}

# Non-neural label features. Embeddings smooth over exact word matches, so
# two labels sharing a rare word ("mesothelioma") look no better than ones
# sharing "other". These catch that.
tokenize <- function(x) strsplit(gsub("[^a-z0-9 ]", " ", x), "\\s+")

lexical_features <- function(a, b, idf) {
  ta <- tokenize(a); tb <- tokenize(b)
  n <- length(a)
  jac <- numeric(n); ovl <- numeric(n); idf_ovl <- numeric(n); first_tok <- numeric(n)
  for (i in seq_len(n)) {
    x <- setdiff(ta[[i]], ""); y <- setdiff(tb[[i]], "")
    if (!length(x) || !length(y)) { jac[i] <- 0; ovl[i] <- 0; idf_ovl[i] <- 0; first_tok[i] <- 0; next }
    inter <- intersect(x, y)
    jac[i] <- length(inter) / length(union(x, y))
    ovl[i] <- length(inter) / min(length(x), length(y))
    # rare words count more than "other"/"disease"
    idf_ovl[i] <- if (length(inter)) sum(idf[inter], na.rm = TRUE) / sum(idf[union(x, y)], na.rm = TRUE) else 0
    first_tok[i] <- as.numeric(x[1] == y[1])
  }
  # Edit distance, pairwise on distinct label pairs only.
  # adist(a, b) builds the full a x b matrix and segfaults at this size.
  key <- paste(a, b, sep = "\r")
  uk <- !duplicated(key)
  ua <- a[uk]; ub <- b[uk]
  ued <- mapply(function(x, y) as.numeric(adist(x, y)), ua, ub, USE.NAMES = FALSE)
  ed <- ued[match(key, key[uk])]
  list(lex_jaccard = jac, lex_overlap = ovl, lex_idf_overlap = idf_ovl,
       lex_first_token = first_tok,
       lex_editdist = ed / pmax(nchar(a), nchar(b), 1))
}

for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  cat(sprintf("\n######## building features for track %s ########\n", tr))

  # --- similarity, all models -----------------------------------------
  sims <- list()
  for (mdl in names(tk$sim)) {
    cat(sprintf("  loading similarity: %s\n", mdl))
    s <- long_similarity(tk$sim[[mdl]])
    # forward: rank each target within its ICD-9 code
    s <- s %>% group_by(ICD_9_CM) %>%
      mutate(sim_rank = rank(-sim, ties.method = "first"),
             sim_rel  = sim / max(sim, na.rm = TRUE),
             sim_z    = (sim - mean(sim, na.rm = TRUE)) / (sd(sim, na.rm = TRUE) + 1e-9)) %>%
      ungroup()

    # reverse: rank each ICD-9 code within its target. Catches generic
    # "other/unspecified" targets that look good to everyone.
    s <- s %>% group_by(target) %>%
      mutate(sim_rank_rev = rank(-sim, ties.method = "first"),
             sim_rel_rev  = sim / max(sim, na.rm = TRUE)) %>%
      ungroup()

    names(s)[names(s) == "sim"]          <- paste0("sim_", mdl)
    names(s)[names(s) == "sim_rank"]     <- paste0("simrank_", mdl)
    names(s)[names(s) == "sim_rel"]      <- paste0("simrel_", mdl)
    names(s)[names(s) == "sim_z"]        <- paste0("simz_", mdl)
    names(s)[names(s) == "sim_rank_rev"] <- paste0("simrankrev_", mdl)
    names(s)[names(s) == "sim_rel_rev"]  <- paste0("simrelrev_", mdl)
    sims[[mdl]] <- s
  }

  # --- co-occurrence ---------------------------------------------------
  cooc_raw <- load_cooccurrence_df(file.path(ORIG_BASE, tk$cooc_file)) %>% as.data.frame()
  names(cooc_raw)[names(cooc_raw) == tk$icd9_col]   <- "ICD_9_CM"
  names(cooc_raw)[names(cooc_raw) == tk$target_col] <- "target"
  cooc <- cooc_raw %>%
    transmute(ICD_9_CM = as.character(ICD_9_CM), target = as.character(target),
              cooc_freq = Co_Occurrence_Frequency) %>%
    group_by(ICD_9_CM) %>%
    mutate(cooc_rank  = rank(-cooc_freq, ties.method = "first"),
           cooc_share = cooc_freq / sum(cooc_freq, na.rm = TRUE),
           cooc_rel   = cooc_freq / max(cooc_freq, na.rm = TRUE)) %>%
    ungroup()

  # --- maximal candidate pool -----------------------------------------
  pool <- bind_rows(
    lapply(names(sims), function(m)
      sims[[m]] %>% filter(.data[[paste0("simrank_", m)]] <= MAX_SIM_K) %>% select(ICD_9_CM, target)),
    cooc %>% filter(cooc_rank <= MAX_COOC_N) %>% select(ICD_9_CM, target)
  ) %>% distinct()
  cat(sprintf("  maximal pool: %d pairs over %d ICD-9 codes (%.1f per code)\n",
              nrow(pool), n_distinct(pool$ICD_9_CM), nrow(pool)/n_distinct(pool$ICD_9_CM)))

  feat <- pool
  for (m in names(sims)) feat <- feat %>% left_join(sims[[m]], by = c("ICD_9_CM", "target"))
  feat <- feat %>% left_join(cooc, by = c("ICD_9_CM", "target"))

  # a candidate present in the co-occurrence table but absent from the
  # similarity matrix (or vice versa) is informative, so record presence
  # explicitly and give the missing numeric a neutral floor rather than NA
  feat <- feat %>%
    mutate(has_cooc = as.integer(!is.na(cooc_freq)),
           cooc_freq  = ifelse(is.na(cooc_freq), 0, cooc_freq),
           cooc_share = ifelse(is.na(cooc_share), 0, cooc_share),
           cooc_rel   = ifelse(is.na(cooc_rel), 0, cooc_rel),
           cooc_rank  = ifelse(is.na(cooc_rank), 999, cooc_rank))
  for (m in names(sims)) {
    for (p in c("simrank_", "simrankrev_")) {
      col <- paste0(p, m); feat[[col]] <- ifelse(is.na(feat[[col]]), 999, feat[[col]])
    }
    for (p in c("sim_", "simrel_", "simz_", "simrelrev_")) {
      col <- paste0(p, m); feat[[col]] <- ifelse(is.na(feat[[col]]), 0, feat[[col]])
    }
  }

  # --- ensemble features ("mix the models") ----------------------------
  # Rank fusion rather than raw-score averaging: cosine scales differ between
  # models, ranks do not, so reciprocal rank fusion combines them without
  # needing calibration.
  rank_cols    <- paste0("simrank_",    names(sims))
  rel_cols     <- paste0("simrel_",     names(sims))
  rankrev_cols <- paste0("simrankrev_", names(sims))
  relrev_cols  <- paste0("simrelrev_",  names(sims))
  feat <- feat %>%
    mutate(ens_rrf       = rowSums(1 / (60 + as.matrix(across(all_of(rank_cols))))),
           ens_mean_rel  = rowMeans(across(all_of(rel_cols))),
           ens_max_rel   = do.call(pmax, c(across(all_of(rel_cols)), na.rm = TRUE)),
           ens_min_rel   = do.call(pmin, c(across(all_of(rel_cols)), na.rm = TRUE)),
           ens_best_rank = do.call(pmin, c(across(all_of(rank_cols)), na.rm = TRUE)),
           # models disagreeing is itself informative
           ens_rel_sd    = apply(across(all_of(rel_cols)), 1, sd),

           # reverse direction
           ens_mean_relrev  = rowMeans(across(all_of(relrev_cols))),
           ens_best_rankrev = do.call(pmin, c(across(all_of(rankrev_cols)), na.rm = TRUE)),
           ens_rrf_rev      = rowSums(1 / (60 + as.matrix(across(all_of(rankrev_cols))))),

           # mutual agreement. geometric mean, so one strong direction can't
           # mask a weak one the way a sum would
           ens_mutual_rel  = sqrt(ens_mean_rel * ens_mean_relrev),
           ens_mutual_rank = ens_best_rank + ens_best_rankrev,
           is_mutual_top1  = as.integer(ens_best_rank == 1 & ens_best_rankrev == 1),
           # positive when the target wants this code more than vice versa
           ens_direction_gap = ens_mean_rel - ens_mean_relrev)

  # --- chapter features ------------------------------------------------
  ch9  <- find_icd9cm_chapter(feat$ICD_9_CM)
  cht  <- tk$find_fn(feat$target)
  feat <- feat %>%
    mutate(chapter_icd9 = ch9, chapter_target = cht,
           chapter_distance = compute_chapter_distance(ch9, cht, tk$align),
           chapter_ok = as.integer(!is.na(chapter_distance) & chapter_distance < 1),
           chapter_pair = paste(ch9, cht, sep = "_"))

  # --- lexical features ------------------------------------------------
  cat("  computing lexical features...\n")
  feat <- feat %>%
    left_join(icd9_labels, by = "ICD_9_CM") %>%
    left_join(tk$target_labels, by = "target") %>%
    mutate(icd9_label = ifelse(is.na(icd9_label), "", icd9_label),
           target_label = ifelse(is.na(target_label), "", target_label))
  all_tokens <- c(tokenize(unique(feat$icd9_label)), tokenize(unique(feat$target_label)))
  df_counts <- table(unlist(lapply(all_tokens, unique)))
  idf <- log(length(all_tokens) / (1 + as.numeric(df_counts)))
  names(idf) <- names(df_counts)
  lf <- lexical_features(feat$icd9_label, feat$target_label, idf)
  for (nm in names(lf)) feat[[nm]] <- lf[[nm]]

  # --- per-code context features ---------------------------------------
  feat <- feat %>% group_by(ICD_9_CM) %>%
    mutate(n_candidates = n(),
           ens_rrf_rank = rank(-ens_rrf, ties.method = "first"),
           ens_rrf_rel  = ens_rrf / max(ens_rrf, na.rm = TRUE),
           lex_idf_rank = rank(-lex_idf_overlap, ties.method = "first"),
           # gap to the runner-up: clear winner vs crowded tie
           ens_margin   = ens_mean_rel - max(ens_mean_rel[ens_mean_rel < max(ens_mean_rel)], -1)) %>%
    ungroup()

  # how contested each target is within the pool
  feat <- feat %>% group_by(target) %>%
    mutate(target_n_suitors    = n(),
           target_rank_here    = rank(-ens_mean_rel, ties.method = "first"),
           target_is_best_here = as.integer(target_rank_here == 1)) %>%
    ungroup()

  # --- label (target variable only; never used as a feature) ------------
  manual <- read_excel(VAL_XLSX, sheet = tk$manual_sheet)
  excl   <- read_excel(VAL_XLSX, sheet = tk$excl_sheet)
  excluded <- as.character(excl$`ICD-9-CM`)
  tp_df <- manual %>%
    transmute(ICD_9_CM = as.character(`ICD-9-CM`), target = as.character(.data[[tk$manual_target_col]])) %>%
    filter(!(ICD_9_CM %in% excluded)) %>% distinct() %>% mutate(y = 1L)

  # Keep every code, including ones with no known answer. has_truth marks
  # which codes can be scored; codes without it are what the finished tool
  # actually has to map. Codes on the exclusion list are dropped.
  eval_codes <- unique(tp_df$ICD_9_CM)
  feat <- feat %>%
    filter(!(ICD_9_CM %in% excluded)) %>%
    left_join(tp_df, by = c("ICD_9_CM", "target")) %>%
    mutate(y = ifelse(is.na(y), 0L, 1L),
           has_truth = as.integer(ICD_9_CM %in% eval_codes))

  n_unl <- n_distinct(feat$ICD_9_CM[feat$has_truth == 0])
  cat(sprintf("  codes with truth: %d | without truth (to be predicted): %d\n",
              length(eval_codes), n_unl))

  hit <- feat %>% filter(y == 1) %>% nrow()
  cat(sprintf("  feature table: %d rows, %d positives\n", nrow(feat), hit))
  cat(sprintf("  recall ceiling of maximal pool: %.4f (%d / %d true pairs)\n",
              hit / nrow(tp_df), hit, nrow(tp_df)))
  cat(sprintf("  candidates per code: %.1f\n", nrow(feat) / n_distinct(feat$ICD_9_CM)))

  attr(feat, "n_true_pairs") <- nrow(tp_df)
  saveRDS(feat, file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr)))
  cat(sprintf("  wrote %s\n", file.path(OUT_DIR, sprintf("rerank_features_%s.rds", tr))))
}

cat("\nDone.\n")
