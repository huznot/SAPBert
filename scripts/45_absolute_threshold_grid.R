source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

ORIG <- "data/original"
GEN  <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- out_path("absolute_threshold_grid.csv")

MODELS <- list(
  ClinicalBERT = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
  SapBERT      = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
  mpnet        = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx"
)

abs_thresholds <- c(0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99)
rel_thresholds <- c(0.95, 0.99, 0.995, 0.999)
top_ns     <- c(3, 5, 10, 15, 20, 25, 30)
flags      <- 1:4

SCENARIO <- c(
  "1" = "top cosine or top co-occurrence",
  "2" = "top cosine or in both lists",
  "3" = "top co-occurrence or in both lists",
  "4" = "any of the three"
)

sim_long <- function(sheets_list, target_col_name) {
  purrr::map_dfr(sheets_list, function(df) {
    id_col <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id_col), function(col) {
      tibble(ICD_9_CM = as.character(col),
             !!target_col_name := as.character(df[[id_col]]),
             Similarity = df[[col]])
    })
  }) %>% filter(!is.na(Similarity))
}

manual <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
excl   <- read_excel(VAL, sheet = "Validation_ICD9_ICD10_Excld")
ccs    <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
cooc   <- load_cooccurrence_df(file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"))

rows <- list()
for (mdl in names(MODELS)) {
  cat("\n", mdl, "\n", sep = "")
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[mdl]]))
  long <- sim_long(sheets, "ICD_10_CA")
  rng <- range(long$Similarity, na.rm = TRUE)
  cat(sprintf("  score range %.4f to %.4f, median %.4f\n", rng[1], rng[2], median(long$Similarity)))

  colmax <- long %>% group_by(ICD_9_CM) %>% mutate(mx = max(Similarity, na.rm = TRUE)) %>% ungroup()

  grid <- rbind(data.frame(mode = "absolute", thr = abs_thresholds),
                data.frame(mode = "relative", thr = rel_thresholds))

  for (g in seq_len(nrow(grid))) {
    md <- grid$mode[g]; thr <- grid$thr[g]
    sim <- if (md == "absolute") long %>% filter(Similarity >= thr) else
             colmax %>% filter(Similarity >= thr * mx) %>% select(-mx)
    kept_pairs <- nrow(sim)
    kept_cols  <- length(unique(sim$ICD_9_CM))
    cat(sprintf("  %-8s %.3f keeps %d pairs across %d of %d icd-9 codes\n",
                md, thr, kept_pairs, kept_cols, length(unique(long$ICD_9_CM))))
    if (kept_pairs == 0) next
    for (tn in top_ns) {
      co <- get_cooccurrence_codes_from_df(cooc, tn, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
      m  <- merge_and_flag(sim, co, "ICD_10_CA", find_icd10ca_chapter, chapter_alignment_10)
      for (fl in flags) {
        auto <- select_rows_by_flags(m, fl)
        fin  <- validate_mapping(manual, auto, ccs, excl, "ICD-10-CA", "ICD_10_CA")
        mt   <- calculate_performance_metrics(fin)$overall
        rows[[length(rows) + 1]] <- data.frame(
          model = mdl, mode = md, threshold = thr, top_n = tn, flag = fl,
          scenario = unname(SCENARIO[as.character(fl)]),
          precision = mt$overall_precision, recall = mt$overall_recall,
          f1 = mt$overall_f1_score, accuracy = mt$overall_accuracy,
          codes_emitted = nrow(auto),
          tp = sum(fin$`True Positive`, na.rm = TRUE),
          fp = sum(fin$`False Positive`, na.rm = TRUE),
          fn = sum(fin$`False Negative`, na.rm = TRUE),
          sim_pairs_kept = kept_pairs, icd9_with_any_sim = kept_cols,
          stringsAsFactors = FALSE)
      }
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, OUT, row.names = FALSE)
cat("\nwrote", OUT, "with", nrow(res), "rows\n\n")

best <- do.call(rbind, lapply(split(res, list(res$model, res$mode)), function(d)
  if (nrow(d)) d[which.max(d$f1), ] else NULL))
best <- best[order(best$model, best$mode), ]
cat("best f1 per model, absolute cutoff against the relative multiplier, same label only matrices\n")
print(best[, c("model","mode","threshold","top_n","scenario","precision","recall","f1","tp","fp","fn")],
      row.names = FALSE)

cat("\nbest scenario at each absolute threshold, per model\n")
a <- res[res$mode == "absolute", ]
bt <- do.call(rbind, lapply(split(a, list(a$model, a$threshold)), function(d)
  if (nrow(d)) d[which.max(d$f1), ] else NULL))
bt <- bt[order(bt$model, bt$threshold), ]
print(bt[, c("model","threshold","top_n","flag","precision","recall","f1","icd9_with_any_sim")],
      row.names = FALSE)
