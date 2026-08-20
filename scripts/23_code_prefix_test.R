# does putting the code number in the embedding text help or hurt
#
# we were embedding "250 diabetes mellitus" rather than "diabetes mellitus". the
# numbers mean nothing across coding systems, icd-9 250 and icd-10 E11 are the
# same concept with nothing in common lexically, so they may just be noise
#
# compares recall ceiling of the pool at a few top-K values plus top-1 hit rate
# on similarity alone
source("pipeline_lib.R")

ORIG <- "../data/original"
GEN  <- "../data/generated"
SAP  <- "../data/sapbert"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

TRACKS <- list(
  `10_9` = list(sheet = "Validation_ICD9_ICD10", excl = "Validation_ICD9_ICD10_Excld",
                tcol = "ICD-10-CA"),
  `8_9`  = list(sheet = "Validaion_ICD9_ICD8", excl = "Validation_ICD9_ICD8_Excld",
                tcol = "ICDA-8")
)

ARMS <- list(
  `sapbert with code`  = "cosine_similarity_matrices_%s_SapBERT.xlsx",
  `sapbert label only` = "cosine_similarity_matrices_%s_sapbert_base_nocode.xlsx",
  `mpnet with code`    = "cosine_similarity_matrices_%s_mpnet_base.xlsx",
  `mpnet label only`   = "cosine_similarity_matrices_%s_mpnet_base_nocode.xlsx"
)

long_sim <- function(path) {
  sheets <- load_similarity_sheets(path)
  purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]; tgt <- as.character(df[[id]])
    purrr::map_dfr(setdiff(names(df), id), function(col)
      tibble(ICD_9_CM = as.character(col), target = tgt, sim = df[[col]]))
  }) %>% filter(!is.na(sim))
}

rows <- list()
for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  manual <- read_excel(VAL, sheet = tk$sheet)
  excl <- as.character(read_excel(VAL, sheet = tk$excl)$`ICD-9-CM`)
  truth <- manual %>%
    transmute(ICD_9_CM = as.character(`ICD-9-CM`), target = as.character(.data[[tk$tcol]])) %>%
    filter(!(ICD_9_CM %in% excl)) %>% distinct()
  n_true <- nrow(truth)

  cat(sprintf("\n######## %s (%d true pairs) ########\n", tr, n_true))
  for (arm in names(ARMS)) {
    p <- file.path(if (grepl("^sapbert with", arm)) SAP else GEN, sprintf(ARMS[[arm]], tr))
    if (!file.exists(p)) { cat(sprintf("  %-20s missing\n", arm)); next }
    s <- long_sim(p) %>% group_by(ICD_9_CM) %>%
      mutate(rk = rank(-sim, ties.method = "first")) %>% ungroup()

    hit1 <- s %>% filter(rk == 1) %>% semi_join(truth, by = c("ICD_9_CM","target")) %>%
      nrow() / n_distinct(truth$ICD_9_CM)
    ceilings <- sapply(c(10, 25, 50), function(k)
      nrow(semi_join(truth, s %>% filter(rk <= k), by = c("ICD_9_CM","target"))) / n_true)

    cat(sprintf("  %-20s top-1 hit %.3f | pool recall K10 %.3f K25 %.3f K50 %.3f\n",
                arm, hit1, ceilings[1], ceilings[2], ceilings[3]))
    rows[[length(rows)+1]] <- tibble(track = tr, arm = arm, top1_hit = round(hit1, 4),
      recall_k10 = round(ceilings[1], 4), recall_k25 = round(ceilings[2], 4),
      recall_k50 = round(ceilings[3], 4))
  }
}

res <- bind_rows(rows)
write.csv(res, "../results/code_prefix_test.csv", row.names = FALSE)

cat("\n===== label only minus with code =====\n")
for (m in c("sapbert", "mpnet")) {
  for (tr in unique(res$track)) {
    a <- res %>% filter(track == tr, arm == paste(m, "with code"))
    b <- res %>% filter(track == tr, arm == paste(m, "label only"))
    if (!nrow(a) || !nrow(b)) next
    cat(sprintf("  %-8s %-5s top-1 %+.4f | K10 %+.4f | K50 %+.4f\n", m, tr,
                b$top1_hit - a$top1_hit, b$recall_k10 - a$recall_k10,
                b$recall_k50 - a$recall_k50))
  }
}
