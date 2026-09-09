source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

ORIG <- "data/original"; GEN <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
COOC <- file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx")

man  <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
excl <- read_excel(VAL, sheet = "Validation_ICD9_ICD10_Excld")
ccs  <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)
cooc <- load_cooccurrence_df(COOC)
res  <- read.csv(out_path("absolute_threshold_grid.csv"), stringsAsFactors = FALSE)

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
            mpnet        = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx")

out <- list()

## 1.
cat("\n\n=== 1. co-occurrence depth, at each model's best threshold ===\n")
depth <- do.call(rbind, lapply(names(MODELS), function(m) {
  d <- res[res$model == m & res$mode == "absolute", ]
  b <- d[which.max(d$f1), ]
  s <- d[d$threshold == b$threshold & d$flag == 4, ]
  s <- s[order(s$top_n), ]
  data.frame(Model = m, `Top N` = s$top_n,
             Precision = sprintf("%.3f", s$precision),
             Recall    = sprintf("%.3f", s$recall),
             F1        = sprintf("%.3f", s$f1),
             check.names = FALSE)
}))
print(depth, row.names = FALSE)
out$cooccurrence_depth <- depth

## 2. 
cat("\n\n=== 2. where the correct answers come from ===\n")
src <- list()
for (m in names(MODELS)) {
  d <- res[res$model == m & res$mode == "absolute", ]
  b <- d[which.max(d$f1), ]
  long <- purrr::map_dfr(load_similarity_sheets(file.path(GEN, MODELS[[m]])), function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble(ICD_9_CM = as.character(c1), ICD_10_CA = as.character(df[[id]]), Similarity = df[[c1]]))
  })
  sim <- long %>% filter(!is.na(Similarity), Similarity >= b$threshold)
  co  <- get_cooccurrence_codes_from_df(cooc, b$top_n, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
  mg  <- merge_and_flag(sim, co, "ICD_10_CA", find_icd10ca_chapter, chapter_alignment_10)
  au  <- select_rows_by_flags(mg, b$flag)
  au$branch <- ifelse(au$highest_similarity_flag == 1 & au$highest_cooccurrence_flag == 1, "both branches",
               ifelse(au$highest_similarity_flag == 1, "cosine only",
               ifelse(au$highest_cooccurrence_flag == 1, "co-occurrence only", "in both lists")))
  key <- paste(as.character(man$`ICD-9-CM`), as.character(man$`ICD-10-CA`))
  au$correct <- paste(as.character(au$ICD_9_CM), as.character(au$ICD_10_CA)) %in% key
  t <- table(au$branch, ifelse(au$correct, "correct", "wrong"))
  src[[m]] <- data.frame(Model = m, Source = rownames(t),
                         Correct = as.integer(t[, "correct"]),
                         Wrong   = as.integer(t[, "wrong"]),
                         `Hit rate` = sprintf("%.0f%%", 100 * t[, "correct"] / rowSums(t)),
                         check.names = FALSE)
}
srcd <- do.call(rbind, src)
print(srcd, row.names = FALSE)
out$answer_source <- srcd

## 3. 
cat("\n\n=== 3. correct codes per icd-9 code, the recall ceiling ===\n")
n <- table(as.character(man$`ICD-9-CM`))
tg <- data.frame(`Correct codes for one ICD-9 code` = c("1", "2", "3", "4 or more"),
                 `ICD-9 codes` = c(sum(n == 1), sum(n == 2), sum(n == 3), sum(n >= 4)),
                 Share = sprintf("%.0f%%", 100 * c(mean(n == 1), mean(n == 2), mean(n == 3), mean(n >= 4))),
                 check.names = FALSE)
print(tg, row.names = FALSE)
cat(sprintf("\n%d pairs over %d codes, median %d. emitting one correct code each caps recall at %.3f\n",
            nrow(man), length(n), median(n), length(n) / nrow(man)))
out$targets_per_code <- tg

## 4. 
cat("\n\n=== 4. rank of the correct code ===\n")
pc <- read.csv(out_path("per_code_rank.csv"), stringsAsFactors = FALSE)
band <- function(r) ifelse(is.na(r), "not in column",
                    ifelse(r == 1, "1", ifelse(r <= 5, "2 to 5",
                    ifelse(r <= 20, "6 to 20", "over 20"))))
pc$band <- factor(band(pc$best_rank), levels = c("1", "2 to 5", "6 to 20", "over 20", "not in column"))
t4 <- table(pc$model, pc$band)
rk <- data.frame(Model = rownames(t4), as.data.frame.matrix(t4), check.names = FALSE)
rk$`rank 1 share` <- sprintf("%.0f%%", 100 * t4[, "1"] / rowSums(t4))
print(rk, row.names = FALSE)
out$correct_rank <- rk

## 5.
cat("\n\n=== 5. top scoring code deleted by the chapter filter ===\n")
t5 <- table(pc$model, ifelse(pc$top_survives_chapter_filter, "kept", "deleted"))
cf <- data.frame(Model = rownames(t5), Deleted = t5[, "deleted"], Kept = t5[, "kept"],
                 Share = sprintf("%.0f%%", 100 * t5[, "deleted"] / rowSums(t5)),
                 check.names = FALSE)
print(cf, row.names = FALSE)
out$chapter_filter <- cf

for (nm in names(out)) write.csv(out[[nm]], out_path(paste0("trend_", nm, ".csv")), row.names = FALSE)
cat("\nwrote", length(out), "tables to results/tables/\n")
