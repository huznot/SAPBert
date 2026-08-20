# Breaks performance down by how many correct targets a code has, and checks
# how much of the ICDA-8 task is identity mapping.

suppressMessages({library(readxl); library(dplyr)})
V <- "data/original/ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"
m10 <- read_excel(V, sheet = "Validation_ICD9_ICD10")
m8  <- read_excel(V, sheet = "Validaion_ICD9_ICD8")

a <- m10 %>% transmute(icd9 = as.character(`ICD-9-CM`), tgt = as.character(`ICD-10-CA`))
b <- m8  %>% transmute(icd9 = as.character(`ICD-9-CM`), tgt = as.character(`ICDA-8`))

cat("=== how much is the target just the same number? ===\n")
cat(sprintf("  ICDA-8 : %.1f%% of pairs have target == source code\n",
            100*mean(b$icd9 == b$tgt)))
cat(sprintf("  ICDA-8 : %.1f%% of codes have at least one identical-code target\n",
            100*mean(b %>% group_by(icd9) %>% summarise(x=any(icd9==tgt)) %>% pull(x))))
cat(sprintf("  ICD-10 : %.1f%% identical (impossible, alphanumeric)\n",
            100*mean(a$icd9 == a$tgt)))

cat("\n=== target space size ===\n")
cat(sprintf("  ICD-10-CA distinct targets used: %d\n", n_distinct(a$tgt)))
cat(sprintf("  ICDA-8   distinct targets used: %d\n", n_distinct(b$tgt)))

cat("\n=== targets per source code ===\n")
for (nm in c("ICD-10-CA", "ICDA-8")) {
  d <- if (nm == "ICD-10-CA") a else b
  x <- d %>% count(icd9)
  cat(sprintf("  %-10s mean %.2f, median %g, max %d, %% with >1: %.1f%%\n",
              nm, mean(x$n), median(x$n), max(x$n), 100*mean(x$n > 1)))
}

cat("\n=== how complete are the auto-accepted mappings per code? ===\n")
p <- readRDS("results/cv_rerank_predictions.rds")
ops <- read.csv("results/precision_coverage_operating_points.csv")
for (tr in c("10_9", "8_9")) {
  tau <- (ops %>% filter(track == tr, precision_target == 0.95))$tau
  d <- p %>% filter(track == tr)
  per <- d %>% group_by(ICD_9_CM) %>%
    summarise(n_true = sum(y), n_hit = sum(y == 1 & .p_score >= tau), .groups = "drop") %>%
    filter(n_true > 0)
  cat(sprintf("\n  %s (tau %.2f)\n", tr, tau))
  cat(sprintf("    codes where ALL correct targets auto-accepted : %.1f%%\n",
              100*mean(per$n_hit == per$n_true)))
  cat(sprintf("    codes where SOME but not all                  : %.1f%%\n",
              100*mean(per$n_hit > 0 & per$n_hit < per$n_true)))
  cat(sprintf("    codes where NONE                              : %.1f%%\n",
              100*mean(per$n_hit == 0)))
  s <- per %>% filter(n_true == 1)
  cat(sprintf("    single-target codes fully solved              : %.1f%% (n=%d)\n",
              100*mean(s$n_hit == 1), nrow(s)))
  mt <- per %>% filter(n_true > 1)
  if (nrow(mt)) cat(sprintf("    multi-target codes fully solved              : %.1f%% (n=%d)\n",
              100*mean(mt$n_hit == mt$n_true), nrow(mt)))
}
