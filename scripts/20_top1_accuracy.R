# top-1 accuracy, for context against the published gpt-4 results
#
# that study reports whether the single predicted code matched the reference.
# this system predicts a set per code and is scored with precision/recall/f1 so
# the numbers arent comparable as they stand. this computes the closest thing:
# for each icd-9 code, is the top scoring target one of its correct ones
#
# read the caveats in the output before putting this next to their 84.7%,
# direction and reference standard and code subset all differ
#
# uses the out of fold predictions from 12_ so it is held out

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

OUT_DIR <- "results"
preds <- readRDS(file.path(OUT_DIR, "cv_rerank_predictions.rds"))

rows <- list()
for (tr in unique(preds$track)) {
  d <- preds %>% filter(track == tr)

  top1 <- d %>% group_by(ICD_9_CM) %>%
    slice_max(.p_score, n = 1, with_ties = FALSE) %>% ungroup()
  acc1 <- mean(top1$y == 1)

  # how often a correct target appears in the top k
  topk <- sapply(c(1, 3, 5), function(k) {
    d %>% group_by(ICD_9_CM) %>%
      slice_max(.p_score, n = k, with_ties = FALSE) %>%
      summarise(hit = as.integer(any(y == 1)), .groups = "drop") %>%
      summarise(mean(hit)) %>% pull()
  })

  n_codes <- n_distinct(d$ICD_9_CM)
  n_tgt <- d %>% filter(y == 1) %>% count(ICD_9_CM) %>% summarise(mean(n)) %>% pull()

  cat(sprintf("\n=== %s ===\n", tr))
  cat(sprintf("  codes evaluated              : %d\n", n_codes))
  cat(sprintf("  correct targets per code     : %.2f on average\n", n_tgt))
  cat(sprintf("  top-1 accuracy               : %.1f%%\n", 100*acc1))
  cat(sprintf("  correct target within top 3  : %.1f%%\n", 100*topk[2]))
  cat(sprintf("  correct target within top 5  : %.1f%%\n", 100*topk[3]))

  rows[[length(rows)+1]] <- tibble(track = tr, n_codes = n_codes,
    mean_targets_per_code = round(n_tgt, 2),
    top1_accuracy = round(acc1, 4),
    top3_accuracy = round(topk[2], 4),
    top5_accuracy = round(topk[3], 4))
}

res <- bind_rows(rows)
write.csv(res, file.path(OUT_DIR, "top1_accuracy.csv"), row.names = FALSE)

cat("\n=== why this is not a head-to-head with the GPT-4 study ===\n")
cat("  direction : that study translates ICD-10-CA -> ICD-9-CM, this goes ICD-9-CM -> ICD-10-CA\n")
cat("  reference : the GPT-4 study scores against the CIHI crosswalk, this against the local manual one\n")
cat("  codes     : it uses 1,272 chronic-disease (Elixhauser) codes, this uses 345/302\n")
cat("  structure : the CIHI crosswalk is one-to-one; here a code averages >1 correct target,\n")
cat("              so top-1 is an easier target for them and a harder one here\n")
cat("  contamination : public crosswalks may sit in GPT-4 training data, which they flag.\n")
cat("              Co-occurrence counts from local records cannot leak that way.\n")
cat("\nWritten to results/top1_accuracy.csv\n")
