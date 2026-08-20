# stitches the per track cv outputs back together after running both in parallel
#
#   Rscript 12_cv_rerank.R 10_9 &   Rscript 12_cv_rerank.R 8_9 &
#   Rscript 12b_merge_cv_results.R
#
# running 12_ with no arguments does both tracks and writes these files directly,
# in which case this isnt needed

suppressMessages(library(dplyr))
OUT_DIR <- "../results"

parts <- list.files(OUT_DIR, pattern = "^cv_rerank_part_.*\\.rds$", full.names = TRUE)
if (!length(parts)) stop("no cv_rerank_part_*.rds found in ", OUT_DIR)
cat(sprintf("Merging %d track part(s): %s\n", length(parts),
            paste(basename(parts), collapse = ", ")))

p <- lapply(parts, readRDS)
results <- bind_rows(lapply(p, `[[`, "results"))
folds   <- bind_rows(lapply(p, `[[`, "folds"))
imp     <- bind_rows(lapply(p, `[[`, "imp"))
preds   <- bind_rows(lapply(p, `[[`, "preds"))

write.csv(results, file.path(OUT_DIR, "cv_rerank_results.csv"), row.names = FALSE)
write.csv(folds,   file.path(OUT_DIR, "cv_rerank_folds.csv"),   row.names = FALSE)
saveRDS(preds, file.path(OUT_DIR, "cv_rerank_predictions.rds"))

importance <- imp %>% group_by(track, Feature) %>%
  summarise(gain = mean(Gain), .groups = "drop") %>% arrange(track, desc(gain))
write.csv(importance, file.path(OUT_DIR, "cv_rerank_importance.csv"), row.names = FALSE)

cat(sprintf("Merged: %d result rows, %d fold rows, %d prediction rows, tracks: %s\n",
            nrow(results), nrow(folds), nrow(preds),
            paste(unique(results$track), collapse = ", ")))
