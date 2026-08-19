# Assembles the per-condition full-grid CSVs written by
# scripts/07_full_grid_comparison.R into the summary tables the report reads.
#
# Outputs (all in results/):
#   full_grid_all.csv          every (condition, track, threshold, top_n, flag) point
#   full_grid_best.csv         best point per condition x track, by F1
#   full_grid_stripping.csv    base vs stripped, paired within family x track (task 1)
#   full_grid_family.csv       best per model family x track (task 2)
#   full_grid_sensitivity.csv  how much the ranking depends on the operating point
#
# The sensitivity table is the reason this exists at all: the earlier
# single-operating-point comparison could not distinguish "model A is better"
# from "model A happens to be better at the one point we scored", so alongside
# each condition best F1 it records that condition F1 spread across the grid,
# and how it ranks when every condition is scored at one shared point instead.

source("pipeline_lib.R")

GRID_DIR <- "../results/full_grid"
OUT_DIR  <- "../results"

files <- list.files(GRID_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no per-condition grid CSVs in ", GRID_DIR,
                             " -- run 07_full_grid_comparison.R first")
cat(sprintf("Assembling %d condition file(s): %s\n", length(files),
            paste(basename(files), collapse = ", ")))

all_grid <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))

# Guard against assembling a partial set: every condition must have run the
# identical grid, otherwise "best F1 per condition" is comparing different
# amounts of search and the comparison is not fair.
pts <- all_grid %>% count(model, track, name = "n_points")
if (length(unique(pts$n_points)) != 1) {
  cat("\nWARNING: conditions did not all run the same number of grid points --\n")
  cat("best-F1 comparisons below are NOT apples-to-apples.\n")
  print(as.data.frame(pts))
}
cat(sprintf("%d rows, %d conditions, %d points per condition x track\n",
            nrow(all_grid), length(unique(all_grid$model)), max(pts$n_points)))

write.csv(all_grid, file.path(OUT_DIR, "full_grid_all.csv"), row.names = FALSE)

# --- best point per condition x track ---------------------------------
best <- all_grid %>%
  group_by(track, model) %>%
  slice_max(order_by = f1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(track, desc(f1))
write.csv(best, file.path(OUT_DIR, "full_grid_best.csv"), row.names = FALSE)

# --- task 1: does filler-word stripping help? -------------------------
# Paired within family and track, using only the regenerated arms, where base
# and stripped share a generation method and so differ ONLY by the stripping.
paired <- all_grid %>%
  filter(model %in% c("ClinicalBERT-base", "ClinicalBERT-stripped",
                      "SapBERT-base-regen", "SapBERT-stripped",
                      "mpnet-base", "mpnet-stripped"))

stripping <- paired %>%
  group_by(track, family, stripping) %>%
  summarise(best_f1 = max(f1), best_accuracy = accuracy[which.max(f1)],
            mean_f1 = round(mean(f1), 4), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = stripping, values_from = c(best_f1, best_accuracy, mean_f1)) %>%
  mutate(delta_best_f1 = round(best_f1_stripped - best_f1_base, 4),
         delta_mean_f1 = round(mean_f1_stripped - mean_f1_base, 4)) %>%
  arrange(track, family)
write.csv(stripping, file.path(OUT_DIR, "full_grid_stripping.csv"), row.names = FALSE)

# Point-by-point win rate: whether stripping helps at a matched
# (threshold, top_n, flag) point, not just at each arm own best point. This
# is the more honest read on a change whose effect is small -- one arm can
# win on best-F1 purely by finding a slightly luckier corner of the grid.
paired_points <- paired %>%
  select(track, family, stripping, similarity_threshold, top_n, flag_combination, f1) %>%
  tidyr::pivot_wider(names_from = stripping, values_from = f1) %>%
  filter(!is.na(base), !is.na(stripped)) %>%
  group_by(track, family) %>%
  summarise(n_points = n(),
            stripped_better = sum(stripped > base),
            tied            = sum(stripped == base),
            base_better     = sum(stripped < base),
            mean_delta_f1   = round(mean(stripped - base), 4),
            .groups = "drop")
write.csv(paired_points, file.path(OUT_DIR, "full_grid_stripping_paired.csv"), row.names = FALSE)

# --- task 2: general-purpose (mpnet) vs domain-specific ---------------
family_best <- all_grid %>%
  group_by(track, family) %>%
  slice_max(order_by = f1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(track, family, model, similarity_threshold, top_n, flag_combination,
         precision, recall, f1, accuracy) %>%
  arrange(track, desc(f1))
write.csv(family_best, file.path(OUT_DIR, "full_grid_family.csv"), row.names = FALSE)

# --- how operating-point-dependent is the ranking? --------------------
# For each track, rank conditions by best-over-grid F1, and separately by F1
# at the single shared point 06_extended_comparison.R used. If the two
# rankings disagree, the single-point comparison was measuring the point.
SHARED_POINT <- list(`10_9` = list(thr = 0.995, tn = 30),
                     `8_9`  = list(thr = 0.99,  tn = 5))

shared <- bind_rows(lapply(names(SHARED_POINT), function(tr) {
  sp <- SHARED_POINT[[tr]]
  all_grid %>%
    filter(track == tr, similarity_threshold == sp$thr, top_n == sp$tn) %>%
    group_by(track, model) %>%
    slice_max(order_by = f1, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(track, model, shared_point_f1 = f1,
              shared_point = sprintf("thr=%.3f, top_n=%d", sp$thr, sp$tn))
}))

sensitivity <- best %>%
  select(track, model, family, best_f1 = f1, best_accuracy = accuracy,
         best_threshold = similarity_threshold, best_top_n = top_n,
         best_flag = flag_combination) %>%
  left_join(all_grid %>% group_by(track, model) %>%
              summarise(f1_min = min(f1), f1_max = max(f1),
                        f1_median = round(median(f1), 4), .groups = "drop"),
            by = c("track", "model")) %>%
  left_join(shared, by = c("track", "model")) %>%
  group_by(track) %>%
  mutate(rank_by_best     = rank(-best_f1, ties.method = "min"),
         rank_at_shared   = rank(-shared_point_f1, ties.method = "min"),
         rank_changed     = rank_by_best != rank_at_shared) %>%
  ungroup() %>%
  arrange(track, rank_by_best)
write.csv(sensitivity, file.path(OUT_DIR, "full_grid_sensitivity.csv"), row.names = FALSE)

cat("\n=== Best per condition x track (full grid) ===\n")
print(as.data.frame(best %>% select(track, model, similarity_threshold, top_n,
                                    flag_combination, precision, recall, f1, accuracy)))
cat("\n=== Filler-word stripping, best-point and paired-point (task 1) ===\n")
print(as.data.frame(stripping))
print(as.data.frame(paired_points))
cat("\n=== Best per model family (task 2) ===\n")
print(as.data.frame(family_best))
cat("\n=== Operating-point sensitivity ===\n")
print(as.data.frame(sensitivity %>% select(track, model, best_f1, f1_median, f1_min, f1_max,
                                           shared_point_f1, rank_by_best, rank_at_shared,
                                           rank_changed)))
cat("\nDone.\n")
