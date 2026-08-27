# combines the per condition grids from 07_ into summary tables
#   full_grid_all.csv          every point
#   full_grid_best.csv         best per condition x track
#   full_grid_stripping.csv    base vs stripped
#   full_grid_family.csv       best per model family
#   full_grid_sensitivity.csv  how much the ranking depends on the operating point
#
# the sensitivity table is the point, comparing models at one arbitrary operating
# point cant tell "A is better" from "A wins at this one spot"
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

GRID_DIR <- "results/full_grid"
OUT_DIR  <- "results"

files <- list.files(GRID_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no per-condition grid CSVs in ", GRID_DIR,
                             " -- run 07_full_grid_comparison.R first")
cat(sprintf("Assembling %d condition file(s): %s\n", length(files),
            paste(basename(files), collapse = ", ")))

all_grid <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))

# every condition must have run the same grid or best-f1 is comparing different
# amounts of search
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

# paired within family and track, regenerated arms only, so base and stripped
# differ only by the stripping
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

# win rate at matched points. better read than best-f1 for a small effect since
# one arm can win just by finding a luckier corner of the grid
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

# rank by best over grid f1 vs f1 at the single point 06_ used. if they disagree
# the single point comparison was measuring the point
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

# figure 1 in the report. built here rather than in 03_ so it comes from the
# same table as table 1 and covers all three models
suppressMessages(library(ggplot2))
KEEP <- c("ClinicalBERT-original" = "ClinicalBERT",
          "SapBERT-base"          = "SapBERT",
          "mpnet-base"            = "all-mpnet-base-v2")
fig <- best %>% filter(model %in% names(KEEP)) %>%
  mutate(Model = factor(unname(KEEP[model]), levels = unname(KEEP)),
         track_label = factor(ifelse(track == "10_9",
                                     "ICD-9-CM to ICD-10-CA", "ICD-9-CM to ICDA-8"),
                              levels = c("ICD-9-CM to ICD-10-CA", "ICD-9-CM to ICDA-8"))) %>%
  select(track_label, Model, F1 = f1, Accuracy = accuracy) %>%
  tidyr::pivot_longer(c(F1, Accuracy), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("F1", "Accuracy")))
if (nrow(fig)) {
  g <- ggplot(fig, aes(metric, value, fill = Model)) +
    geom_col(position = position_dodge(0.8), width = 0.72) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(0.8), vjust = -0.45, size = 3.1) +
    facet_wrap(~track_label) +
    scale_fill_manual(values = c("ClinicalBERT" = "#4E79A7", "SapBERT" = "#E1873C",
                                 "all-mpnet-base-v2" = "#59A14F")) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Score", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, "plot_f1_accuracy_comparison.png"), g,
         width = 9, height = 4.6, dpi = 150)
  cat("\nwrote results/plot_f1_accuracy_comparison.png\n")
}
