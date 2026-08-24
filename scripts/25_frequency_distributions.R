suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})

ORIG <- "data/original"
OUT  <- "results"

TRACKS <- list(
  `10_9` = list(
    label = "ICD-9-CM to ICD-10-CA",
    sim   = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
    cooc  = file.path(ORIG, "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"),
    icd9  = "ICD_9_CM_Code3", tgt = "ICD_10_CA_Code3"),
  `8_9` = list(
    label = "ICD-9-CM to ICDA-8",
    sim   = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx"),
    cooc  = file.path(ORIG, "Co_occurrence/icd_8_9_co_occurrence_3d.xSlsx"),
    icd9  = "ICD_9_CM_Code", tgt = "ICDA_8_Code")
)

max_sim_per_code <- function(path) {
  bind_rows(lapply(excel_sheets(path), function(s) {
    d <- read_excel(path, sheet = s)
    names(d)[1] <- "target"
    d %>% mutate(target = as.character(target)) %>%
      pivot_longer(-target, names_to = "icd9", values_to = "sim")
  })) %>%
    group_by(icd9) %>%
    slice_max(sim, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(icd9, best_target = target, max_similarity = sim)
}

top_cooc_per_code <- function(path, icd9_col, tgt_col) {
  read_excel(path) %>%
    transmute(icd9 = as.character(.data[[icd9_col]]),
              target = as.character(.data[[tgt_col]]),
              freq = Co_Occurrence_Frequency) %>%
    group_by(icd9) %>%
    slice_max(freq, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(top_target = target, top_frequency = freq)
}

describe <- function(x, name) {
  cat(sprintf("  %-18s n %d | min %.4g | q1 %.4g | median %.4g | mean %.4g | q3 %.4g | max %.4g | sd %.4g\n",
              name, length(x), min(x), quantile(x, .25), median(x),
              mean(x), quantile(x, .75), max(x), sd(x)))
}

sim_all <- list(); cooc_all <- list()

for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  cat(sprintf("\n===== %s =====\n", tk$label))

  ms <- max_sim_per_code(tk$sim) %>% mutate(track = tr, track_label = tk$label)
  cat(sprintf("\nmaximum cosine similarity, %d ICD-9-CM codes\n", nrow(ms)))
  describe(ms$max_similarity, "max similarity")
  # cosine similarity can come back a hair over 1 from floating point, which
  # would fall outside the bins and be dropped. clamp before binning
  n_perfect <- sum(ms$max_similarity >= 1 - 1e-6)
  cat(sprintf("  identical-label matches (similarity ~= 1): %d (%.1f%%)\n",
              n_perfect, 100 * n_perfect / nrow(ms)))
  brk <- seq(0, 1, by = 0.05)
  tb <- table(cut(pmin(ms$max_similarity, 1), breaks = brk, include.lowest = TRUE))
  tb <- tb[tb > 0]
  for (i in seq_along(tb))
    cat(sprintf("    %-14s %4d  %5.1f%%  %s\n", names(tb)[i], tb[i],
                100 * tb[i] / nrow(ms), strrep("#", round(50 * tb[i] / max(tb)))))

  tc <- top_cooc_per_code(tk$cooc, tk$icd9, tk$tgt) %>%
    mutate(track = tr, track_label = tk$label)
  cat(sprintf("\nfrequency of the single most co-occurring target, %d ICD-9-CM codes\n", nrow(tc)))
  describe(tc$top_frequency, "top frequency")
  cat(sprintf("  codes with no co-occurrence data: %d\n", nrow(ms) - nrow(tc)))
  qb <- c(0, 10, 25, 50, 100, 250, 500, 1000, 5000, 10000, Inf)
  tb2 <- table(cut(tc$top_frequency, breaks = qb, include.lowest = TRUE))
  tb2 <- tb2[tb2 > 0]
  for (i in seq_along(tb2))
    cat(sprintf("    %-16s %4d  %5.1f%%  %s\n", names(tb2)[i], tb2[i],
                100 * tb2[i] / nrow(tc), strrep("#", round(50 * tb2[i] / max(tb2)))))

  sim_all[[tr]] <- ms; cooc_all[[tr]] <- tc
}

sim_all  <- bind_rows(sim_all)
cooc_all <- bind_rows(cooc_all)
write.csv(sim_all,  file.path(OUT, "freq_dist_max_similarity.csv"), row.names = FALSE)
write.csv(cooc_all, file.path(OUT, "freq_dist_top_cooccurrence.csv"), row.names = FALSE)

sim_all <- sim_all %>% mutate(max_similarity_plot = pmin(max_similarity, 1))

g1 <- ggplot(sim_all, aes(max_similarity_plot)) +
  geom_histogram(binwidth = 0.02, fill = "#3182bd", colour = "white") +
  facet_wrap(~track_label, ncol = 1, scales = "free_y") +
  labs(x = "maximum cosine similarity for an ICD-9-CM code", y = "number of codes",
       title = "Distribution of maximum similarity score (ClinicalBERT)") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT, "plot_freq_dist_max_similarity.png"), g1, width = 7, height = 6, dpi = 150)

g2 <- ggplot(cooc_all, aes(top_frequency)) +
  geom_histogram(bins = 40, fill = "#31a354", colour = "white") +
  scale_x_log10() +
  facet_wrap(~track_label, ncol = 1, scales = "free_y") +
  labs(x = "co-occurrence count of the most frequent target (log scale)",
       y = "number of codes",
       title = "Distribution of the most frequent co-occurring code") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT, "plot_freq_dist_top_cooccurrence.png"), g2, width = 7, height = 6, dpi = 150)

cat("\nwrote results/freq_dist_max_similarity.csv\n")
cat("wrote results/freq_dist_top_cooccurrence.csv\n")
cat("wrote results/plot_freq_dist_max_similarity.png\n")
cat("wrote results/plot_freq_dist_top_cooccurrence.png\n")
