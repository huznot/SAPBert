source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

library(dplyr)
library(tidyr)
library(ggplot2)

OUT_DIR <- "results"

best_by_model <- read.csv(file.path(OUT_DIR, "best_by_model.csv"))
ccs_10_9 <- read.csv(file.path(OUT_DIR, "ccs_breakdown_10_9.csv"))
ccs_8_9  <- read.csv(file.path(OUT_DIR, "ccs_breakdown_8_9.csv"))

track_labels <- c("10_9" = "ICD-9-CM -> ICD-10-CA", "8_9" = "ICD-9-CM -> ICDA-8")

comparison_long <- best_by_model %>%
  select(track, model, f1, accuracy) %>%
  pivot_longer(cols = c(f1, accuracy), names_to = "metric", values_to = "value") %>%
  mutate(
    track = track_labels[track],
    metric = recode(metric, f1 = "F1", accuracy = "Accuracy")
  )

p1 <- ggplot(comparison_long, aes(x = metric, y = value, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.5) +
  facet_wrap(~track) +
  scale_fill_manual(values = c("ClinicalBERT" = "#4C72B0", "SapBERT" = "#DD8452")) +
  labs(title = "SapBERT vs ClinicalBERT: best-combination F1 and Accuracy",
       subtitle = "Best parameter combination (similarity threshold x top-N x flag rule) selected by F1 per model/track",
       x = NULL, y = "Score", fill = "Model") +
  ylim(0, max(comparison_long$value) * 1.2) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "plot_f1_accuracy_comparison.png"), p1, width = 9, height = 5.5, dpi = 150)

make_ccs_breakdown_plot <- function(ccs_df, track_label, top_n_categories = 30) {
  wide <- ccs_df %>%
    select(CCS_ID, model, F1) %>%
    pivot_wider(names_from = model, values_from = F1) %>%
    filter(!is.na(ClinicalBERT), !is.na(SapBERT)) %>%
    mutate(delta = SapBERT - ClinicalBERT,
           direction = ifelse(delta >= 0, "Improved (SapBERT >= ClinicalBERT)", "Regressed (SapBERT < ClinicalBERT)")) %>%
    arrange(desc(abs(delta))) %>%
    slice_head(n = top_n_categories) %>%
    arrange(delta) %>%
    mutate(CCS_ID = factor(CCS_ID, levels = CCS_ID))

  n_improved <- sum(wide$delta >= 0)
  n_regressed <- sum(wide$delta < 0)

  ggplot(wide, aes(x = CCS_ID, y = delta, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Improved (SapBERT >= ClinicalBERT)" = "#55A868",
                                  "Regressed (SapBERT < ClinicalBERT)" = "#C44E52")) +
    labs(title = sprintf("Per-CCS-category F1 change: SapBERT minus ClinicalBERT (%s)", track_label),
         subtitle = sprintf("Top %d categories by |change|. %d improved, %d regressed (of %d categories with both models scored)",
                             top_n_categories, n_improved, n_regressed, nrow(wide)),
         x = "CCS category", y = "F1 (SapBERT - ClinicalBERT)", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
}

p2 <- make_ccs_breakdown_plot(ccs_10_9, track_labels["10_9"])
p3 <- make_ccs_breakdown_plot(ccs_8_9, track_labels["8_9"])

ggsave(file.path(OUT_DIR, "plot_ccs_breakdown_10_9.png"), p2, width = 8, height = 8, dpi = 150)
ggsave(file.path(OUT_DIR, "plot_ccs_breakdown_8_9.png"), p3, width = 8, height = 8, dpi = 150)

summarize_regressions <- function(ccs_df) {
  wide <- ccs_df %>%
    select(CCS_ID, model, F1) %>%
    pivot_wider(names_from = model, values_from = F1) %>%
    filter(!is.na(ClinicalBERT), !is.na(SapBERT)) %>%
    mutate(delta = SapBERT - ClinicalBERT)
  tibble(
    n_categories = nrow(wide),
    n_improved = sum(wide$delta > 0),
    n_unchanged = sum(wide$delta == 0),
    n_regressed = sum(wide$delta < 0),
    mean_delta = mean(wide$delta),
    worst_regression = min(wide$delta)
  )
}

regression_summary <- bind_rows(
  summarize_regressions(ccs_10_9) %>% mutate(track = track_labels["10_9"]),
  summarize_regressions(ccs_8_9) %>% mutate(track = track_labels["8_9"])
)
write.csv(regression_summary, file.path(OUT_DIR, "regression_summary.csv"), row.names = FALSE)

cat("\n=== Per-CCS-category regression summary ===\n")
print(regression_summary)
cat("\nCharts and summary saved to", OUT_DIR, "\n")
