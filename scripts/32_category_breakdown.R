# performance broken down over all 130 ccs categories, for the two stage system
# and for the original rules, both measured on held out codes
#
# the caution that goes with this: most categories contain one or two icd-9
# codes, so a per category f1 is often computed from a handful of pairs and
# jumps between 0 and 1 for reasons that have nothing to do with the system.
# the counts are kept in the output so that is visible
#
# needs heldout_counts_<track>.rds from 31_variability.R
# usage: Rscript 32_category_breakdown.R
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
library(ggplot2)

ORIG_BASE <- "data/original"
OUT_DIR   <- "results"
SYSTEMS   <- c("two stage system", "original rules")
TRACKS    <- list(`10_9` = "ICD-9-CM to ICD-10-CA", `8_9` = "ICD-9-CM to ICDA-8")

ccs_full <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                       sheet = "CCS ICD-9-CM-3Level") %>%
  mutate(ICD_9_CM = as.character(ICD_9_CM))
ccs_index <- ccs_full %>%
  group_by(CCS_ID) %>%
  summarise(description = first(CCS_CATEGORY_DESCRIPTION), n_codes = n(), .groups = "drop")
cat(sprintf("%d ccs categories over %d icd-9 codes\n", nrow(ccs_index), nrow(ccs_full)))

breakdown <- list()
for (tr in names(TRACKS)) {
  path <- file.path(OUT_DIR, sprintf("heldout_counts_%s.rds", tr))
  if (!file.exists(path)) stop("missing ", path, ", run 31_variability.R first")
  counts <- readRDS(path)
  for (sys in SYSTEMS) {
    per_cat <- counts[[sys]] %>%
      mutate(ICD_9_CM = as.character(code)) %>%
      left_join(ccs_full %>% select(ICD_9_CM, CCS_ID), by = "ICD_9_CM") %>%
      group_by(CCS_ID) %>%
      summarise(n_codes_evaluated = n_distinct(ICD_9_CM),
                TP = sum(TP), FP = sum(FP), FN = sum(FN), .groups = "drop") %>%
      mutate(n_true_pairs = TP + FN,
             precision = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
             recall    = ifelse(TP + FN > 0, TP / (TP + FN), NA_real_),
             f1 = ifelse(!is.na(precision) & !is.na(recall) & precision + recall > 0,
                         2 * precision * recall / (precision + recall),
                         ifelse(TP + FP + FN > 0, 0, NA_real_)),
             accuracy = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_))
    breakdown[[length(breakdown) + 1]] <- ccs_index %>%
      left_join(per_cat, by = "CCS_ID") %>%
      mutate(track = tr, system = sys)
  }
}
breakdown <- bind_rows(breakdown) %>%
  mutate(across(c(precision, recall, f1, accuracy), ~ round(.x, 3)))

for (tr in names(TRACKS)) {
  breakdown %>% filter(track == tr) %>%
    select(CCS_ID, description, n_codes, system, n_codes_evaluated, n_true_pairs,
           TP, FP, FN, precision, recall, f1, accuracy) %>%
    arrange(CCS_ID, system) %>%
    write.csv(file.path(OUT_DIR, sprintf("ccs_all_categories_%s.csv", tr)), row.names = FALSE)
}

# the unweighted mean treats a category holding one code the same as one holding
# eighteen, so the pooled f1 over all codes is reported next to it
category_summary <- breakdown %>%
  group_by(track, system) %>%
  summarise(n_categories = n(), n_scored = sum(!is.na(f1)),
            mean_f1 = round(mean(f1, na.rm = TRUE), 3),
            sd_f1 = round(sd(f1, na.rm = TRUE), 3),
            median_f1 = round(median(f1, na.rm = TRUE), 3),
            q1_f1 = round(quantile(f1, 0.25, na.rm = TRUE, names = FALSE), 3),
            q3_f1 = round(quantile(f1, 0.75, na.rm = TRUE, names = FALSE), 3),
            n_perfect = sum(f1 == 1, na.rm = TRUE),
            n_zero = sum(f1 == 0, na.rm = TRUE),
            pooled_f1 = round({
              p <- sum(TP, na.rm = TRUE) / sum(TP + FP, na.rm = TRUE)
              r <- sum(TP, na.rm = TRUE) / sum(TP + FN, na.rm = TRUE)
              2 * p * r / (p + r)
            }, 3), .groups = "drop")
write.csv(category_summary, file.path(OUT_DIR, "ccs_category_summary.csv"), row.names = FALSE)

size_summary <- breakdown %>%
  filter(!is.na(f1)) %>%
  mutate(size_group = cut(n_codes, breaks = c(0, 1, 2, 4, Inf),
                          labels = c("1 code", "2 codes", "3-4 codes", "5+ codes"))) %>%
  group_by(track, system, size_group) %>%
  summarise(n_categories = n(), mean_f1 = round(mean(f1), 3),
            sd_f1 = round(sd(f1), 3), .groups = "drop")
write.csv(size_summary, file.path(OUT_DIR, "ccs_category_by_size.csv"), row.names = FALSE)

# --- charts -----------------------------------------------------------

# 130 bars do not fit in one readable column, so they are ranked and dealt into
# three columns of about 44
deal_columns <- function(df, n_col = 3) {
  per_col <- ceiling(nrow(df) / n_col)
  df %>% arrange(desc(value), CCS_ID) %>%
    mutate(rank = row_number(),
           column = sprintf("ranked %d to %d",
                            ((rank - 1) %/% per_col) * per_col + 1,
                            pmin(n(), ((rank - 1) %/% per_col + 1) * per_col)),
           label = sprintf("%d %s (%d)", CCS_ID, str_trunc(description, 30), n_codes),
           label = factor(label, levels = rev(unique(label))))
}

f1_chart <- function(tr, sys) {
  d <- breakdown %>% filter(track == tr, system == sys, !is.na(f1)) %>%
    rename(value = f1) %>% deal_columns()
  ggplot(d, aes(x = label, y = value, fill = value)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_wrap(~column, scales = "free_y", nrow = 1) +
    scale_fill_gradient(low = "#C44E52", high = "#55A868", limits = c(0, 1)) +
    labs(title = sprintf("F1 by CCS category, %s, %s", TRACKS[[tr]], sys),
         subtitle = sprintf("held out. all %d categories with a score, ranked. number in brackets is how many ICD-9-CM codes the category holds",
                            nrow(d)),
         x = NULL, y = "F1", fill = "F1") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_blank())
}

pair_deltas <- function(tr) {
  breakdown %>% filter(track == tr) %>%
    select(CCS_ID, description, n_codes, system, f1) %>%
    tidyr::pivot_wider(names_from = system, values_from = f1) %>%
    filter(!is.na(`two stage system`), !is.na(`original rules`)) %>%
    mutate(delta = `two stage system` - `original rules`)
}

delta_chart <- function(tr) {
  lv <- c("two stage system better", "no change", "original rules better")
  wide <- pair_deltas(tr) %>%
    mutate(direction = factor(ifelse(delta > 0, lv[1], ifelse(delta < 0, lv[3], lv[2])),
                              levels = lv)) %>%
    rename(value = delta) %>% deal_columns() %>% rename(delta = value)
  counts <- table(wide$direction)
  ggplot(wide, aes(x = label, y = delta, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_wrap(~column, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = setNames(c("#55A868", "#B0B0B0", "#C44E52"), lv), drop = FALSE) +
    labs(title = sprintf("F1 change by CCS category, two stage system minus original rules, %s",
                         TRACKS[[tr]]),
         subtitle = sprintf("held out. all %d categories. %d better, %d unchanged, %d worse",
                            nrow(wide), counts[[1]], counts[[2]], counts[[3]]),
         x = NULL, y = "F1 difference", fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_blank())
}

sys_tag <- c(`two stage system` = "twostage", `original rules` = "rules")
for (tr in names(TRACKS)) {
  for (sys in SYSTEMS)
    ggsave(file.path(OUT_DIR, sprintf("plot_ccs_f1_all_%s_%s.png", tr, sys_tag[[sys]])),
           f1_chart(tr, sys), width = 15, height = 10, dpi = 150)
  ggsave(file.path(OUT_DIR, sprintf("plot_ccs_delta_twostage_vs_rules_%s.png", tr)),
         delta_chart(tr), width = 15, height = 10, dpi = 150)
}

delta_summary <- bind_rows(lapply(names(TRACKS), function(tr) {
  d <- pair_deltas(tr)
  tibble(track = tr, n_categories = nrow(d),
         n_better = sum(d$delta > 0), n_same = sum(d$delta == 0), n_worse = sum(d$delta < 0),
         mean_delta = round(mean(d$delta), 3),
         largest_gain = round(max(d$delta), 3), largest_loss = round(min(d$delta), 3))
}))
write.csv(delta_summary, file.path(OUT_DIR, "ccs_delta_summary.csv"), row.names = FALSE)

dist <- breakdown %>% filter(!is.na(f1)) %>%
  mutate(track_label = ifelse(track == "10_9", TRACKS[["10_9"]], TRACKS[["8_9"]]))
means <- dist %>% group_by(track_label, system) %>%
  summarise(mean_f1 = mean(f1), .groups = "drop")
cols <- c("two stage system" = "#DD8452", "original rules" = "#4C72B0")
p <- ggplot(dist, aes(x = f1, fill = system)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, colour = NA) +
  geom_vline(data = means, aes(xintercept = mean_f1, colour = system),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  facet_wrap(~track_label) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(title = "How F1 is spread across the CCS categories",
       subtitle = "held out. dashed lines are the unweighted mean over categories",
       x = "F1 within a category", y = "Number of categories", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "plot_ccs_f1_distribution.png"), p, width = 10, height = 5, dpi = 150)

cat("\n=== per category f1, spread over categories ===\n")
print(as.data.frame(category_summary))
cat("\n=== by category size ===\n")
print(as.data.frame(size_summary))
cat("\n=== two stage system against the rules, category by category ===\n")
print(as.data.frame(delta_summary))
cat("\nDone. Wrote ccs_all_categories_*.csv, ccs_category_summary.csv,",
    "ccs_category_by_size.csv, ccs_delta_summary.csv and five charts\n")
