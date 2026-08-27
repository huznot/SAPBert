# performance broken down over all 130 ccs categories, not just the ones with
# the biggest change. 03_ only charted the top 30 by size of change, which hides
# how uneven the rest is
#
# the caution that goes with this: 59 of the 130 categories contain a single
# icd-9 code, so a per category f1 is often computed from one or two pairs and
# jumps between 0 and 1 for reasons that have nothing to do with the model.
# the counts are kept in the output so that is visible
#
# usage: Rscript 31_category_breakdown.R
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
library(ggplot2)

ORIG_BASE    <- "data/original"
SAPBERT_BASE <- "data/sapbert"
GEN_BASE     <- "data/generated"
OUT_DIR      <- "results"

gen <- function(tag) file.path(GEN_BASE, sprintf("cosine_similarity_matrices_%%s_%s.xlsx", tag))

# no model has been settled on yet, so all three are broken down. the four
# clinicalbert text preparations are broken down too, to check whether the text
# cleaning helps in particular clinical areas or only in aggregate
MODELS <- list(
  `ClinicalBERT-original` = file.path(ORIG_BASE, "Cosine_Similarity_Matrices/cosine_similarity_matrices_%s_ClinicalBERT.xlsx"),
  `SapBERT-base`          = file.path(SAPBERT_BASE, "cosine_similarity_matrices_%s_SapBERT.xlsx"),
  `mpnet-base`            = gen("mpnet_base"),
  `ClinicalBERT-base`             = gen("clinicalbert_base"),
  `ClinicalBERT-base-nocode`      = gen("clinicalbert_base_nocode"),
  `ClinicalBERT-stopwords`        = gen("clinicalbert_stopwords"),
  `ClinicalBERT-stopwords-nocode` = gen("clinicalbert_stopwords_nocode")
)

# models that get their own per category chart, the rest appear in the tables
CHARTED <- c("ClinicalBERT-original", "SapBERT-base", "mpnet-base")

# pairs charted as a per category difference
DELTA_PAIRS <- list(
  sapbert_vs_clinicalbert = c("SapBERT-base", "ClinicalBERT-original"),
  textclean_vs_base       = c("ClinicalBERT-stopwords-nocode", "ClinicalBERT-base")
)

TRACKS <- list(
  `10_9` = list(label = "ICD-9-CM to ICD-10-CA",
                target_col_name = "ICD_10_CA", icd9_col = "ICD_9_CM_Code3", target_col = "ICD_10_CA_Code3",
                find_target_chapter_fn = find_icd10ca_chapter, chapter_alignment = chapter_alignment_10,
                manual_target_col = "ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"),
  `8_9`  = list(label = "ICD-9-CM to ICDA-8",
                target_col_name = "ICDA_8", icd9_col = "ICD_9_CM_Code", target_col = "ICDA_8_Code",
                find_target_chapter_fn = find_icda8_chapter, chapter_alignment = chapter_alignment_8,
                manual_target_col = "ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx")
)

VAL_XLSX <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
ccs_full <- read_excel(file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                       sheet = "CCS ICD-9-CM-3Level") %>%
  mutate(ICD_9_CM = as.character(ICD_9_CM))
ccs_df <- ccs_full %>% select(ICD_9_CM, CCS_ID)

# one row per ccs category, with the description and how many icd-9 codes it holds
ccs_index <- ccs_full %>%
  group_by(CCS_ID) %>%
  summarise(description = first(CCS_CATEGORY_DESCRIPTION),
            n_codes = n(), .groups = "drop")
cat(sprintf("%d ccs categories over %d icd-9 codes\n", nrow(ccs_index), nrow(ccs_full)))

for (tr in names(TRACKS)) {
  TRACKS[[tr]]$manual <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$manual_sheet)
  TRACKS[[tr]]$excl   <- read_excel(VAL_XLSX, sheet = TRACKS[[tr]]$excl_sheet)
  TRACKS[[tr]]$cooc   <- load_cooccurrence_df(file.path(ORIG_BASE, TRACKS[[tr]]$cooc_file))
}

best <- read.csv(file.path(OUT_DIR, "full_grid_best.csv"), stringsAsFactors = FALSE)

per_category <- function(track_name, sheets, thr, tn, fc) {
  tk  <- TRACKS[[track_name]]
  tcn <- tk$target_col_name
  similarity_df   <- get_similarity_scores_from_sheets(sheets, thr, tcn)
  cooccurrence_df <- get_cooccurrence_codes_from_df(tk$cooc, tn, tk$icd9_col, tk$target_col, tcn)
  merged_df <- merge_and_flag(similarity_df, cooccurrence_df, tcn,
                              tk$find_target_chapter_fn, tk$chapter_alignment)
  auto_df   <- select_rows_by_flags(merged_df, fc)
  final_valid_df <- validate_mapping(tk$manual, auto_df, ccs_df, tk$excl,
                                     manual_target_col = tk$manual_target_col,
                                     target_col_name = tcn)
  final_valid_df %>%
    group_by(CCS_ID) %>%
    summarise(n_codes_evaluated = n_distinct(`ICD-9-CM`),
              TP = sum(`True Positive`, na.rm = TRUE),
              FP = sum(`False Positive`, na.rm = TRUE),
              FN = sum(`False Negative`, na.rm = TRUE), .groups = "drop") %>%
    mutate(n_manual_pairs = TP + FN,
           precision = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
           recall    = ifelse(TP + FN > 0, TP / (TP + FN), NA_real_),
           f1        = ifelse(!is.na(precision) & !is.na(recall) & precision + recall > 0,
                              2 * precision * recall / (precision + recall),
                              ifelse(TP + FP + FN > 0, 0, NA_real_)),
           accuracy  = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_))
}

breakdown <- list()
for (track_name in names(TRACKS)) {
  for (model in names(MODELS)) {
    row <- best %>% filter(track == track_name, model == !!model)
    if (nrow(row) != 1) stop("no unique best row for ", model, " on ", track_name)
    path <- sprintf(MODELS[[model]], track_name)
    cat(sprintf("%s / %s at thr=%.3f top_n=%d flags=%d\n", track_name, model,
                row$similarity_threshold, row$top_n, row$flag_combination))
    res <- per_category(track_name, load_similarity_sheets(path),
                        row$similarity_threshold, row$top_n, row$flag_combination)
    breakdown[[length(breakdown) + 1]] <- ccs_index %>%
      left_join(res, by = "CCS_ID") %>%
      mutate(track = track_name, model = model)
  }
}
breakdown <- bind_rows(breakdown) %>%
  mutate(across(c(precision, recall, f1, accuracy), ~ round(.x, 3)))

for (track_name in names(TRACKS)) {
  out <- breakdown %>% filter(track == track_name) %>%
    select(CCS_ID, description, n_codes, model, n_codes_evaluated,
           n_manual_pairs, TP, FP, FN, precision, recall, f1, accuracy) %>%
    arrange(CCS_ID, model)
  write.csv(out, file.path(OUT_DIR, sprintf("ccs_all_categories_%s.csv", track_name)),
            row.names = FALSE)
}

# spread of the per category scores. the unweighted mean treats a category with
# one code the same as one with eighteen, so the overall f1 from the pooled
# counts is reported next to it
category_summary <- breakdown %>%
  group_by(track, model) %>%
  summarise(n_categories = n(),
            n_scored = sum(!is.na(f1)),
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
            }, 3),
            .groups = "drop")
write.csv(category_summary, file.path(OUT_DIR, "ccs_category_summary.csv"), row.names = FALSE)

# same spread split by how many codes a category holds, to show that the low
# scoring tail is mostly the small categories
size_summary <- breakdown %>%
  filter(!is.na(f1)) %>%
  mutate(size_group = cut(n_codes, breaks = c(0, 1, 2, 4, Inf),
                          labels = c("1 code", "2 codes", "3-4 codes", "5+ codes"))) %>%
  group_by(track, model, size_group) %>%
  summarise(n_categories = n(), mean_f1 = round(mean(f1), 3),
            sd_f1 = round(sd(f1), 3), .groups = "drop")
write.csv(size_summary, file.path(OUT_DIR, "ccs_category_by_size.csv"), row.names = FALSE)

# --- charts -----------------------------------------------------------

# 130 bars do not fit in one readable column, so they are ranked and dealt into
# three columns of about 44
deal_columns <- function(df, n_col = 3) {
  per_col <- ceiling(nrow(df) / n_col)
  df %>% arrange(desc(f1), CCS_ID) %>%
    mutate(rank = row_number(),
           column = sprintf("ranked %d to %d",
                            ((rank - 1) %/% per_col) * per_col + 1,
                            pmin(n(), ((rank - 1) %/% per_col + 1) * per_col)),
           label = sprintf("%d %s (%d)", CCS_ID, str_trunc(description, 30), n_codes),
           label = factor(label, levels = rev(unique(label))))
}

f1_chart <- function(track_name, model) {
  d <- breakdown %>% filter(track == track_name, model == !!model, !is.na(f1)) %>% deal_columns()
  ggplot(d, aes(x = label, y = f1, fill = f1)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_wrap(~column, scales = "free_y", nrow = 1) +
    scale_fill_gradient(low = "#C44E52", high = "#55A868", limits = c(0, 1)) +
    labs(title = sprintf("F1 by CCS category, %s, %s", TRACKS[[track_name]]$label, model),
         subtitle = sprintf("all %d categories with a score, ranked. number in brackets is how many ICD-9-CM codes the category holds",
                            nrow(d)),
         x = NULL, y = "F1", fill = "F1") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_blank())
}

pair_deltas <- function(track_name, pair) {
  breakdown %>%
    filter(track == track_name, model %in% pair) %>%
    select(CCS_ID, description, n_codes, model, f1) %>%
    tidyr::pivot_wider(names_from = model, values_from = f1) %>%
    filter(!is.na(.data[[pair[1]]]), !is.na(.data[[pair[2]]])) %>%
    mutate(delta = .data[[pair[1]]] - .data[[pair[2]]])
}

delta_chart <- function(track_name, pair) {
  a <- pair[1]; b <- pair[2]
  lv <- c(sprintf("%s better", a), "no change", sprintf("%s better", b))
  wide <- pair_deltas(track_name, pair) %>%
    mutate(direction = factor(ifelse(delta > 0, lv[1], ifelse(delta < 0, lv[3], lv[2])),
                              levels = lv)) %>%
    rename(f1 = delta) %>%
    deal_columns() %>%
    rename(delta = f1)

  counts <- table(wide$direction)
  cols <- setNames(c("#55A868", "#B0B0B0", "#C44E52"), lv)
  ggplot(wide, aes(x = label, y = delta, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_wrap(~column, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(title = sprintf("F1 change by CCS category, %s minus %s, %s",
                         a, b, TRACKS[[track_name]]$label),
         subtitle = sprintf("all %d categories. %d better, %d unchanged, %d worse",
                            nrow(wide), counts[[1]], counts[[2]], counts[[3]]),
         x = NULL, y = "F1 difference", fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_blank())
}

model_tag <- c(`ClinicalBERT-original` = "clinicalbert", `SapBERT-base` = "sapbert",
               `mpnet-base` = "mpnet")
for (track_name in names(TRACKS)) {
  for (model in CHARTED) {
    ggsave(file.path(OUT_DIR, sprintf("plot_ccs_f1_all_%s_%s.png", track_name, model_tag[[model]])),
           f1_chart(track_name, model), width = 15, height = 10, dpi = 150)
  }
  for (nm in names(DELTA_PAIRS)) {
    ggsave(file.path(OUT_DIR, sprintf("plot_ccs_delta_%s_%s.png", nm, track_name)),
           delta_chart(track_name, DELTA_PAIRS[[nm]]), width = 15, height = 10, dpi = 150)
  }
}

# how the two charted differences land across categories, in numbers
delta_summary <- bind_rows(lapply(names(DELTA_PAIRS), function(nm) {
  p <- DELTA_PAIRS[[nm]]
  bind_rows(lapply(names(TRACKS), function(tr) {
    d <- pair_deltas(tr, p)
    tibble(comparison = nm, model = p[1], reference = p[2], track = tr,
           n_categories = nrow(d),
           n_better = sum(d$delta > 0), n_same = sum(d$delta == 0),
           n_worse = sum(d$delta < 0),
           mean_delta = round(mean(d$delta), 3),
           largest_gain = round(max(d$delta), 3),
           largest_loss = round(min(d$delta), 3))
  }))
}))
write.csv(delta_summary, file.path(OUT_DIR, "ccs_delta_summary.csv"), row.names = FALSE)

# the spread itself, which is what the per category standard deviation measures
dist <- breakdown %>%
  filter(!is.na(f1), model %in% CHARTED) %>%
  mutate(track_label = ifelse(track == "10_9", TRACKS[["10_9"]]$label, TRACKS[["8_9"]]$label))
means <- dist %>% group_by(track_label, model) %>%
  summarise(mean_f1 = mean(f1), .groups = "drop")

p <- ggplot(dist, aes(x = f1, fill = model)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, colour = NA) +
  geom_vline(data = means, aes(xintercept = mean_f1, colour = model),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  facet_wrap(~track_label) +
  scale_fill_manual(values = c("ClinicalBERT-original" = "#4C72B0", "SapBERT-base" = "#DD8452",
                               "mpnet-base" = "#55A868")) +
  scale_colour_manual(values = c("ClinicalBERT-original" = "#4C72B0", "SapBERT-base" = "#DD8452",
                               "mpnet-base" = "#55A868")) +
  labs(title = "How F1 is spread across the CCS categories",
       subtitle = "dashed lines are the unweighted mean over categories",
       x = "F1 within a category", y = "Number of categories", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "plot_ccs_f1_distribution.png"), p, width = 10, height = 5, dpi = 150)

cat("\n=== per category f1, spread over categories ===\n")
print(as.data.frame(category_summary))
cat("\n=== by category size ===\n")
print(as.data.frame(size_summary))
cat("\n=== the charted differences, category by category ===\n")
print(as.data.frame(delta_summary))
cat("\nDone. Wrote ccs_all_categories_*.csv, ccs_category_summary.csv,",
    "ccs_category_by_size.csv and four charts\n")
