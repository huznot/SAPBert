# descriptive analysis of the icd-9 codes that have no reference match
#
# 9 icd-9-cm codes have no icd-10-ca match and 52 have no icda-8 match. this
# describes the two signals the pipeline runs on, cosine similarity and
# co-occurrence frequency, for those codes on their own, with the codes that do
# have a match as a comparison group.
#
# nothing here changes the pipeline. it is a description of what the two
# signals look like for these codes, which is the step before deciding whether
# any rule could tell them apart.
#
# usage: Rscript 36_unmatched_descriptives.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})
options(width = 170)

ORIG <- "data/original"; SAP <- "data/sapbert"; GEN <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
OUT  <- "results"

TRACKS <- list(
  `10_9` = list(label = "ICD-9-CM to ICD-10-CA", target = "ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                cooc_icd9 = "ICD_9_CM_Code3", cooc_target = "ICD_10_CA_Code3",
                target_labels = "ICD-10-CA-3Level"),
  `8_9`  = list(label = "ICD-9-CM to ICDA-8", target = "ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                cooc_icd9 = "ICD_9_CM_Code", cooc_target = "ICDA_8_Code",
                target_labels = "ICDA-8-3Level")
)

SIM <- list(
  ClinicalBERT = list(`10_9` = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
                      `8_9`  = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")),
  SapBERT      = list(`10_9` = file.path(SAP, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                      `8_9`  = file.path(SAP, "cosine_similarity_matrices_8_9_SapBERT.xlsx")),
  mpnet        = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_mpnet_base.xlsx"))
)

# the similarity workbook is one sheet per ccs category, all target codes as
# rows and that category's icd-9 codes as columns. stacking every sheet gives
# one row per icd-9 code and target pair, which is what a distribution needs
long_similarity <- function(path) {
  bind_rows(lapply(excel_sheets(path), function(s) {
    d <- read_excel(path, sheet = s)
    names(d)[1] <- "target"
    d %>% mutate(target = as.character(target)) %>%
      pivot_longer(-target, names_to = "icd9", values_to = "similarity")
  }))
}

describe <- function(x) {
  tibble(n = length(x), mean = mean(x), sd = sd(x), min = min(x),
         q1 = quantile(x, .25), median = median(x), q3 = quantile(x, .75), max = max(x))
}

r3 <- function(d) d %>% mutate(across(where(is.numeric), ~round(.x, 3)))

# all 52 icda-8 unmatched codes are absent from the co-occurrence file, so a
# whole group can be empty and min/max would return Inf
q <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE))

icd9_labels <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>%
  transmute(icd9 = as.character(ICD_9_CM), icd9_label = ICD_9_CM_LABEL,
            ccs_id = CCS_ID, ccs_category = CCS_CATEGORY_DESCRIPTION)

sim_by_code <- list(); cooc_by_code <- list(); top_pairs <- list(); group_rows <- list()

for (tr in names(TRACKS)) {
  tk  <- TRACKS[[tr]]
  unmatched <- unique(as.character(read_excel(VAL, sheet = tk$excl_sheet)$`ICD-9-CM`))
  matched   <- unique(as.character(read_excel(VAL, sheet = tk$manual_sheet)$`ICD-9-CM`))
  tgt_lab   <- read_excel(LAB, sheet = tk$target_labels) %>% setNames(c("target", "target_label")) %>%
    mutate(target = as.character(target))

  cat(sprintf("\n\n=================== %s ===================\n", tk$label))
  cat(sprintf("%d icd-9 codes with no %s match, %d with at least one\n",
              length(unmatched), tk$target, length(matched)))

  group_of <- function(x) ifelse(x %in% unmatched, "no reference match", "has a reference match")

  # ---- cosine similarity ----
  for (mdl in names(SIM)) {
    L <- long_similarity(SIM[[mdl]][[tr]]) %>% mutate(group = group_of(icd9))

    per_code <- L %>% group_by(icd9, group) %>%
      reframe(describe(similarity)) %>%
      left_join(L %>% group_by(icd9) %>% slice_max(similarity, n = 1, with_ties = FALSE) %>%
                  ungroup() %>% transmute(icd9, nearest_target = target), by = "icd9") %>%
      left_join(tgt_lab, by = c("nearest_target" = "target")) %>%
      rename(nearest_target_label = target_label) %>%
      mutate(track = tr, track_label = tk$label, model = mdl, .before = 1)
    sim_by_code[[paste(tr, mdl)]] <- per_code

    # the plainest form of the separation question. a cutoff on the maximum
    # similarity set just above the highest of the codes with no answer clears
    # all of them, and cutoff_cost is how many codes with an answer it takes
    ceiling_no_match <- max(per_code$max[per_code$group == "no reference match"])

    grp <- per_code %>% group_by(group) %>%
      summarise(codes = n(),
                mean_similarity = mean(mean),
                mean_of_max = mean(max), min_of_max = min(max),
                q1_of_max = quantile(max, .25), median_of_max = median(max),
                q3_of_max = quantile(max, .75), max_of_max = max(max), .groups = "drop") %>%
      mutate(cutoff_clearing_all = ceiling_no_match,
             cutoff_cost = sum(per_code$max[per_code$group == "has a reference match"] <= ceiling_no_match)) %>%
      mutate(track = tr, track_label = tk$label, model = mdl, signal = "cosine similarity", .before = 1)
    group_rows[[paste(tr, mdl)]] <- grp

    if (mdl == "ClinicalBERT") {
      cat("\n--- cosine similarity against every target code, ClinicalBERT ---\n")
      cat("each icd-9 code is scored against all", per_code$n[1], tk$target, "codes\n\n")
      cat("averaged over the codes in each group:\n")
      print(as.data.frame(r3(grp %>% select(group, codes, mean_similarity, mean_of_max,
                                            min_of_max, q1_of_max, median_of_max,
                                            q3_of_max, max_of_max))), row.names = FALSE)
      cat("\nper code, the ones with no reference match:\n")
      print(as.data.frame(r3(per_code %>% filter(group == "no reference match") %>%
        arrange(desc(max)) %>%
        select(icd9, mean, sd, min, q1, median, q3, max, nearest_target, nearest_target_label))),
        row.names = FALSE)

      # how far the two groups overlap, against the quartiles of the codes that
      # do have a match. icda-8 reuses the icd-9 numbering, so a code with no
      # reference match can still have a target carrying its own number
      no_m  <- per_code %>% filter(group == "no reference match")
      has_m <- per_code %>% filter(group == "has a reference match")
      cat(sprintf("\n%d of the %d codes with no match reach a higher maximum than the first quartile\n",
                  sum(no_m$max > quantile(has_m$max, .25)), nrow(no_m)))
      cat(sprintf("of the codes that do have one, and %d beat that group's median\n",
                  sum(no_m$max > median(has_m$max))))
      cat(sprintf("a target code carrying the same three digit number exists for %d of them, and is\n",
                  sum(no_m$icd9 %in% tgt_lab$target)))
      cat(sprintf("the nearest target for %d\n", sum(no_m$nearest_target == no_m$icd9)))
      cat(sprintf("\na cutoff on the maximum above %.3f is the lowest that clears all %d codes with\n",
                  max(no_m$max), nrow(no_m)))
      cat(sprintf("no match, and it also removes %d of the %d codes that have one\n",
                  sum(has_m$max <= max(no_m$max)), nrow(has_m)))
    }
  }

  # ---- co-occurrence ----
  cooc <- read_excel(file.path(ORIG, tk$cooc_file))
  names(cooc)[names(cooc) == tk$cooc_icd9]   <- "icd9"
  names(cooc)[names(cooc) == tk$cooc_target] <- "target"
  cooc <- cooc %>% transmute(icd9 = as.character(icd9), target = as.character(target),
                             freq = Co_Occurrence_Frequency)

  all_codes <- tibble(icd9 = union(unmatched, matched)) %>% mutate(group = group_of(icd9))
  per_cooc <- all_codes %>%
    left_join(cooc %>% group_by(icd9) %>%
                reframe(partners = n(), describe(freq) %>% select(-n), total = sum(freq)),
              by = "icd9") %>%
    mutate(track = tr, track_label = tk$label, .before = 1)
  cooc_by_code[[tr]] <- per_cooc

  cat("\n--- co-occurrence frequency ---\n")
  cat("a code is in the co-occurrence file only if it shares a record with at\n")
  cat("least one target code, so a missing code has no empirical data at all\n\n")
  cooc_grp <- per_cooc %>% group_by(group) %>%
    summarise(codes = n(), no_cooccurrence_data = sum(is.na(partners)),
              median_partners = q(partners, .5),
              min_of_max = q(max, 0), q1_of_max = q(max, .25),
              median_of_max = q(max, .5), q3_of_max = q(max, .75),
              max_of_max = q(max, 1), .groups = "drop")
  print(as.data.frame(cooc_grp), row.names = FALSE)
  group_rows[[paste(tr, "cooc")]] <- cooc_grp %>%
    mutate(track = tr, track_label = tk$label, model = NA_character_,
           signal = "co-occurrence frequency", .before = 1)

  present <- per_cooc %>% filter(group == "no reference match", !is.na(partners))
  if (nrow(present)) {
    cat("\nper code, the ones with no reference match, by highest single frequency:\n")
    print(as.data.frame(present %>% arrange(desc(max)) %>%
      select(icd9, partners, min, q1, median, q3, max, total)), row.names = FALSE)

    tp <- cooc %>% filter(icd9 %in% present$icd9) %>%
      group_by(icd9) %>% slice_max(freq, n = 3, with_ties = FALSE) %>% ungroup() %>%
      left_join(tgt_lab, by = "target") %>%
      left_join(icd9_labels %>% select(icd9, icd9_label), by = "icd9") %>%
      mutate(track = tr, track_label = tk$label, .before = 1)
    top_pairs[[tr]] <- tp
    cat("\nthe three most frequent partners of each, this is the check for an\n")
    cat("unexpectedly high co-occurrence:\n")
    print(as.data.frame(tp %>% select(icd9, icd9_label, target, target_label, freq)), row.names = FALSE)
  } else {
    cat("\nnone of the", length(unmatched), "codes with no reference match appear in the\n")
    cat("co-occurrence file at all\n")
  }
}

sim_by_code  <- bind_rows(sim_by_code) %>% left_join(icd9_labels, by = "icd9")
cooc_by_code <- bind_rows(cooc_by_code) %>% left_join(icd9_labels, by = "icd9")
group_summary <- bind_rows(group_rows)

write.csv(r3(sim_by_code),  file.path(OUT, "unmatched_similarity_by_code.csv"), row.names = FALSE)
write.csv(cooc_by_code,     file.path(OUT, "unmatched_cooccurrence_by_code.csv"), row.names = FALSE)
write.csv(r3(group_summary), file.path(OUT, "unmatched_descriptives_summary.csv"), row.names = FALSE)
if (length(top_pairs))
  write.csv(bind_rows(top_pairs), file.path(OUT, "unmatched_top_cooccurrence_pairs.csv"), row.names = FALSE)

# ---- charts ----
# groups differ in size, 9 against 345 and 52 against 302, so counts on a
# shared axis would say nothing. one point per code with a box behind it shows
# every code and still puts the two spreads on the same scale
GRP_COL <- c(`no reference match` = "#C44E52", `has a reference match` = "#8C8C8C")
lvl <- function(d) d %>% mutate(group = factor(group, levels = names(GRP_COL)))

g1 <- ggplot(lvl(sim_by_code %>% filter(model == "ClinicalBERT")), aes(group, max, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey35") +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 1.1, show.legend = FALSE) +
  facet_wrap(~track_label, scales = "free_x") +
  coord_flip() +
  scale_colour_manual(values = GRP_COL) +
  labs(title = "Highest cosine similarity reached by each ICD-9-CM code",
       subtitle = "one point per code, scored against every target code, ClinicalBERT",
       x = NULL, y = "maximum cosine similarity") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "plot_unmatched_similarity_max.png"), g1, width = 9.5, height = 4, dpi = 150)

# the whole distribution for each code with no match, which is what was asked
# for. the box is min to max with the quartiles inside
spread <- sim_by_code %>% filter(model == "ClinicalBERT", group == "no reference match") %>%
  arrange(track_label, max) %>% mutate(pos = factor(row_number()))
g2 <- ggplot(spread, aes(pos)) +
  geom_linerange(aes(ymin = min, ymax = max), colour = "#BFBFBF", linewidth = 0.8) +
  geom_linerange(aes(ymin = q1, ymax = q3), colour = "#C44E52", linewidth = 2.2) +
  geom_point(aes(y = median), colour = "white", size = 0.9) +
  geom_point(aes(y = max), colour = "#31374A", size = 1.4) +
  coord_flip() +
  # free on both scales so each panel drops the other track's codes, the two
  # lists overlap on five codes and a shared axis interleaves them wrongly
  facet_wrap(~track_label, scales = "free") +
  scale_x_discrete(labels = setNames(spread$icd9, as.character(spread$pos))) +
  labs(title = "Cosine similarity of each unmatched ICD-9-CM code against every target code",
       subtitle = "grey line min to max, red bar the quartiles, white dot the median, dark dot the maximum (ClinicalBERT)",
       x = "ICD-9-CM code", y = "cosine similarity") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "plot_unmatched_similarity_spread.png"), g2, width = 10, height = 7.5, dpi = 150)

cooc_plot <- cooc_by_code %>% filter(!is.na(max))
missing_lab <- cooc_by_code %>% group_by(track_label, group) %>%
  summarise(missing = sum(is.na(max)), n = n(), .groups = "drop") %>%
  mutate(txt = sprintf("%d of %d codes have no co-occurrence data", missing, n))
g3 <- ggplot(lvl(cooc_plot), aes(group, max, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey35") +
  geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 1.1, show.legend = FALSE) +
  geom_text(data = lvl(missing_lab), aes(x = group, y = 8, label = txt),
            hjust = 0, nudge_x = 0.44, size = 2.9, colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~track_label) +
  coord_flip() +
  scale_y_log10() +
  scale_colour_manual(values = GRP_COL) +
  labs(title = "Frequency of each ICD-9-CM code's most co-occurring target",
       subtitle = "one point per code, log scale, codes absent from the co-occurrence file cannot be plotted",
       x = NULL, y = "co-occurrence count of the most frequent partner") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "plot_unmatched_cooccurrence.png"), g3, width = 9.5, height = 5, dpi = 150)

cat("\n\nwrote results/unmatched_similarity_by_code.csv\n")
cat("wrote results/unmatched_cooccurrence_by_code.csv\n")
cat("wrote results/unmatched_descriptives_summary.csv\n")
cat("wrote results/unmatched_top_cooccurrence_pairs.csv\n")
cat("wrote results/plot_unmatched_similarity_max.png\n")
cat("wrote results/plot_unmatched_similarity_spread.png\n")
cat("wrote results/plot_unmatched_cooccurrence.png\n")
