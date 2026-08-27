# how the pipeline handles icd-9 codes that have no correct answer
#
# 9 icd-9-cm codes have no icd-10-ca reference match and 52 have no icda-8
# match. they used to be dropped before scoring, so a wrong mapping on them
# cost nothing. they are scored now.
#
# f1 catches that when it goes wrong, because a mapping on one of them is a
# false positive. what f1 cannot do is credit getting them right. a code with
# no correct answer has no true positive to earn, so staying silent scores the
# same as never being tested. this reports that half directly.
#
# it also counts the codes that do have answers and got nothing, because the
# same silence causes both and reporting only the first is misleading
#
# usage: Rscript 35_unmatched_codes.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")
suppressMessages(library(ggplot2))

ORIG <- "data/original"; SAP <- "data/sapbert"; GEN <- "data/generated"
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")
OUT  <- "results"

ccs <- read_excel(file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"),
                  sheet = "CCS ICD-9-CM-3Level") %>% select(ICD_9_CM, CCS_ID)

TRACKS <- list(
  `10_9` = list(tcn = "ICD_10_CA", label = "ICD-9-CM to ICD-10-CA",
                manual_sheet = "Validation_ICD9_ICD10", excl_sheet = "Validation_ICD9_ICD10_Excld",
                cooc_file = "Co_occurrence/icd_10_9_co_occurrence_3c.xlsx",
                manual_target_col = "ICD-10-CA", run = run_pipeline_10_9_cached),
  `8_9`  = list(tcn = "ICDA_8", label = "ICD-9-CM to ICDA-8",
                manual_sheet = "Validaion_ICD9_ICD8", excl_sheet = "Validation_ICD9_ICD8_Excld",
                cooc_file = "Co_occurrence/icd_8_9_co_occurrence_3d.xlsx",
                manual_target_col = "ICDA-8", run = run_pipeline_8_9_cached)
)

SIM <- list(
  ClinicalBERT = list(`10_9` = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"),
                      `8_9`  = file.path(ORIG, "Cosine_Similarity_Matrices/cosine_similarity_matrices_8_9_ClinicalBERT.xlsx")),
  SapBERT      = list(`10_9` = file.path(SAP, "cosine_similarity_matrices_10_9_SapBERT.xlsx"),
                      `8_9`  = file.path(SAP, "cosine_similarity_matrices_8_9_SapBERT.xlsx")),
  mpnet        = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_mpnet_base.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_mpnet_base.xlsx"))
)

best <- read.csv(file.path(OUT, "full_grid_best.csv"), stringsAsFactors = FALSE)
MODEL_KEY <- c(ClinicalBERT = "ClinicalBERT-original", SapBERT = "SapBERT-base", mpnet = "mpnet-base")

rows <- list()
for (tr in names(TRACKS)) {
  tk   <- TRACKS[[tr]]
  man  <- read_excel(VAL, sheet = tk$manual_sheet)
  exc  <- read_excel(VAL, sheet = tk$excl_sheet)
  cooc <- load_cooccurrence_df(file.path(ORIG, tk$cooc_file))
  unmatched <- unique(as.character(exc$`ICD-9-CM`))
  matched   <- unique(as.character(man$`ICD-9-CM`))

  for (mdl in names(SIM)) {
    b <- best[best$track == tr & best$model == MODEL_KEY[mdl], ]
    if (!nrow(b)) { cat(sprintf("  [no grid row for %s / %s, skipped]\n", tr, mdl)); next }

    res <- tk$run(load_similarity_sheets(SIM[[mdl]][[tr]]), cooc, man, ccs, exc,
                  similarity_threshold = b$similarity_threshold[1],
                  top_n = b$top_n[1], flag_combination = b$flag_combination[1])
    emitted_for <- unique(as.character(res$auto_df$ICD_9_CM))

    rows[[length(rows) + 1]] <- tibble(
      track = tr, track_label = tk$label, model = mdl,
      n_unmatched = length(unmatched),
      correctly_silent = length(setdiff(unmatched, emitted_for)),
      wrongly_mapped   = length(intersect(unmatched, emitted_for)),
      n_matched = length(matched),
      wrongly_silent = length(setdiff(matched, emitted_for)),
      correct_rejection_rate = round(length(setdiff(unmatched, emitted_for)) / length(unmatched), 3))
  }
}

out <- bind_rows(rows)
write.csv(out, file.path(OUT, "unmatched_code_handling.csv"), row.names = FALSE)

cat("\n=== codes with no correct answer, and the cost of the silence that gets them right ===\n")
print(as.data.frame(out %>% select(-track_label)), row.names = FALSE)

# the tradeoff, as two bars per model. staying quiet is what gets the unmatched
# codes right and it is also what loses real mappings, so both are plotted
fig <- out %>%
  transmute(track_label, model,
            `codes with no answer, wrongly mapped` = wrongly_mapped / n_unmatched,
            `codes with an answer, nothing emitted` = wrongly_silent / n_matched) %>%
  tidyr::pivot_longer(-c(track_label, model), names_to = "kind", values_to = "share")

lab <- out %>%
  transmute(track_label, model,
            `codes with no answer, wrongly mapped` = sprintf("%d of %d", wrongly_mapped, n_unmatched),
            `codes with an answer, nothing emitted` = sprintf("%d of %d", wrongly_silent, n_matched)) %>%
  tidyr::pivot_longer(-c(track_label, model), names_to = "kind", values_to = "txt")

fig <- fig %>% left_join(lab, by = c("track_label", "model", "kind"))

g <- ggplot(fig, aes(model, share, fill = kind)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = txt), position = position_dodge(0.75),
            vjust = -0.4, size = 3, show.legend = FALSE) +
  facet_wrap(~track_label) +
  scale_fill_manual(values = c("codes with no answer, wrongly mapped" = "#C44E52",
                               "codes with an answer, nothing emitted" = "#8C8C8C")) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Two ways of being wrong about whether a code has an answer",
       subtitle = "the same caution avoids the red bar and causes the grey one",
       x = NULL, y = "Share of codes", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "plot_unmatched_code_handling.png"), g, width = 9.5, height = 5, dpi = 150)

cat("\nwrote results/unmatched_code_handling.csv and results/plot_unmatched_code_handling.png\n")
