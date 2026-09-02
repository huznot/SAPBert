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

# the label only matrices. the ICD code used to be pasted onto the front of
# every label before embedding, which broke the similarity scores. see
# 41_identical_labels.R for the proof
SIM <- list(
  ClinicalBERT = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_clinicalbert_base_nocode.xlsx")),
  SapBERT      = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_sapbert_base_nocode.xlsx")),
  mpnet        = list(`10_9` = file.path(GEN, "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx"),
                      `8_9`  = file.path(GEN, "cosine_similarity_matrices_8_9_mpnet_base_nocode.xlsx"))
)

best <- read.csv(out_path("full_grid_best.csv"), stringsAsFactors = FALSE)
MODEL_KEY <- c(ClinicalBERT = "ClinicalBERT-base-nocode", SapBERT = "SapBERT-base-nocode", mpnet = "mpnet-base-nocode")

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
write.csv(out, out_path("unmatched_code_handling.csv"), row.names = FALSE)

cat("\n=== codes with no correct answer, and the cost of the silence that gets them right ===\n")
print(as.data.frame(out %>% select(-track_label)), row.names = FALSE)

# counts, not shares. the two groups have different sizes, 52 codes with no
# answer against 302 with one, so a percentage axis put every red bar near 100%
# and every grey bar near zero and compared nothing
fig <- out %>%
  transmute(track_label, model,
            `mapped a code that has no answer` = wrongly_mapped,
            `emitted nothing for a code that has one` = wrongly_silent) %>%
  tidyr::pivot_longer(-c(track_label, model), names_to = "kind", values_to = "n") %>%
  mutate(kind = factor(kind, levels = c("mapped a code that has no answer",
                                        "emitted nothing for a code that has one")))

g <- ggplot(fig, aes(model, n, fill = kind)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = n), position = position_dodge(0.75),
            vjust = -0.4, size = 3.2, show.legend = FALSE) +
  facet_wrap(~track_label, scales = "free_y") +
  scale_fill_manual(values = c("mapped a code that has no answer" = "#C44E52",
                               "emitted nothing for a code that has one" = "#8C8C8C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Codes handled wrongly at each model's best setting",
       subtitle = "9 ICD-10-CA and 52 ICDA-8 codes have no correct answer, 345 and 302 have one",
       x = NULL, y = "Number of ICD-9-CM codes", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(out_path("plot_unmatched_code_handling.png"), g, width = 9.5, height = 5, dpi = 150)

cat("\nwrote results/unmatched_code_handling.csv and results/plot_unmatched_code_handling.png\n")
