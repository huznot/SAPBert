source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(source("scripts/pipeline_lib.R"))
suppressMessages({library(openxlsx); library(readxl)})

XLS  <- out_path("9_codes_review_UPDATED.xlsx")
LAB  <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"
COOC <- "data/original/Co_occurrence/icd_10_9_co_occurrence_3c.xlsx"
GEN  <- "data/generated"

arms <- list(
  list(tab = "ClinicalBERT", file = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx", thr = 0.995, top_n = 25),
  list(tab = "SapBERT",      file = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",      thr = 0.95,  top_n = 25),
  list(tab = "mpnet",        file = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx",        thr = 0.95,  top_n = 30)
)

NOMAP <- "no correct ICD-10-CA mapping"

base <- data.frame(
  `ICD-9-CM` = c("339","338","175","515","327","249","239","445","209"),
  `ICD-9-CM label` = c("other headache syndromes","pain not elsewhere classified",
    "malignant neoplasm of male breast","postinflammatory pulmonary fibrosis",
    "organic sleep disorders","secondary diabetes mellitus",
    "neoplasms of unspecified nature","atheroembolism","neuroendocrine tumors"),
  `Correct ICD-10-CA` = c("G44","none","C50","J84","G47","none","none","none","none"),
  `ICD-10-CA label` = c("Other headache syndromes", NOMAP, "Malignant neoplasm of breast",
    "Other interstitial pulmonary diseases","Other sleep disorders",
    NOMAP, NOMAP, NOMAP, NOMAP),
  Verdict = c("validation data wrong","validation data correct","validation data wrong",
    "validation data wrong","validation data wrong",
    "validation data correct","validation data correct","validation data correct",
    "validation data correct"),
  Notes = c(
    "G44 is the same label as the ICD-9 code, word for word. The pipeline gets it now. It used to miss it because the ICD code was pasted onto the front of the label before the model saw it, which pushed G44 down to 0.9050 and 11th place. With the label alone it scores 1.0 and comes first.",
    "There is no three character code for this. R52 is the top similarity match out of all 2038 but it sits in chapter 18 and this ICD-9 code is only allowed chapters 6, 7 and 8, so the chapter filter drops it. Anything the pipeline returns here is a false positive.",
    "The pipeline gets C50. C79 comes in from co-occurrence, not from meaning.",
    "J84 is right. I went through the ICD-10-CA three character list and J84, other interstitial pulmonary diseases, is the only code in the respiratory chapter that covers lung fibrosis. The ones near it do not fit. J70 is respiratory conditions due to external agents, J80 is adult respiratory distress syndrome, J98 is other respiratory disorders. The CMS general equivalence mappings send 515 to J84.10, pulmonary fibrosis unspecified, or J84.89, other specified interstitial pulmonary diseases. Both are J84 at three characters, so it does not matter which one you take. That is a US mapping, ICD-10-CM not ICD-10-CA, but the two share the same three character rubrics here. Cosine returns nothing for this code, J84 comes from co-occurrence.",
    "G47 is right. ICD-10-CA splits sleep disorders by cause and both halves are in our three character list. G47, other sleep disorders, is in the nervous system chapter. F51, nonorganic sleep disorders, is in the mental health chapter. ICD-9 327 is the organic one, it pairs with 307.4 for the nonorganic. So 327 goes to G47 and 307.4 goes to F51. The CMS general equivalence mappings agree, every 327 subcode goes to a G47 subcode, for example 327.00 to G47.01, 327.20 to G47.30 and 327.35 to G47.25. Both branches return G47.",
    "There is no three character code for this, so E13 is not the answer. All three models still rank E13 near the top with the label only matrices, third on SapBERT and fourth on mpnet, so the models are reading it fine. The mapping just does not exist at three characters.",
    "There is no three character code for this, so D48 is not the answer. All three models rank D48 second, so the pipeline is not failing to see it, there is nothing correct to return.",
    "There is no three character code for atheroembolism. I74 is the nearest thing to it but nearest is not correct, so the pipeline returning I74 is a false positive.",
    "ICD-10-CA has no neuroendocrine code so there is nothing to map to. The pipeline cannot answer no match. Co-occurrence always hands back its top N, so something comes out even when nothing fits."),
  check.names = FALSE, stringsAsFactors = FALSE)

nine <- base$`ICD-9-CM`

t10 <- read_excel(LAB, sheet = "ICD-10-CA-3Level")
lookup <- setNames(as.character(t10[[2]]), as.character(t10[[1]]))
cooc <- load_cooccurrence_df(COOC)

cell <- function(cs) {
  if (!length(cs)) return("nothing")
  paste(vapply(cs, function(c1) paste(c1,
    if (c1 %in% names(lookup)) lookup[[c1]] else "not in the ICD-10-CA code list"),
    character(1)), collapse = " / ")
}

column_for <- function(path, icd9) {
  for (sh in excel_sheets(path)) {
    if (icd9 %in% names(read_excel(path, sheet = sh, n_max = 0))) {
      d <- read_excel(path, sheet = sh)
      return(data.frame(code = as.character(d[[1]]), score = as.numeric(d[[icd9]]),
                        stringsAsFactors = FALSE))
    }
  }
  NULL
}

why_missing <- function(path, thr, top_n, icd9, want, sim_c, co_c) {
  emitted <- unique(c(sim_c, co_c))
  if (want == "none") {
    if (!length(emitted)) return("there is no correct ICD-10-CA code and the pipeline returns nothing, which is the right answer")
    return(sprintf("there is no correct ICD-10-CA code, so every code returned here is a false positive. cosine returns %s and co-occurrence returns %s. the pipeline has no way to answer no match",
                   if (length(sim_c)) paste(sim_c, collapse = " and ") else "nothing",
                   if (length(co_c)) paste(co_c, collapse = " and ") else "nothing"))
  }
  if (want %in% emitted) return("")
  x <- column_for(path, icd9)
  co <- get_cooccurrence_codes_from_df(cooc, top_n, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
  co_codes <- as.character(co$ICD_10_CA[as.character(co$ICD_9_CM) == icd9])
  ch9 <- find_icd9cm_chapter(icd9)
  if (!want %in% x$code) return(paste0(want, " is not in the ICD-10-CA code list the pipeline searches"))
  s <- x$score[x$code == want]
  dist <- compute_chapter_distance(ch9, find_icd10ca_chapter(want), chapter_alignment_10)
  if (is.na(dist) || dist >= 1)
    return(sprintf("%s is removed by the chapter filter. it sits in ICD-10-CA chapter %s and ICD-9 chapter %s is only allowed to match chapters one step away or closer",
                   want, find_icd10ca_chapter(want), ch9))
  cut <- thr * max(x$score)
  if (s < cut)
    return(sprintf("%s scores %.4f and the cutoff is %.4f, which is %s of the highest score in the column, so it is dropped before anything else%s",
                   want, s, cut, thr,
                   if (want %in% co_codes) ", although it is in the top co-occurrence codes" else ", and it is not in the top co-occurrence codes either"))
  sprintf("%s clears the cutoff at %.4f but %s scores higher, and only the single highest scoring code is handed on%s",
          want, s, if (length(sim_c)) sim_c[1] else "another code",
          if (want %in% co_codes) ", even though it is in the top co-occurrence codes" else ", and it is not in the top co-occurrence codes either")
}

wb <- createWorkbook()
head_style <- createStyle(textDecoration = "bold", valign = "bottom")
body_style <- createStyle(valign = "top", wrapText = TRUE)

for (a in arms) {
  cat("running", a$tab, "\n")
  path <- file.path(GEN, a$file)
  sim <- get_similarity_scores_from_sheets(load_similarity_sheets(path), a$thr, "ICD_10_CA")
  co  <- get_cooccurrence_codes_from_df(cooc, a$top_n, "ICD_9_CM_Code3", "ICD_10_CA_Code3", "ICD_10_CA")
  m <- merge_and_flag(sim, co, "ICD_10_CA", find_icd10ca_chapter, chapter_alignment_10)
  m$ICD_9_CM <- as.character(m$ICD_9_CM)

  sim_l <- lapply(nine, function(k) sort(unique(as.character(m$ICD_10_CA[m$ICD_9_CM == k & m$highest_similarity_flag == 1]))))
  co_l  <- lapply(nine, function(k) sort(unique(as.character(m$ICD_10_CA[m$ICD_9_CM == k & m$highest_cooccurrence_flag == 1]))))

  d <- base
  d$`Cosine output`        <- vapply(sim_l, cell, character(1))
  d$`Co-occurrence output` <- vapply(co_l,  cell, character(1))
  d$`Why the correct code is missing` <- vapply(seq_along(nine), function(i)
    why_missing(path, a$thr, a$top_n, nine[i], base$`Correct ICD-10-CA`[i], sim_l[[i]], co_l[[i]]),
    character(1))

  cols <- c("ICD-9-CM","ICD-9-CM label","Correct ICD-10-CA","ICD-10-CA label",
            "Cosine output","Co-occurrence output","Why the correct code is missing")
  widths <- c(10, 28, 16, 32, 42, 42, 70)
  if (a$tab == "ClinicalBERT") { cols <- c(cols, "Verdict", "Notes"); widths <- c(widths, 26, 80) }
  d <- d[, cols]

  addWorksheet(wb, a$tab)
  writeData(wb, a$tab, d, headerStyle = head_style)
  addStyle(wb, a$tab, body_style, rows = 2:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE, stack = TRUE)
  setColWidths(wb, a$tab, cols = seq_along(widths), widths = widths)
  freezePane(wb, a$tab, firstActiveRow = 2)
}

saveWorkbook(wb, XLS, overwrite = TRUE)
cat("wrote", XLS, "\n")
