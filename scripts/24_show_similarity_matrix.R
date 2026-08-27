source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})

SIM   <- "data/sapbert/cosine_similarity_matrices_10_9_SapBERT.xlsx"
LAB   <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"
VAL   <- "data/original/ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx"
OUT   <- "results"
N_SHOW <- 8

lab9  <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level") %>%
  transmute(code = as.character(ICD_9_CM), label = ICD_9_CM_LABEL)
lab10 <- read_excel(LAB, sheet = "ICD-10-CA-3Level") %>%
  transmute(code = as.character(ICD_10_CA), label = ICD_10_CA_LABEL)
truth <- read_excel(VAL, sheet = "Validation_ICD9_ICD10") %>%
  transmute(icd9 = as.character(`ICD-9-CM`), target = as.character(`ICD-10-CA`))

# the workbook is one sheet per block of icd-9 codes. rows are icd-10 codes,
# columns are icd-9 codes, each cell is the cosine similarity between the two
# label embeddings
sheets <- excel_sheets(SIM)
long <- bind_rows(lapply(sheets, function(s) {
  d <- read_excel(SIM, sheet = s)
  names(d)[1] <- "target"
  d %>% mutate(target = as.character(target)) %>%
    pivot_longer(-target, names_to = "icd9", values_to = "sim")
}))
cat(sprintf("matrix is %d icd-9 codes x %d icd-10 codes = %s cells\n\n",
            n_distinct(long$icd9), n_distinct(long$target),
            format(nrow(long), big.mark = ",")))

# a few codes that have more than one correct answer, those are the interesting
# ones to look at
multi <- truth %>% count(icd9) %>% filter(n >= 2) %>% pull(icd9)
picks <- intersect(c("250", "410", "486", "820"), multi)
if (length(picks) < 3) picks <- head(sort(multi), 4)

for (p in picks) {
  tr <- truth$target[truth$icd9 == p]
  top <- long %>% filter(icd9 == p) %>% arrange(desc(sim)) %>% head(N_SHOW) %>%
    left_join(lab10, by = c("target" = "code")) %>%
    mutate(correct = ifelse(target %in% tr, "<- correct", ""),
           rank = row_number())
  cat(sprintf("--- ICD-9 %s  %s ---\n", p, lab9$label[lab9$code == p][1]))
  cat(sprintf("    %d correct answers: %s\n", length(tr), paste(tr, collapse = ", ")))
  for (i in seq_len(nrow(top)))
    cat(sprintf("    %2d. %-5s %.4f  %-55s %s\n", top$rank[i], top$target[i],
                top$sim[i], substr(top$label[i], 1, 55), top$correct[i]))
  cat(sprintf("    correct answers found in top %d: %d of %d\n\n",
              N_SHOW, sum(top$target %in% tr), length(tr)))
}

# heatmap of one corner. take a handful of icd-9 codes and the union of their
# best candidates so the picture is small enough to read
sub9  <- head(picks, 4)
sub10 <- long %>% filter(icd9 %in% sub9) %>% group_by(icd9) %>%
  slice_max(sim, n = 6) %>% ungroup() %>% pull(target) %>% unique()
hm <- long %>% filter(icd9 %in% sub9, target %in% sub10) %>%
  left_join(truth %>% mutate(is_true = TRUE), by = c("icd9", "target")) %>%
  mutate(is_true = !is.na(is_true))

g <- ggplot(hm, aes(icd9, target, fill = sim)) +
  geom_tile(colour = "white") +
  geom_point(data = filter(hm, is_true), shape = 4, size = 3, colour = "black") +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "cosine\nsimilarity") +
  labs(x = "ICD-9-CM code", y = "ICD-10-CA candidate",
       title = "Cosine similarity between code label embeddings",
       subtitle = "x marks a pair the manual crosswalk says is correct") +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT, "plot_similarity_matrix_example.png"), g,
       width = 7, height = 6, dpi = 150)

write.csv(long %>% filter(icd9 %in% picks) %>% arrange(icd9, desc(sim)) %>%
            group_by(icd9) %>% slice_head(n = 25) %>% ungroup(),
          file.path(OUT, "similarity_matrix_example.csv"), row.names = FALSE)

cat("wrote results/plot_similarity_matrix_example.png\n")
cat("wrote results/similarity_matrix_example.csv\n")
