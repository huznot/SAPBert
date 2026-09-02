source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages({library(readxl); library(ggplot2)})

with_code <- "data/original/Cosine_Similarity_Matrices/cosine_similarity_matrices_10_9_ClinicalBERT.xlsx"
no_code   <- "data/generated/cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx"
labels    <- "data/original/ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx"

t9  <- read_excel(labels, sheet = "CCS ICD-9-CM-3Level")
t10 <- read_excel(labels, sheet = "ICD-10-CA-3Level")
l9  <- setNames(tolower(trimws(as.character(t9[[2]]))),  as.character(t9[[1]]))
l10 <- setNames(tolower(trimws(as.character(t10[[2]]))), as.character(t10[[1]]))
same <- l9[l9 %in% l10]

collect <- function(path, tag) {
  rows <- list()
  biggest <- 0
  for (sh in excel_sheets(path)) {
    d <- read_excel(path, sheet = sh)
    codes <- as.character(d[[1]])
    biggest <- max(biggest, max(as.numeric(as.matrix(d[, -1])), na.rm = TRUE))
    for (a in intersect(names(same), names(d))) {
      b <- names(l10)[l10 == same[a]][1]
      v <- as.numeric(d[[a]])
      rows[[length(rows) + 1]] <- data.frame(
        version = tag, icd9 = a, icd10 = b,
        label9  = trimws(as.character(t9[[2]])[as.character(t9[[1]]) == a][1]),
        label10 = trimws(as.character(t10[[2]])[as.character(t10[[1]]) == b][1]),
        score = v[codes == b],
        rank = which(codes[order(-v)] == b))
    }
  }
  r <- do.call(rbind, rows)
  attr(r, "biggest") <- biggest
  r
}

report <- function(r) {
  cat("\n---", r$version[1], "---\n")
  cat("pairs with word for word identical labels:", nrow(r), "\n")
  cat("how many score exactly 1.0:", sum(r$score > 0.9999), "of", nrow(r), "\n")
  cat("how many are the top scoring code in their own column:", sum(r$rank == 1), "of", nrow(r), "\n")
  cat("highest score anywhere in the entire matrix:", round(attr(r, "biggest"), 4), "\n")
  print(summary(r$score))
}

a <- collect(with_code, "code pasted in front of the label")
b <- collect(no_code,   "label only")

report(a)
report(b)

cat("\nwhat the model was actually given, five examples\n\n")
show <- head(a[order(a$score), ], 5)
for (i in 1:nrow(show)) {
  cat(sprintf("%s and %s   with code %.4f   label only %.4f\n",
      show$icd9[i], show$icd10[i], show$score[i],
      b$score[b$icd9 == show$icd9[i] & b$icd10 == show$icd10[i]]))
  cat(sprintf("   \"%s %s\"  vs  \"%s %s\"\n\n",
      show$icd9[i], tolower(show$label9[i]), show$icd10[i], tolower(show$label10[i])))
}

both <- rbind(a, b)
both$version <- factor(both$version, levels = c("code pasted in front of the label", "label only"))

ggplot(both, aes(score)) +
  geom_histogram(binwidth = 0.005, fill = "grey40", colour = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~version, ncol = 1) +
  coord_cartesian(xlim = c(0.84, 1.03)) +
  labs(title = paste(nrow(a), "pairs of ICD-9 and ICD-10-CA codes with word for word identical labels"),
       subtitle = "cosine similarity of each pair, same model, only the input text differs",
       x = "cosine similarity", y = "number of pairs") +
  theme_minimal()

ggsave(out_path("identical_labels.png"), width = 8, height = 6, dpi = 150)
cat("saved results/review/identical_labels.png\n")
