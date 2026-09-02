source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(library(openxlsx))

XLS <- out_path("9_codes_review_UPDATED.xlsx")

patch <- list(
  "515" = "J84 is right. I went through the ICD-10-CA three character list and J84, other interstitial pulmonary diseases, is the only code in the respiratory chapter that covers lung fibrosis. The ones near it do not fit. J70 is respiratory conditions due to external agents, J80 is adult respiratory distress syndrome, J98 is other respiratory disorders. The CMS general equivalence mappings send 515 to J84.10, pulmonary fibrosis unspecified, or J84.89, other specified interstitial pulmonary diseases. Both are J84 at three characters, so it does not matter which one you take. That is a US mapping, ICD-10-CM not ICD-10-CA, but the two share the same three character rubrics here. Cosine returns nothing for this code, J84 comes from co-occurrence.",
  "327" = "G47 is right. ICD-10-CA splits sleep disorders by cause and both halves are in our three character list. G47, other sleep disorders, is in the nervous system chapter. F51, nonorganic sleep disorders, is in the mental health chapter. ICD-9 327 is the organic one, it pairs with 307.4 for the nonorganic. So 327 goes to G47 and 307.4 goes to F51. The CMS general equivalence mappings agree, every 327 subcode goes to a G47 subcode, for example 327.00 to G47.01, 327.20 to G47.30 and 327.35 to G47.25. Both branches return G47."
)

wb <- loadWorkbook(XLS)
sh <- names(wb)[1]
d  <- read.xlsx(XLS, sheet = 1)
codes <- as.character(d[[1]])
notes_col <- which(names(d) == "Notes")

for (k in names(patch)) {
  r <- which(codes == k)
  if (!length(r)) { cat("skipped", k, "not found\n"); next }
  writeData(wb, sh, patch[[k]], startRow = r + 1, startCol = notes_col, colNames = FALSE)
  cat("patched", k, "at row", r + 1, "\n")
}

saveWorkbook(wb, XLS, overwrite = TRUE)
cat("saved", XLS, "\n")
