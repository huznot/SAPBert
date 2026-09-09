source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

# per icd-9 code, not aggregated. where does the correct icd-10-ca code actually
# sit in each model's ranking, and what kind of code is it. everything else in
# results/ collapses this to one number per condition.

ORIG <- "data/original"
GEN  <- "data/generated"
LAB  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")
VAL  <- file.path(ORIG, "ICD_Codes_Files_and_Validation_Data/Validation_Data .xlsx")

MODELS <- c(ClinicalBERT = "cosine_similarity_matrices_10_9_clinicalbert_base_nocode.xlsx",
            SapBERT      = "cosine_similarity_matrices_10_9_sapbert_base_nocode.xlsx",
            mpnet        = "cosine_similarity_matrices_10_9_mpnet_base_nocode.xlsx")

man  <- read_excel(VAL, sheet = "Validation_ICD9_ICD10")
lab9 <- read_excel(LAB, sheet = "CCS ICD-9-CM-3Level")
name9 <- setNames(as.character(lab9[[2]]), as.character(lab9[[1]]))
ccs   <- setNames(as.character(lab9$CCS_ID), as.character(lab9$ICD_9_CM))

targets <- split(as.character(man$`ICD-10-CA`), as.character(man$`ICD-9-CM`))

# a residual code is a bucket rather than a named condition. these are the ones
# the failures kept landing on, so it gets its own column instead of a guess
RESID <- "other|unspecified|not elsewhere classified|nec\\b|nos\\b|ill-defined|without specification"

rows <- list()
for (mdl in names(MODELS)) {
  cat("reading", mdl, "\n")
  sheets <- load_similarity_sheets(file.path(GEN, MODELS[[mdl]]))
  long <- purrr::map_dfr(sheets, function(df) {
    id <- names(df)[1]
    purrr::map_dfr(setdiff(names(df), id), function(c1)
      tibble::tibble(icd9 = as.character(c1), t10 = as.character(df[[id]]), s = df[[c1]]))
  })
  long <- long[!is.na(long$s), ]

  for (k in names(targets)) {
    x <- long[long$icd9 == k, ]
    if (!nrow(x)) next
    x <- x[order(-x$s), ]
    want <- targets[[k]]
    hit <- which(x$t10 %in% want)
    lab <- if (k %in% names(name9)) name9[[k]] else ""
    ch9 <- find_icd9cm_chapter(k)
    topok <- {
      d <- compute_chapter_distance(ch9, find_icd10ca_chapter(x$t10[1]), chapter_alignment_10)
      !is.na(d) && d < 1
    }
    rows[[length(rows) + 1]] <- data.frame(
      model = mdl, icd9 = k, icd9_label = lab, ccs = ccs[[k]],
      n_valid_targets = length(want),
      best_rank  = if (length(hit)) min(hit) else NA_integer_,
      best_score = if (length(hit)) x$s[min(hit)] else NA_real_,
      top_score  = x$s[1],
      margin     = if (length(hit)) x$s[1] - x$s[min(hit)] else NA_real_,
      top_survives_chapter_filter = topok,
      label_words = lengths(strsplit(trimws(lab), "\\s+")),
      residual_label = grepl(RESID, tolower(lab)),
      stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
write.csv(d, out_path("per_code_rank.csv"), row.names = FALSE)

band <- function(r) ifelse(is.na(r), "not in the column",
                    ifelse(r == 1, "1", ifelse(r <= 5, "2 to 5",
                    ifelse(r <= 20, "6 to 20", ifelse(r <= 100, "21 to 100", "over 100")))))
d$band <- factor(band(d$best_rank),
                 levels = c("1", "2 to 5", "6 to 20", "21 to 100", "over 100", "not in the column"))

show <- function(tab, title) {
  cat("\n", title, "\n", sep = "")
  pct <- round(100 * prop.table(tab, 1), 1)
  for (i in rownames(tab))
    cat(sprintf("  %-22s %s\n", i,
        paste(sprintf("%s %d (%.0f%%)", colnames(tab), tab[i, ], pct[i, ]), collapse = "   ")))
}

cat("\n", nrow(d) / length(MODELS), " icd-9 codes with a validated target, per model\n", sep = "")
show(table(d$model, d$band), "where the correct code ranks")

cat("\nsplit by whether the icd-9 label is a residual bucket\n")
for (m in names(MODELS)) {
  s <- d[d$model == m, ]
  cat("\n ", m, "\n")
  show(table(ifelse(s$residual_label, "residual", "specific"), s$band), "")
}

cat("\nhow often the top scoring code is deleted by the chapter filter\n")
t2 <- table(d$model, ifelse(d$top_survives_chapter_filter, "survives", "deleted"))
show(t2, "")

cat("\nlabel length against rank, quartiles of word count\n")
q <- quantile(d$label_words, c(0, .25, .5, .75, 1))
d$len <- cut(d$label_words, unique(q), include.lowest = TRUE)
show(table(d$len, d$band), "")

cat("\nwrote", out_path("per_code_rank.csv"), "\n")
