# checks what stripping filler words does to the label text before embedding
#
# the risk with any stopword list on icd labels is two codes with different
# meanings ending up with identical text. "with complication" and "without
# complication" are separate codes, if the list removes with/without they
# collapse and no model can tell them apart afterwards
#
# compares the hand written list against standard english stopword lists

source("pipeline_lib.R")

ORIG_BASE <- "../data/original"
LABELS <- file.path(ORIG_BASE, "ICD_Codes_Files_and_Validation_Data/ICD_Codes_Labels.xlsx")

sheets <- list(
  `ICD-9-CM`  = list(sheet = "CCS ICD-9-CM-3Level", code = "ICD_9_CM",  lab = "ICD_9_CM_LABEL"),
  `ICD-10-CA` = list(sheet = "ICD-10-CA-3Level",    code = "ICD_10_CA", lab = "ICD_10_CA_LABEL"),
  `ICDA-8`    = list(sheet = "ICDA-8-3Level",       code = "ICDA_8",    lab = "ICDA_8_LABEL")
)

# the hand written list currently in use
ours <- load_filler_words("filler_words.json")

# standard lists. snowball/smart come from the stopwords package if installed,
# otherwise fall back to copies so this runs anywhere
nltk_like <- c("i","me","my","myself","we","our","ours","ourselves","you","your","yours",
  "he","him","his","she","her","hers","it","its","they","them","their","what","which","who",
  "whom","this","that","these","those","am","is","are","was","were","be","been","being",
  "have","has","had","having","do","does","did","doing","a","an","the","and","but","if","or",
  "because","as","until","while","of","at","by","for","with","about","against","between",
  "into","through","during","before","after","above","below","to","from","up","down","in",
  "out","on","off","over","under","again","further","then","once","here","there","when",
  "where","why","how","all","any","both","each","few","more","most","other","some","such",
  "no","nor","not","only","own","same","so","than","too","very","s","t","can","will","just",
  "don","should","now")

LISTS <- list(`filler_words.json` = ours, `NLTK english` = nltk_like)
if (requireNamespace("stopwords", quietly = TRUE)) {
  for (s in c("snowball", "smart", "stopwords-iso")) {
    LISTS[[s]] <- tryCatch(stopwords::stopwords("en", source = s), error = function(e) NULL)
  }
  LISTS <- LISTS[!vapply(LISTS, is.null, logical(1))]
}

# words that change clinical meaning and should survive any list
CLINICAL <- c("with","without","not","no","nor","other","due","in","of","against","between",
              "during","before","after","above","below","under","over","further")

rows <- list()
for (ln in names(LISTS)) {
  words <- LISTS[[ln]]
  cat(sprintf("\n===== %s (%d words) =====\n", ln, length(words)))
  cat(sprintf("  clinically loaded words it removes: %s\n",
              paste(intersect(CLINICAL, words), collapse = ", ")))

  for (sn in names(sheets)) {
    sp <- sheets[[sn]]
    d <- read_excel(LABELS, sheet = sp$sheet)
    code <- as.character(d[[sp$code]]); lab <- as.character(d[[sp$lab]])
    stripped <- strip_filler_words(lab, words)

    empty <- sum(nchar(stripped) == 0)
    before <- length(unique(tolower(lab)))
    after  <- length(unique(stripped))
    # codes whose stripped text now matches another code's stripped text
    dup <- sum(duplicated(stripped) | duplicated(stripped, fromLast = TRUE)) -
           (sum(duplicated(tolower(lab)) | duplicated(tolower(lab), fromLast = TRUE)))

    cat(sprintf("  %-10s %4d codes | %3d labels emptied | %3d new collisions\n",
                sn, nrow(d), empty, dup))

    if (dup > 0) {
      coll <- tibble(code = code, lab = tolower(lab), stripped = stripped) %>%
        group_by(stripped) %>% filter(n() > 1, n_distinct(lab) > 1) %>%
        summarise(codes = paste(code, collapse = "/"),
                  labels = paste(unique(lab), collapse = "  ||  "), .groups = "drop") %>%
        head(3)
      for (i in seq_len(nrow(coll)))
        cat(sprintf("      %s -> \"%s\"\n         %s\n",
                    coll$codes[i], coll$stripped[i], coll$labels[i]))
    }

    rows[[length(rows)+1]] <- tibble(
      list = ln, n_words = length(words), system = sn, n_codes = nrow(d),
      labels_emptied = empty, new_collisions = dup,
      unique_before = before, unique_after = after,
      clinical_words_removed = length(intersect(CLINICAL, words)))
  }
}

out <- bind_rows(rows)
write.csv(out, "../results/stopword_collision_check.csv", row.names = FALSE)
cat("\n\n===== summary =====\n")
print(as.data.frame(out %>% group_by(list, n_words) %>%
  summarise(labels_emptied = sum(labels_emptied),
            new_collisions = sum(new_collisions),
            clinical_words_removed = first(clinical_words_removed),
            .groups = "drop") %>% arrange(new_collisions)))
cat("\nWritten to results/stopword_collision_check.csv\n")

# the list in use must not merge any two codes. run this after changing
# filler_words.json
ours_coll <- sum(out$new_collisions[out$list == "filler_words.json"])
if (ours_coll > 0) {
  cat(sprintf("\nFAIL: filler_words.json merges %d code(s). Fix the list.\n", ours_coll))
  quit(status = 1)
}
cat("\nPASS: filler_words.json merges no codes.\n")
