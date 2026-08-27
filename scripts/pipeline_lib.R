library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(purrr)

# low content words stripped from label text before embedding. kept in
# filler_words.json so R and generate_embeddings.py use the same list
load_filler_words <- function(path = "filler_words.json") {
  cfg <- jsonlite::fromJSON(path)
  cfg$filler_words
}

strip_filler_words <- function(text, filler_words) {
  text <- tolower(as.character(text))
  # longest phrases first so e.g. "due to" is removed before its parts
  words <- filler_words[order(-nchar(filler_words))]
  for (w in words) {
    pattern <- paste0("\\b", str_replace_all(w, " ", "\\\\s+"), "\\b")
    text <- str_replace_all(text, pattern, " ")
  }
  str_squish(text)
}

icd9cm_chapters <- tribble(
  ~chapter, ~start, ~end,
  1,  "001", "139",
  2,  "140", "239",
  3,  "240", "279",
  4,  "280", "289",
  5,  "290", "319",
  6,  "320", "389",
  7,  "390", "459",
  8,  "460", "519",
  9,  "520", "579",
  10, "580", "629",
  11, "630", "679",
  12, "680", "709",
  13, "710", "739",
  14, "740", "759",
  15, "760", "779",
  16, "780", "799",
  17, "800", "999"
)

icd10ca_chapters <- tribble(
  ~chapter, ~start, ~end,
  1,  "A00", "B99",
  2,  "C00", "D48",
  3,  "D50", "D89",
  4,  "E00", "E90",
  5,  "F00", "F99",
  6,  "G00", "G99",
  7,  "H00", "H59",
  8,  "H60", "H95",
  9,  "I00", "I99",
  10, "J00", "J99",
  11, "K00", "K93",
  12, "L00", "L99",
  13, "M00", "M99",
  14, "N00", "N99",
  15, "O00", "O99",
  16, "P00", "P96",
  17, "Q00", "Q99",
  18, "R00", "R99",
  19, "S00", "T98",
  20, "V01", "Y98",
  21, "Z00", "Z99",
  22, "U00", "U99"
)

icda8_chapters <- tribble(
  ~chapter, ~start, ~end,
  1,  "000", "136",
  2,  "140", "239",
  3,  "240", "279",
  4,  "280", "289",
  5,  "290", "315",
  6,  "320", "389",
  7,  "390", "458",
  8,  "460", "519",
  9,  "520", "577",
  10, "580", "629",
  11, "630", "678",
  12, "680", "709",
  13, "710", "738",
  14, "740", "759",
  15, "760", "779",
  16, "780", "796",
  17, "800", "999"
)

chapter_alignment_10 <- list(
  `1` = 1, `2` = 2, `3` = 4, `4` = 3, `5` = 5, `6` = c(6, 7, 8), `7` = 9,
  `8` = 10, `9` = 11, `10` = 14, `11` = 15, `12` = 12, `13` = 13, `14` = 17,
  `15` = 16, `16` = 18, `17` = 19
)

chapter_alignment_8 <- setNames(as.list(1:17), as.character(1:17))

# vectorized but still fine for length 1 input so existing scalar calls work.
# these used to run per row under rowwise, one dplyr::filter per candidate, which
# cost ~20-30s per grid point. both are lookups over a <=22 row table so compute
# once per distinct code and recycle with match. output is identical to the old
# version, verify_vectorized_equivalence.R checks it
find_chapter <- function(code, chapter_table, pad_width = 3, alpha = FALSE) {
  code <- as.character(code)
  uniq <- unique(code)
  key <- if (alpha) substr(uniq, 1, 3) else str_pad(str_split_i(uniq, "\\.", 1), pad_width, pad = "0")

  hit <- outer(key, chapter_table$start, ">=") & outer(key, chapter_table$end, "<=")
  hit[is.na(hit)] <- FALSE

  ch <- rep(NA_integer_, length(uniq))
  any_hit <- rowSums(hit) > 0
  # ties.method "first" reproduces the old hit$chapter[1]
  ch[any_hit] <- as.integer(chapter_table$chapter[max.col(hit * 1, ties.method = "first")[any_hit]])
  ch[match(code, uniq)]
}

find_icd9cm_chapter <- function(code) find_chapter(code, icd9cm_chapters, pad_width = 3, alpha = FALSE)
find_icd10ca_chapter <- function(code) find_chapter(code, icd10ca_chapters, alpha = TRUE)
find_icda8_chapter   <- function(code) find_chapter(code, icda8_chapters, pad_width = 3, alpha = FALSE)

compute_chapter_distance <- function(chapter_icd9cm, chapter_target, alignment) {
  n <- max(length(chapter_icd9cm), length(chapter_target))
  if (n == 0) return(integer(0))
  a <- rep_len(as.integer(chapter_icd9cm), n)
  b <- rep_len(as.integer(chapter_target), n)

  # at most 17 x 22 distinct source/target chapter pairs
  key <- paste(a, b, sep = "|")
  uniq <- unique(key)
  ua <- a[match(uniq, key)]
  ub <- b[match(uniq, key)]

  val <- vapply(seq_along(uniq), function(i) {
    if (is.na(ua[i]) || is.na(ub[i])) return(NA_integer_)
    allowed <- alignment[[as.character(ua[i])]]
    if (!is.null(allowed) && ub[i] %in% allowed) return(0L)
    if (abs(ua[i] - ub[i]) == 1) return(1L)
    2L
  }, integer(1))

  val[match(key, uniq)]
}

load_similarity_sheets <- function(path) {
  sheets <- excel_sheets(path)
  setNames(lapply(sheets, function(s) read_excel(path, sheet = s)), sheets)
}

get_similarity_scores_from_sheets <- function(sheets_list, threshold, target_col_name) {
  purrr::map_dfr(sheets_list, function(df) {
    id_col <- names(df)[1]
    icd9_cols <- setdiff(names(df), id_col)

    purrr::map_dfr(icd9_cols, function(col) {
      scores <- df[[col]]
      target <- df[[id_col]]
      max_score <- suppressWarnings(max(scores, na.rm = TRUE))
      if (!is.finite(max_score)) return(tibble())
      cutoff <- threshold * max_score
      keep <- which(scores >= cutoff)
      if (length(keep) == 0) return(tibble())
      tibble(
        ICD_9_CM = as.character(col),
        !!target_col_name := as.character(target[keep]),
        Similarity = scores[keep]
      )
    })
  })
}

load_cooccurrence_df <- function(path) read_excel(path)

get_cooccurrence_codes_from_df <- function(df, top_n, icd9_col, target_col, target_col_name) {
  df <- as.data.frame(df)
  names(df)[names(df) == icd9_col]  <- ".icd9"
  names(df)[names(df) == target_col] <- ".target"

  out <- df %>%
    group_by(.icd9) %>%
    slice_max(order_by = Co_Occurrence_Frequency, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      ICD_9_CM = as.character(.icd9),
      .target = as.character(.target),
      Co_Occurrence_Frequency = Co_Occurrence_Frequency
    )
  names(out)[names(out) == ".target"] <- target_col_name
  out
}

merge_and_flag <- function(similarity_df, cooccurrence_df, target_col_name,
                            find_target_chapter_fn, chapter_alignment) {
  merged <- full_join(
    similarity_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM),
                              !!target_col_name := as.character(.data[[target_col_name]])),
    cooccurrence_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM),
                                !!target_col_name := as.character(.data[[target_col_name]])),
    by = c("ICD_9_CM", target_col_name)
  )

  merged <- merged %>%
    mutate(
      chapter_icd9cm = find_icd9cm_chapter(ICD_9_CM),
      chapter_target = find_target_chapter_fn(.data[[target_col_name]]),
      chapter_distance = compute_chapter_distance(chapter_icd9cm, chapter_target, chapter_alignment)
    ) %>%
    ungroup() %>%
    filter(!is.na(chapter_distance), chapter_distance < 1)

  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

  merged <- merged %>%
    group_by(ICD_9_CM) %>%
    mutate(
      highest_similarity_flag   = as.integer(!is.na(Similarity) & Similarity == safe_max(Similarity)),
      highest_cooccurrence_flag = as.integer(!is.na(Co_Occurrence_Frequency) &
                                                Co_Occurrence_Frequency == safe_max(Co_Occurrence_Frequency)),
      exists_in_both_flag       = as.integer(!is.na(Similarity) & !is.na(Co_Occurrence_Frequency))
    ) %>%
    ungroup()

  merged
}

select_rows_by_flags <- function(merged_df, flag_combination) {
  switch(as.character(flag_combination),
    "1" = merged_df %>% filter(highest_similarity_flag == 1 | highest_cooccurrence_flag == 1),
    "2" = merged_df %>% filter(highest_similarity_flag == 1 | exists_in_both_flag == 1),
    "3" = merged_df %>% filter(highest_cooccurrence_flag == 1 | exists_in_both_flag == 1),
    "4" = merged_df %>% filter(highest_similarity_flag == 1 | highest_cooccurrence_flag == 1 | exists_in_both_flag == 1),
    stop("flag_combination must be 1-4")
  )
}

calc_precision <- function(tp, fp) if ((tp + fp) != 0) tp / (tp + fp) else 0
calc_recall    <- function(tp, fn) if ((tp + fn) != 0) tp / (tp + fn) else 0
calc_f1        <- function(p, r) if ((p + r) != 0) 2 * p * r / (p + r) else 0
calc_accuracy  <- function(tp, fp, fn) if ((tp + fp + fn) != 0) tp / (tp + fp + fn) else 0

# icd-9 codes with no reference match used to be dropped here before tp/fp/fn
# were counted, so the pipeline was never charged for mapping a code that has no
# correct answer. they are kept now, which is what the crosswalk is asked to
# handle in practice. drop_unmatched = TRUE restores the old behaviour
validate_mapping <- function(manual_df, auto_df, ccs_df, valid_excluding_df,
                              manual_target_col, target_col_name,
                              drop_unmatched = FALSE) {
  manual_df <- manual_df %>% mutate(`ICD-9-CM` = as.character(`ICD-9-CM`),
                                     !!manual_target_col := as.character(.data[[manual_target_col]]))
  auto_df <- auto_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM),
                                 !!target_col_name := as.character(.data[[target_col_name]]))
  ccs_df <- ccs_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM))

  valid_df <- full_join(
    manual_df, auto_df,
    by = c("ICD-9-CM" = "ICD_9_CM", setNames(target_col_name, manual_target_col)),
    keep = TRUE
  )

  valid_df <- valid_df %>%
    mutate(`ICD-9-CM` = coalesce(`ICD-9-CM`, ICD_9_CM)) %>%
    rename(Manual = !!manual_target_col, Auto = !!target_col_name)

  valid_df <- left_join(valid_df, ccs_df %>% select(ICD_9_CM, CCS_ID),
                         by = c("ICD-9-CM" = "ICD_9_CM"))

  if (drop_unmatched) {
    excluded <- as.character(valid_excluding_df$`ICD-9-CM`)
    valid_df <- valid_df %>% filter(!(`ICD-9-CM` %in% excluded))
  }

  final <- valid_df %>%
    mutate(
      `True Positive`  = as.integer(!is.na(Manual) & !is.na(Auto) & Manual == Auto),
      `False Negative` = as.integer(!is.na(Manual) & is.na(Auto)),
      `False Positive` = as.integer(is.na(Manual) & !is.na(Auto))
    )

  final
}

calculate_performance_metrics <- function(final_valid_df) {
  grouped <- final_valid_df %>%
    group_by(CCS_ID) %>%
    summarise(
      TP = sum(`True Positive`, na.rm = TRUE),
      FP = sum(`False Positive`, na.rm = TRUE),
      FN = sum(`False Negative`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Precision = round(mapply(calc_precision, TP, FP), 3),
      Recall    = round(mapply(calc_recall, TP, FN), 3),
      F1        = round(mapply(calc_f1, Precision, Recall), 3),
      Accuracy  = round(mapply(calc_accuracy, TP, FP, FN), 3)
    )

  tp <- sum(final_valid_df$`True Positive`, na.rm = TRUE)
  fp <- sum(final_valid_df$`False Positive`, na.rm = TRUE)
  fn <- sum(final_valid_df$`False Negative`, na.rm = TRUE)
  p  <- round(calc_precision(tp, fp), 3)
  r  <- round(calc_recall(tp, fn), 3)
  f1 <- round(calc_f1(p, r), 3)
  acc <- round(calc_accuracy(tp, fp, fn), 3)

  list(grouped = grouped,
       overall = list(overall_precision = p, overall_recall = r,
                       overall_f1_score = f1, overall_accuracy = acc))
}

run_pipeline_10_9_cached <- function(sheets_dict, cooc_df, manual_df, ccs_df, valid_excluding_df,
                                      similarity_threshold, top_n, flag_combination) {
  target_col_name <- "ICD_10_CA"
  similarity_df <- get_similarity_scores_from_sheets(sheets_dict, similarity_threshold, target_col_name)
  cooccurrence_df <- get_cooccurrence_codes_from_df(cooc_df, top_n,
                                                     icd9_col = "ICD_9_CM_Code3",
                                                     target_col = "ICD_10_CA_Code3",
                                                     target_col_name = target_col_name)
  merged_df <- merge_and_flag(similarity_df, cooccurrence_df, target_col_name,
                               find_icd10ca_chapter, chapter_alignment_10)
  auto_df <- select_rows_by_flags(merged_df, flag_combination)
  final_valid_df <- validate_mapping(manual_df, auto_df, ccs_df, valid_excluding_df,
                                      manual_target_col = "ICD-10-CA", target_col_name = target_col_name)
  metrics <- calculate_performance_metrics(final_valid_df)
  list(grouped = metrics$grouped, overall = metrics$overall, final_valid_df = final_valid_df, auto_df = auto_df)
}

run_pipeline_8_9_cached <- function(sheets_dict, cooc_df, manual_df, ccs_df, valid_excluding_df,
                                     similarity_threshold, top_n, flag_combination) {
  target_col_name <- "ICDA_8"
  similarity_df <- get_similarity_scores_from_sheets(sheets_dict, similarity_threshold, target_col_name)
  cooccurrence_df <- get_cooccurrence_codes_from_df(cooc_df, top_n,
                                                     icd9_col = "ICD_9_CM_Code",
                                                     target_col = "ICDA_8_Code",
                                                     target_col_name = target_col_name)
  merged_df <- merge_and_flag(similarity_df, cooccurrence_df, target_col_name,
                               find_icda8_chapter, chapter_alignment_8)
  auto_df <- select_rows_by_flags(merged_df, flag_combination)
  final_valid_df <- validate_mapping(manual_df, auto_df, ccs_df, valid_excluding_df,
                                      manual_target_col = "ICDA-8", target_col_name = target_col_name)
  metrics <- calculate_performance_metrics(final_valid_df)
  list(grouped = metrics$grouped, overall = metrics$overall, final_valid_df = final_valid_df, auto_df = auto_df)
}

# reverse direction, target -> icd-9-cm
# same pipeline but normalized per target instead of per icd-9 code, with an
# inverted chapter table. cosine is symmetric so the same sheets are reused

invert_alignment <- function(alignment) {
  rev <- list()
  for (src_chapter in names(alignment)) {
    for (tgt_chapter in alignment[[src_chapter]]) {
      key <- as.character(tgt_chapter)
      rev[[key]] <- unique(c(rev[[key]], as.integer(src_chapter)))
    }
  }
  rev
}

chapter_alignment_10_rev <- invert_alignment(chapter_alignment_10)
chapter_alignment_8_rev  <- invert_alignment(chapter_alignment_8)

get_similarity_scores_from_sheets_reverse <- function(sheets_list, threshold, target_col_name) {
  long_df <- purrr::map_dfr(sheets_list, function(df) {
    id_col <- names(df)[1]
    icd9_cols <- setdiff(names(df), id_col)
    purrr::map_dfr(icd9_cols, function(col) {
      tibble(ICD_9_CM = as.character(col),
             !!target_col_name := as.character(df[[id_col]]),
             Similarity = df[[col]])
    })
  })
  long_df %>%
    group_by(.data[[target_col_name]]) %>%
    mutate(max_sim = suppressWarnings(max(Similarity, na.rm = TRUE))) %>%
    filter(is.finite(max_sim), Similarity >= threshold * max_sim) %>%
    ungroup() %>%
    select(-max_sim)
}

get_cooccurrence_codes_from_df_reverse <- function(df, top_n, icd9_col, target_col, target_col_name) {
  df <- as.data.frame(df)
  names(df)[names(df) == icd9_col]  <- ".icd9"
  names(df)[names(df) == target_col] <- ".target"

  out <- df %>%
    group_by(.target) %>%
    slice_max(order_by = Co_Occurrence_Frequency, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      ICD_9_CM = as.character(.icd9),
      .target = as.character(.target),
      Co_Occurrence_Frequency = Co_Occurrence_Frequency
    )
  names(out)[names(out) == ".target"] <- target_col_name
  out
}

merge_and_flag_reverse <- function(similarity_df, cooccurrence_df, target_col_name,
                                    find_target_chapter_fn, chapter_alignment_rev) {
  merged <- full_join(
    similarity_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM),
                              !!target_col_name := as.character(.data[[target_col_name]])),
    cooccurrence_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM),
                                !!target_col_name := as.character(.data[[target_col_name]])),
    by = c("ICD_9_CM", target_col_name)
  )

  merged <- merged %>%
    mutate(
      chapter_icd9cm = find_icd9cm_chapter(ICD_9_CM),
      chapter_target = find_target_chapter_fn(.data[[target_col_name]]),
      # note the swapped argument order vs. the forward direction: distance
      # is looked up from the target chapter's allowed set, via the
      # inverted alignment table
      chapter_distance = compute_chapter_distance(chapter_target, chapter_icd9cm, chapter_alignment_rev)
    ) %>%
    ungroup() %>%
    filter(!is.na(chapter_distance), chapter_distance < 1)

  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

  merged <- merged %>%
    group_by(.data[[target_col_name]]) %>%
    mutate(
      highest_similarity_flag   = as.integer(!is.na(Similarity) & Similarity == safe_max(Similarity)),
      highest_cooccurrence_flag = as.integer(!is.na(Co_Occurrence_Frequency) &
                                                Co_Occurrence_Frequency == safe_max(Co_Occurrence_Frequency)),
      exists_in_both_flag       = as.integer(!is.na(Similarity) & !is.na(Co_Occurrence_Frequency))
    ) %>%
    ungroup()

  merged
}

validate_mapping_reverse <- function(manual_df, auto_df, ccs_df, valid_excluding_df,
                                      manual_target_col, target_col_name,
                                      drop_unmatched = FALSE) {
  manual_long <- manual_df %>%
    transmute(TargetCode = as.character(.data[[manual_target_col]]),
              ManualICD9 = as.character(`ICD-9-CM`))
  auto_long <- auto_df %>%
    transmute(TargetCode = as.character(.data[[target_col_name]]),
              AutoICD9 = as.character(ICD_9_CM))

  valid_df <- full_join(manual_long, auto_long, by = "TargetCode", relationship = "many-to-many")
  valid_df <- valid_df %>% mutate(icd9_for_ccs = coalesce(ManualICD9, AutoICD9))

  ccs_df <- ccs_df %>% mutate(ICD_9_CM = as.character(ICD_9_CM))
  valid_df <- left_join(valid_df, ccs_df %>% select(ICD_9_CM, CCS_ID),
                         by = c("icd9_for_ccs" = "ICD_9_CM"))

  if (drop_unmatched) {
    excluded <- as.character(valid_excluding_df$`ICD-9-CM`)
    valid_df <- valid_df %>% filter(!(icd9_for_ccs %in% excluded))
  }

  valid_df %>%
    mutate(
      `True Positive`  = as.integer(!is.na(ManualICD9) & !is.na(AutoICD9) & ManualICD9 == AutoICD9),
      `False Negative` = as.integer(!is.na(ManualICD9) & is.na(AutoICD9)),
      `False Positive` = as.integer(is.na(ManualICD9) & !is.na(AutoICD9))
    )
}

run_pipeline_10_9_reverse_cached <- function(sheets_dict, cooc_df, manual_df, ccs_df, valid_excluding_df,
                                              similarity_threshold, top_n, flag_combination) {
  target_col_name <- "ICD_10_CA"
  similarity_df <- get_similarity_scores_from_sheets_reverse(sheets_dict, similarity_threshold, target_col_name)
  cooccurrence_df <- get_cooccurrence_codes_from_df_reverse(cooc_df, top_n,
                                                             icd9_col = "ICD_9_CM_Code3",
                                                             target_col = "ICD_10_CA_Code3",
                                                             target_col_name = target_col_name)
  merged_df <- merge_and_flag_reverse(similarity_df, cooccurrence_df, target_col_name,
                                       find_icd10ca_chapter, chapter_alignment_10_rev)
  auto_df <- select_rows_by_flags(merged_df, flag_combination)
  final_valid_df <- validate_mapping_reverse(manual_df, auto_df, ccs_df, valid_excluding_df,
                                              manual_target_col = "ICD-10-CA", target_col_name = target_col_name)
  metrics <- calculate_performance_metrics(final_valid_df)
  list(grouped = metrics$grouped, overall = metrics$overall, final_valid_df = final_valid_df, auto_df = auto_df)
}

run_pipeline_8_9_reverse_cached <- function(sheets_dict, cooc_df, manual_df, ccs_df, valid_excluding_df,
                                             similarity_threshold, top_n, flag_combination) {
  target_col_name <- "ICDA_8"
  similarity_df <- get_similarity_scores_from_sheets_reverse(sheets_dict, similarity_threshold, target_col_name)
  cooccurrence_df <- get_cooccurrence_codes_from_df_reverse(cooc_df, top_n,
                                                             icd9_col = "ICD_9_CM_Code",
                                                             target_col = "ICDA_8_Code",
                                                             target_col_name = target_col_name)
  merged_df <- merge_and_flag_reverse(similarity_df, cooccurrence_df, target_col_name,
                                       find_icda8_chapter, chapter_alignment_8_rev)
  auto_df <- select_rows_by_flags(merged_df, flag_combination)
  final_valid_df <- validate_mapping_reverse(manual_df, auto_df, ccs_df, valid_excluding_df,
                                              manual_target_col = "ICDA-8", target_col_name = target_col_name)
  metrics <- calculate_performance_metrics(final_valid_df)
  list(grouped = metrics$grouped, overall = metrics$overall, final_valid_df = final_valid_df, auto_df = auto_df)
}

# round trip, does a target a code maps to map back to that same code
# independent of the manual crosswalk
compute_roundtrip_consistency <- function(forward_auto_df, reverse_auto_df, target_col_name) {
  fwd <- forward_auto_df %>%
    select(ICD_9_CM, !!target_col_name) %>%
    distinct()
  rev <- reverse_auto_df %>%
    select(ICD_9_CM, !!target_col_name) %>%
    distinct() %>%
    rename(ICD_9_CM_back = ICD_9_CM)

  joined <- fwd %>%
    left_join(rev, by = target_col_name, relationship = "many-to-many") %>%
    mutate(consistent = !is.na(ICD_9_CM_back) & ICD_9_CM_back == ICD_9_CM)

  n_total <- nrow(joined)
  n_consistent <- sum(joined$consistent, na.rm = TRUE)

  list(
    rate = if (n_total > 0) round(n_consistent / n_total, 3) else NA_real_,
    n_total = n_total,
    n_consistent = n_consistent,
    detail = joined
  )
}
