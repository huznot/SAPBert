# Are the multiple targets of one source code clustered in the target system?
#
# ICD-10 splits single ICD-9 concepts into adjacent codes. ICD-9 250 (diabetes
# mellitus) becomes E10, E11, E13, E14, which is the defined ICD-10 block for
# diabetes. Coders work from that structure: find the block, take what fits.
#
# The reranker scores every candidate independently and has no notion of a
# block, which is why it emits the strongest target and stops. If true target
# sets really are contiguous, expanding around a confident prediction should
# recover the missing ones cheaply.
#
# This measures how clustered the target sets are, then tests what a simple
# neighbour-expansion rule would recover and what it would cost in precision.

source("pipeline_lib.R")

OUT_DIR <- "../results"
ORIG <- "../data/original/ICD_Codes_Files_and_Validation_Data"

# numeric part of a target code, for measuring adjacency
codenum <- function(x) suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
codealpha <- function(x) gsub("[0-9]", "", x)

TRACKS <- list(
  `10_9` = list(sheet = "Validation_ICD9_ICD10", tcol = "ICD-10-CA"),
  `8_9`  = list(sheet = "Validaion_ICD9_ICD8",   tcol = "ICDA-8")
)

for (tr in names(TRACKS)) {
  tk <- TRACKS[[tr]]
  m <- read_excel(file.path(ORIG, "Validation_Data .xlsx"), sheet = tk$sheet) %>%
    transmute(icd9 = as.character(`ICD-9-CM`), tgt = as.character(.data[[tk$tcol]]))

  multi <- m %>% group_by(icd9) %>% filter(n() > 1)
  stats <- multi %>%
    summarise(n = n(),
              same_letter = n_distinct(codealpha(tgt)) == 1,
              span = max(codenum(tgt)) - min(codenum(tgt)),
              .groups = "drop")

  cat(sprintf("\n######## %s: %d multi-target codes ########\n", tr, nrow(stats)))
  cat(sprintf("  all targets share the same letter prefix : %.1f%%\n",
              100*mean(stats$same_letter)))
  for (s in c(2, 5, 10)) {
    cat(sprintf("  all targets within a numeric span of %-2d  : %.1f%%\n",
                s, 100*mean(stats$span <= s)))
  }
  cat(sprintf("  same letter AND span <= 10               : %.1f%%\n",
              100*mean(stats$same_letter & stats$span <= 10)))
  cat(sprintf("  median span %g, mean targets %.2f\n", median(stats$span), mean(stats$n)))
}

# --- what would neighbour expansion recover? --------------------------
# Take the auto-accepted mappings at the 95% precision point. For each, add
# any pooled candidate sharing the letter prefix and within a numeric span,
# provided the reranker gave it at least a modest score.
cat("\n\n######## neighbour expansion on the auto-accepted set ########\n")
preds <- readRDS(file.path(OUT_DIR, "cv_rerank_predictions.rds"))
ops   <- read.csv(file.path(OUT_DIR, "precision_coverage_operating_points.csv"))

for (tr in names(TRACKS)) {
  d <- preds %>% filter(track == tr)
  tau <- (ops %>% filter(track == tr, precision_target == 0.95))$tau
  n_true_all <- sum(d$y)

  seed <- d %>% filter(.p_score >= tau)
  base_tp <- sum(seed$y); base_fp <- nrow(seed) - base_tp
  cat(sprintf("\n  %s (tau %.2f)\n", tr, tau))
  cat(sprintf("    baseline: %d emitted, %d correct, precision %.3f, recall %.3f\n",
              nrow(seed), base_tp, base_tp/nrow(seed), base_tp/n_true_all))

  for (floor_p in c(0.05, 0.10, 0.20)) {
    for (span in c(2, 5)) {
      anchors <- seed %>% transmute(ICD_9_CM, a_alpha = codealpha(target),
                                    a_num = codenum(target)) %>% distinct()
      cand <- d %>% filter(.p_score >= floor_p, .p_score < tau) %>%
        mutate(c_alpha = codealpha(target), c_num = codenum(target))
      add <- cand %>% inner_join(anchors, by = "ICD_9_CM",
                                 relationship = "many-to-many") %>%
        filter(c_alpha == a_alpha, abs(c_num - a_num) <= span) %>%
        distinct(ICD_9_CM, target, y)
      tp <- base_tp + sum(add$y)
      tot <- nrow(seed) + nrow(add)
      cat(sprintf("    +expand floor %.2f span %d: %d emitted, precision %.3f, recall %.3f\n",
                  floor_p, span, tot, tp/tot, tp/n_true_all))
    }
  }
}
