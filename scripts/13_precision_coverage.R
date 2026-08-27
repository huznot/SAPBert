# precision at coverage, and a review queue
#
# f1 assumes the system has to answer everything. it doesnt, some codes have no
# counterpart at all and what matters in practice is knowing which mappings are
# safe to accept without checking
#
# reports precision vs coverage, the threshold that hits a target precision and
# how much of the crosswalk it builds there, and an auto-accept / review / no
# candidate split. runs off the out of fold predictions from 12_
source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
source("scripts/pipeline_lib.R")

OUT_DIR <- "results"
preds <- readRDS(file.path(OUT_DIR, "cv_rerank_predictions.rds"))

prf <- function(tp, fp, fn) {
  p <- if (tp + fp > 0) tp/(tp+fp) else NA_real_
  r <- if (tp + fn > 0) tp/(tp+fn) else 0
  c(precision = p, recall = r, f1 = if (!is.na(p) && p + r > 0) 2*p*r/(p+r) else 0)
}

curve_rows <- list(); op_rows <- list(); triage_rows <- list()

for (tr in unique(preds$track)) {
  d <- preds %>% filter(track == tr)
  n_codes <- n_distinct(d$ICD_9_CM)
  n_true  <- sum(d$y)   # true pairs reachable in the pool (the ceiling)
  cat(sprintf("\n######## track %s: %d candidates, %d codes, %d reachable true pairs ########\n",
              tr, nrow(d), n_codes, n_true))

  # precision / coverage curve over the emission threshold
  for (tau in c(seq(0.01, 0.95, by = 0.01), 0.97, 0.99)) {
    em <- d %>% filter(.p_score >= tau)
    tp <- sum(em$y); fp <- nrow(em) - tp
    m <- prf(tp, fp, n_true - tp)
    curve_rows[[length(curve_rows)+1]] <- tibble(
      track = tr, tau = tau,
      n_emitted = nrow(em),
      codes_covered = n_distinct(em$ICD_9_CM),
      code_coverage = n_distinct(em$ICD_9_CM) / n_codes,
      precision = m["precision"], recall = m["recall"], f1 = m["f1"])
  }

  cur <- bind_rows(curve_rows) %>% filter(track == tr)

  # operating points at target precision. the practical question is, if we insist
  # on being right 95% of the time how much can we build automatically
  for (target in c(0.80, 0.90, 0.95, 0.99)) {
    ok <- cur %>% filter(!is.na(precision), precision >= target)
    if (!nrow(ok)) {
      cat(sprintf("  precision >= %.2f : NOT ACHIEVABLE at any threshold\n", target)); next
    }
    # among thresholds meeting the precision target, take the one that emits most
    b <- ok %>% slice_max(n_emitted, n = 1, with_ties = FALSE)
    cat(sprintf("  precision >= %.2f : tau %.2f -> %d mappings, %.1f%% of codes covered, actual precision %.3f, recall %.3f\n",
                target, b$tau, b$n_emitted, 100*b$code_coverage, b$precision, b$recall))
    op_rows[[length(op_rows)+1]] <- b %>% mutate(precision_target = target)
  }

  # triage, what a coder would actually see
  # auto-accept   : candidates above the 95%-precision threshold
  # review        : the code has candidates, but none confident enough
  # no candidate  : nothing plausible retrieved at all
  ok95 <- cur %>% filter(!is.na(precision), precision >= 0.95)
  tau95 <- if (nrow(ok95)) (ok95 %>% slice_max(n_emitted, n = 1, with_ties = FALSE))$tau else NA_real_
  if (!is.na(tau95)) {
    auto <- d %>% filter(.p_score >= tau95)
    auto_codes <- unique(auto$ICD_9_CM)
    per_code <- d %>% group_by(ICD_9_CM) %>%
      summarise(best_p = max(.p_score), n_true_here = sum(y), .groups = "drop") %>%
      mutate(bucket = case_when(
        ICD_9_CM %in% auto_codes ~ "auto-accept",
        best_p >= 0.10           ~ "review (candidates, low confidence)",
        TRUE                     ~ "review (no plausible candidate)"))
    tb <- per_code %>% count(bucket) %>% mutate(pct = round(100*n/sum(n), 1))
    cat(sprintf("\n  triage at precision>=0.95 (tau %.2f):\n", tau95))
    for (i in seq_len(nrow(tb))) cat(sprintf("    %-38s %3d codes (%.1f%%)\n", tb$bucket[i], tb$n[i], tb$pct[i]))
    tp_auto <- sum(auto$y)
    cat(sprintf("    auto-accepted mappings: %d, of which correct: %d (%.1f%%)\n",
                nrow(auto), tp_auto, 100*tp_auto/nrow(auto)))
    cat(sprintf("    => %.1f%% of the manual crosswalk built automatically at %.1f%% precision\n",
                100*tp_auto/n_true, 100*tp_auto/nrow(auto)))
    triage_rows[[length(triage_rows)+1]] <- tb %>% mutate(track = tr, tau = tau95)
  }
}

curve <- bind_rows(curve_rows)
write.csv(curve, file.path(OUT_DIR, "precision_coverage_curve.csv"), row.names = FALSE)
write.csv(bind_rows(op_rows), file.path(OUT_DIR, "precision_coverage_operating_points.csv"), row.names = FALSE)
write.csv(bind_rows(triage_rows), file.path(OUT_DIR, "precision_coverage_triage.csv"), row.names = FALSE)

# plot
png(file.path(OUT_DIR, "plot_precision_coverage.png"), width = 1100, height = 500, res = 110)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
for (tr in unique(curve$track)) {
  cc <- curve %>% filter(track == tr, !is.na(precision)) %>% arrange(recall)
  ttl <- if (tr == "10_9") "ICD-9 -> ICD-10-CA" else "ICD-9 -> ICDA-8"
  plot(cc$recall, cc$precision, type = "l", lwd = 2, col = "#2b6cb0",
       xlab = "Recall (share of reachable true pairs emitted)", ylab = "Precision",
       main = ttl, ylim = c(0, 1), xlim = c(0, 1))
  abline(h = c(0.90, 0.95), lty = 3, col = c("#888888", "#c05621"))
  text(0.02, 0.96, "95% precision", cex = 0.7, col = "#c05621", adj = 0)
  grid(col = "#eeeeee")
}
invisible(dev.off())

cat("\nWritten: results/precision_coverage_curve.csv, _operating_points.csv, _triage.csv, plot_precision_coverage.png\n")
