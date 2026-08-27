# mirrors the .rds intermediates as csv so they can be read without r.
# the rds files stay, they are what the scripts actually load, these are for
# anyone checking the numbers by hand
#
# only the out of fold predictions are exported. every held out number in the
# report is computed from that one table, so it is the one worth being able to
# open. the feature tables are 60mb as csv and 11_ rebuilds them in a minute,
# and base_emit_*.rds is an internal cache of 224 grid points, not a result
#
# usage: Rscript 34_export_rds_as_csv.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")
suppressMessages(library(dplyr))

OUT <- "results/csv_export"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FILES <- c(
  "results/cv_rerank_predictions.rds",
  "results/heldout_counts_10_9.rds",
  "results/heldout_counts_8_9.rds"
)

for (f in FILES) {
  if (!file.exists(f)) {
    cat(sprintf("  [missing %s, skipped]\n", f))
    next
  }
  obj <- readRDS(f)
  base <- sub("\\.rds$", "", basename(f))

  if (!is.data.frame(obj)) {
    obj <- tryCatch(as.data.frame(obj), error = function(e) NULL)
    if (is.null(obj)) {
      cat(sprintf("  [%s is not a table, skipped]\n", basename(f)))
      next
    }
  }
  p <- file.path(OUT, paste0(base, ".csv"))
  write.csv(obj, p, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows, %d cols)\n", p, nrow(obj), ncol(obj)))
}

cat("\nDone.\n")
