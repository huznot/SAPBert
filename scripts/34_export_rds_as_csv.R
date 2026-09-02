# the out of fold predictions as csv, so the held out numbers can be checked
# without r. the rds files stay, they are what the scripts load
#
# only the predictions are exported. the feature tables are 60mb as csv and 11_
# rebuilds them in a minute, and base_emit_*.rds is a cache of 224 grid points
#
# usage: Rscript 34_export_rds_as_csv.R

source(if (file.exists("paths.R")) "paths.R" else "scripts/paths.R")

OUT <- "results/rerank/csv_export"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

files <- c(out_path("cv_rerank_predictions.rds"),
           out_path("heldout_counts_10_9.rds"),
           out_path("heldout_counts_8_9.rds"))

for (f in files) {
  if (!file.exists(f)) {
    cat(sprintf("  [missing %s, skipped]\n", f))
    next
  }
  d <- readRDS(f)
  p <- file.path(OUT, sub("\\.rds$", ".csv", basename(f)))
  write.csv(d, p, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows, %d cols)\n", p, nrow(d), ncol(d)))
}

cat("\nDone.\n")
