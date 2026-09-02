# scripts used to assume they were started from scripts/, the newer ones
# assumed the repo root, so half of them broke depending on where you ran them.
# this walks up to the root once, then every path below is root relative
local({
  d <- getwd()
  while (!file.exists(file.path(d, "icd_crosswalk.Rproj")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "icd_crosswalk.Rproj"))) stop("cant find the repo root above ", getwd())
  setwd(d)
})

# results/ used to be one directory with 110 files in it. it is grouped now, but
# the grouping is worked out here from the file name instead of at every call
# site, so scripts still pass a bare name and names built at runtime land in the
# right place too.
results_subdir <- function(name) {
  n <- basename(name)
  if (grepl("^plot_", n))                                            return("figures")
  if (grepl("^full_grid", n))                                        return("grid")
  if (grepl("^unmatched|^validation_pairs", n))                      return("unmatched")
  if (grepl("^(cv_rerank|rerank_features|base_emit|heldout_counts|precision_coverage|candidate_generation|chapter_filter)", n))
                                                                     return("rerank")
  if (grepl("^(9_codes|e13_d48|cutoff_report|absolute_threshold|identical_labels)", n))
                                                                     return("review")
  if (grepl("[.](png|pdf)$", n))                                     return("figures")
  "tables"
}

# out_path("x.csv") -> "results/tables/x.csv", creating the directory.
# a name that already carries a subdirectory is left where it asks to go.
out_path <- function(...) {
  name <- file.path(...)
  name <- sub("^results[/\\]", "", name)
  p <- if (dirname(name) != ".") file.path("results", name)
       else file.path("results", results_subdir(name), name)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
