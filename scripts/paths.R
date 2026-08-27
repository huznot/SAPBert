# every script runs from the repo root so all paths below are root relative.
# scripts used to assume they were started from scripts/, the newer ones assumed
# the root, so half of them broke depending on where you ran them from. this
# walks up to the root and sets the working directory once
local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(d, "icd_crosswalk.Rproj")) && dirname(d) != d) {
    d <- dirname(d)
  }
  if (!file.exists(file.path(d, "icd_crosswalk.Rproj"))) {
    stop("cant find the repo root (looked for icd_crosswalk.Rproj above ", getwd(), ")")
  }
  setwd(d)
})
