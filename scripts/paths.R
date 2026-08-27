# scripts used to assume they were started from scripts/, the newer ones
# assumed the repo root, so half of them broke depending on where you ran them.
# this walks up to the root once, then every path below is root relative
local({
  d <- getwd()
  while (!file.exists(file.path(d, "icd_crosswalk.Rproj")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "icd_crosswalk.Rproj"))) stop("cant find the repo root above ", getwd())
  setwd(d)
})
