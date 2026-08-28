# runs the whole analysis from the committed data, in dependency order.
# about 45 minutes on 8 cores. run from anywhere:
#   Rscript run_all.R
#
# to skip the slow parameter searches and only rebuild the two stage system:
#   Rscript run_all.R --quick
#
# embeddings are not rebuilt here. data/generated/ is committed, and
# regenerating it needs python and a gpu, see scripts/run_all_embeddings.sh

source(if (file.exists("scripts/paths.R")) "scripts/paths.R" else "paths.R")

QUICK <- "--quick" %in% commandArgs(trailingOnly = TRUE)

RSCRIPT <- file.path(R.home("bin"), "Rscript")

run <- function(script, args = character()) {
  cat(sprintf("\n===== %s %s =====\n", script, paste(args, collapse = " ")))
  t0 <- Sys.time()
  st <- system2(RSCRIPT, c(file.path("scripts", script), args))
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (st != 0) stop(script, " failed with status ", st)
  cat(sprintf("----- %s done in %.1f min\n", script, mins))
}

# the four step pipeline, searched over the full grid
if (!QUICK) {
  run("07_full_grid_comparison.R")   # slowest step, all 12 conditions
  run("08_assemble_full_grid.R")     # summary tables + report figure 1
  run("09_error_analysis.R")         # where correct pairs are lost, section 4
  run("10_candidate_generation_study.R")
}

run("11_rerank_features.R")          # candidate pool + features
run("12_cv_rerank.R", "10_9")
run("12_cv_rerank.R", "8_9")
run("12b_merge_cv_results.R")
run("13_precision_coverage.R")       # confidence thresholds and triage

# follow ups that depend on the cv predictions
run("30_variability.R")              # bootstrap sd, section 12
run("31_category_breakdown.R")       # all 130 ccs categories, section 13

# these name their output after the tracks they were given, and 18_ reads the
# per track files back, so pass one track at a time
if (!QUICK) {
  for (tr in c("10_9", "8_9")) {
    run("16_ablation.R", tr)
    run("17_retrieval_sensitivity.R", tr)
    run("18_learning_curve.R", tr)
    run("19_category_holdout.R", tr)
  }
  run("23_code_prefix_test.R")
  run("29_portability.R")
}

# reporting only, no recomputation
run("20_top1_accuracy.R")
run("21_error_by_code_type.R")
run("22_target_block_structure.R")
run("25_frequency_distributions.R")
run("26_stopword_choice.R")
run("28_stopwords_and_codes.R")
run("35_unmatched_codes.R")          # codes with no correct answer, section 2
run("36_unmatched_descriptives.R")   # those codes described on their own, section 11
run("27_show_results.R")             # every headline number

cat("\nall done. results/ is rebuilt.\n")
