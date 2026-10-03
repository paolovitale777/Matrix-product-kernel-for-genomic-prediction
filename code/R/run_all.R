# =============================================================================
# run_all.R - run every R analysis of Part I and Part II, then compare with Python.
#
#   Rscript code/R/run_all.R                  # everything (several hours with R's default BLAS)
#   Rscript code/R/run_all.R --quick          # one small dataset per script (about 30 minutes)
#
# Run from the repository root. Results go to code/R/outputs/. A faster BLAS (e.g. OpenBLAS,
# or R linked to Intel MKL) reduces the run time several-fold; the results do not change.
# =============================================================================
quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
R <- file.path(R.home("bin"), "Rscript")
steps <- if (!quick) list(
  c("code/R/part1/P1_T1_datasets_rank.R"),
  c("code/R/part1/P1_T2_permutation_test.R"),
  c("code/R/part1/P1_S32_psd_diagnostics.R"),
  c("code/R/part1/P1_T3_record_order_accuracy.R"),
  c("code/R/part1/P1_T3_identifier_diagnostics.R"),
  c("code/R/part1/P1_T4_S1_benchmark_5050.R"),
  c("code/R/part1/P1_T5_prospective_year.R"),
  c("code/R/part1/P1_T6_cuevas_triangles.R"),
  c("code/R/part1/P1_S2_connectivity_models.R"),
  c("code/R/part2/P2_T3_cross_covariance.R"),
  c("code/R/part2/P2_T4_identifier_neighbours.R"),
  c("code/R/part2/P2_T5_T6_identifier_split.R"),
  c("code/R/part2/P2_T7_S1_breeding_scenarios.R")
) else list(
  c("code/R/part1/P1_T1_datasets_rank.R", "--datasets", "Indica", "Japonica", "Wheat_6"),
  c("code/R/part1/P1_T2_permutation_test.R", "--datasets", "Wheat_6", "--permutations", "3"),
  c("code/R/part1/P1_S32_psd_diagnostics.R", "--datasets", "Wheat_6"),
  c("code/R/part1/P1_T3_record_order_accuracy.R", "--datasets", "Wheat_6", "--orderings", "1", "--partitions", "1"),
  c("code/R/part1/P1_T3_identifier_diagnostics.R", "--datasets", "Wheat_6", "--partitions", "1"),
  c("code/R/part1/P1_T4_S1_benchmark_5050.R", "--datasets", "Wheat_6", "--partitions", "1"),
  c("code/R/part1/P1_T5_prospective_year.R", "--datasets", "Indica", "Wheat_6"),
  c("code/R/part1/P1_T6_cuevas_triangles.R", "--datasets", "Wheat_599"),
  c("code/R/part1/P1_S2_connectivity_models.R", "--datasets", "Japonica", "--partitions", "1"),
  c("code/R/part2/P2_T3_cross_covariance.R", "--datasets", "Wheat_599"),
  c("code/R/part2/P2_T4_identifier_neighbours.R"),
  c("code/R/part2/P2_T5_T6_identifier_split.R", "--datasets", "Wheat_6"),
  c("code/R/part2/P2_T7_S1_breeding_scenarios.R", "--datasets", "Wheat_6", "Indica")
)
for (s in steps) {
  cat("\n=====", paste(s, collapse = " "), "\n"); t0 <- Sys.time()
  status <- system2(R, s)
  cat(sprintf("----- exit status %d, %.1f min\n", status, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (status != 0) stop("step failed: ", s[1])
}
system2(R, "code/R/compare_with_python.R")
