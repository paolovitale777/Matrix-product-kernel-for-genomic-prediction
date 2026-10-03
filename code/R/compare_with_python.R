# =============================================================================
# compare_with_python.R - put every table cell computed by the R scripts next to the value
# obtained from the published Python outputs, and report the differences.
#
# Run after any of the R scripts (only the tables whose R output exists are compared):
#     Rscript code/R/compare_with_python.R
# Output: code/R/outputs/comparison_R_vs_Python.csv  (table, dataset, model/column, R, Python,
#         difference) and a summary of the largest absolute difference per table.
# Table cells are means over the same partitions / splits in both languages.
# =============================================================================
source("code/R/00_setup.R")
PY1 <- file.path("code", "python", "part1", "outputs")
PY2 <- file.path("code", "python", "part2", "outputs")
rout <- function(f) { p <- file.path(OUT_DIR, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
pyout <- function(dir, pattern) {
  fs <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(fs)) return(NULL)
  xs <- lapply(fs, function(f) { x <- read.csv(f, stringsAsFactors = FALSE); x$.file <- basename(f); x })
  cols <- unique(unlist(lapply(xs, names)))                  # files can differ in their columns
  do.call(rbind, lapply(xs, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
}
cmp <- list()
add <- function(table, key, r, py) {
  common <- intersect(names(r), names(py))
  if (!length(common)) return()
  cmp[[length(cmp) + 1]] <<- data.frame(table = table, cell = common, R = unname(r[common]), Python = unname(py[common]),
                                         difference = unname(r[common] - py[common]), stringsAsFactors = FALSE)
}
cellmeans <- function(df, by, value) { a <- aggregate(df[[value]], df[by], mean); setNames(a$x, do.call(paste, c(a[by], sep = " | "))) }

# ---- Part I Table 2 ------------------------------------------------------------
r <- rout("P1_T2_permutation_summary.csv"); py <- read.csv(file.path(PY1, "order_diagnostics_summary.csv"))
if (!is.null(r)) for (v in c("C0_delta_median", "P0_delta_median")) {
  pv <- sub("0_delta", "_delta", paste0("original_", v))
  add("P1 Table 2", NULL, setNames(r[[v]], paste(r$dataset, v, sep = " | ")), setNames(py[[pv]], paste(py$dataset, v, sep = " | ")))
}
# ---- Part I Table 3 (identifier diagnostics) -----------------------------------
r <- rout("P1_T3_identifier_diagnostics.csv")
if (!is.null(r)) {
  py <- read.csv(file.path(PY1, "id_order_metrics.csv")); py$partition <- as.integer(sub("partition_", "", py$partition))
  py <- merge(py, unique(r[c("dataset", "partition")]))
  for (v in c("COR_within", "PM10_within")) add(paste("P1 Table 3 identifier", v), NULL, cellmeans(r, c("dataset", "model"), v), cellmeans(py, c("dataset", "model"), v))
}
# ---- Part I Table 3 (record orders) ----------------------------------------------
r <- rout("P1_T3_record_order_accuracy.csv")
if (!is.null(r)) {
  py <- merge(read.csv(file.path(PY1, "order_accuracy_sensitivity.csv")), unique(r[c("dataset", "ordering", "partition")]))
  r$o <- ifelse(r$ordering == 0, "file order", "random orders"); py$o <- ifelse(py$ordering == 0, "file order", "random orders")
  add("P1 Table 3 record orders", NULL, cellmeans(r, c("dataset", "model", "o"), "COR_within"), cellmeans(py, c("dataset", "model", "o"), "COR_within"))
}
# ---- Part I Table 4 / S1 ---------------------------------------------------------
r <- rout("P1_T4_S1_benchmark_5050.csv")
if (!is.null(r)) {
  py <- read.csv(file.path(PY1, "v4_ulte_metrics.csv")); py <- py[py$subset == "all", ]
  py$partition <- as.integer(sub("partition_", "", py$split)); py <- merge(py, unique(r[c("dataset", "partition")]))
  for (v in c("COR_within", "PM10_within", "COR_pooled")) add(paste("P1 Table 4/S1", v), NULL, cellmeans(r, c("dataset", "model"), v), cellmeans(py, c("dataset", "model"), v))
}
# ---- Part I Table 5 ----------------------------------------------------------------
r <- rout("P1_T5_prospective_year.csv")
if (!is.null(r)) {
  py <- rbind(transform(read.csv(file.path(PY1, "prospective_gy_metrics.csv")), ec = "obs"),
              transform(read.csv(file.path(PY1, "prospective_gy_metrics_ecprev.csv")), ec = "prev"))
  py <- py[!(py$dataset %in% c("Indica", "Japonica") & py$ec == "obs"), ]
  add("P1 Table 5", NULL, cellmeans(r, c("dataset", "subset", "model"), "COR_within"), cellmeans(py, c("dataset", "subset", "model"), "COR_within"))
}
# ---- Part I Table 6 ----------------------------------------------------------------
r <- rout("P1_T6_cuevas_triangles.csv")
if (!is.null(r)) {
  py <- pyout(PY1, "^cuevas_files_.*\\.csv$"); py <- merge(py, unique(r[c("dataset", "env")]))
  em <- function(x) { a <- aggregate(COR ~ dataset + env + group, x, mean); cellmeans(a, c("dataset", "group"), "COR") }
  add("P1 Table 6", NULL, em(r), em(py))
}
# ---- Part I Table S2 -----------------------------------------------------------------
r <- rout("P1_S2_connectivity_models.csv")
if (!is.null(r)) {
  py <- read.csv(file.path(PY1, "recover_gain_metrics.csv")); py$partition <- as.integer(sub("partition_", "", py$partition))
  py <- merge(py, unique(r[c("dataset", "partition")]))
  add("P1 Table S2", NULL, cellmeans(r, c("dataset", "model"), "COR_within"), cellmeans(py, c("dataset", "model"), "COR_within"))
}
# ---- Part I Section 3.2 ---------------------------------------------------------------
r <- rout("P1_S32_psd_diagnostics.csv")
if (!is.null(r)) {
  py <- merge(read.csv(file.path(PY1, "psd_diagnostics.csv")), unique(r[c("dataset", "kernel")]))
  for (v in c("negative_eigenvalues", "spectral_negativity", "relative_psd_projection_change"))
    add(paste("P1 Section 3.2", v), NULL, cellmeans(r, c("dataset", "kernel"), v), cellmeans(py, c("dataset", "kernel"), v))
}
# ---- Part II Table 3 -------------------------------------------------------------------
r <- rout("P2_T3_cross_covariance.csv")
if (!is.null(r)) {
  py <- pyout(PY2, "^part2_cross_v2_.*\\.csv$"); py <- merge(py, unique(r[c("dataset", "env")]))
  em <- function(x) { a <- aggregate(COR ~ dataset + env + model, x, mean); cellmeans(a, c("dataset", "model"), "COR") }
  add("P2 Table 3", NULL, em(r), em(py))
}
# ---- Part II Table 4 -------------------------------------------------------------------
r <- rout("P2_T4_identifier_neighbours.csv")
if (!is.null(r)) {
  py <- rbind(read.csv(file.path(PY2, "part2_neighbour_cuevas.csv")), read.csv(file.path(PY2, "part2_neighbour_ge.csv")))
  py$ratio <- py$nb_cross / py$same_line_cross_env
  for (v in c("nb_within", "ratio")) add(paste("P2 Table 4", v), NULL, cellmeans(r, c("dataset", "lag"), v), cellmeans(py, c("dataset", "lag"), v))
  r2 <- rout("P2_T4_identifier_relatedness.csv"); p2 <- read.csv(file.path(PY2, "part2_id_relatedness.csv"))
  add("P2 Table 4 Related", NULL, setNames(r2$adjacent_in_top5pct, r2$dataset), setNames(p2$adjacent_in_top5pct, p2$dataset))
}
# ---- Part II Tables 5-6 ------------------------------------------------------------------
r <- rout("P2_T5_T6_identifier_split.csv")
if (!is.null(r)) {
  py <- pyout(PY2, "^part2_id_split_.*\\.csv$"); py <- merge(py, unique(r[c("dataset", "scenario")]))
  for (v in c("COR_within", "PM10_within")) add(paste("P2 Tables 5-6", v), NULL, cellmeans(r, c("dataset", "scenario", "model"), v), cellmeans(py, c("dataset", "scenario", "model"), v))
}
# ---- Part II Table 7 / S1 ------------------------------------------------------------------
r <- rout("P2_T7_S1_breeding_scenarios.csv")
if (!is.null(r)) {
  py <- pyout(PY2, "^part2_scenarios_.*\\.csv$")
  py$ec <- ifelse(grepl("ec_previous", py$.file), "previous", ifelse(grepl("ec_train_mean", py$.file), "train_mean", "observed"))
  py <- merge(py, unique(r[c("dataset", "ec")]))
  py$grp <- ifelse(startsWith(py$split, "new_lines"), "new_lines", py$split); r$grp <- ifelse(startsWith(r$split, "new_lines"), "new_lines", r$split)
  for (v in c("COR_within", "PM10_within")) add(paste("P2 Table 7/S1", v), NULL, cellmeans(r, c("dataset", "ec", "grp", "subset", "model"), v), cellmeans(py, c("dataset", "ec", "grp", "subset", "model"), v))
}

res <- do.call(rbind, cmp)
if (is.null(res)) stop("No R outputs found in ", OUT_DIR, ". Run the R scripts first.")
write_out(res, "comparison_R_vs_Python.csv")
summ <- do.call(rbind, lapply(split(res, res$table), function(x) data.frame(table = x$table[1], cells = nrow(x),
  max_abs_difference = max(abs(x$difference), na.rm = TRUE), mean_abs_difference = mean(abs(x$difference), na.rm = TRUE))))
rownames(summ) <- NULL
cat("\nLargest absolute difference between R and Python, per table:\n")
print(summ, digits = 3)
