# =============================================================================
# Part I - Table 4 and Supplementary Table S1: 50:50 benchmark (untested lines in tested
# environments).
# Python equivalent: code/python/part1/code/benchmark_v4.py --scenario ulte
#   (the models of Table 4: A_main, B_Hadamard, F_original, M_masks = "Proposed").
#
# For each dataset and each of the five cultivar partitions (splits/partitions_5050.csv.xz):
#   1. build the record-level kernels in file order (build_components);
#   2. project every kernel onto the PSD cone using the training block only (psd_blocks);
#   3. fit each model by REML on the training records and predict the test records;
#   4. score within each environment (COR, NRMSE, PM10) and pooled.
# Output: code/R/outputs/P1_T4_S1_benchmark_5050.csv (one row per dataset x partition x model).
#   Table 4  = mean of COR_within (and PM10_within) over the five partitions.
#   Table S1 = mean of COR_pooled over the five partitions.
# Usage:  Rscript code/R/part1/P1_T4_S1_benchmark_5050.R [--datasets ...] [--partitions 1 2 ...]
# =============================================================================
source("code/R/00_setup.R")
datasets   <- arg_values("--datasets", GE_DATASETS)
partitions <- as.integer(arg_values("--partitions", 1:5))
parts <- read_split("partitions_5050.csv.xz")

model_list <- function(info) {                     # Python benchmark_v4.model_list (Table 4 models)
  m <- list(A_main = c("E", "G"), B_Hadamard = c("E", "G", "H"))
  if (info$degenerate) return(m)
  m$F_original <- c("E", "G", "H", "C0", "P0")
  m$M_masks <- if (info$directed) c("E", "G", "H", "CSP", "PSP", "TSP") else c("E", "G", "H", "SGE")
  m
}

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  bc <- build_components(d); comps <- bc$comps; info <- bc$info
  models <- model_list(info)
  log_msg(nm, "degenerate =", info$degenerate, " directed =", info$directed)
  y <- ph$GY
  for (fi in partitions) {
    is_te <- test_mask_5050(parts, nm, fi, ph$GID)
    tr <- which(!is_te); te <- which(is_te)
    needed <- unique(unlist(models))
    blocks <- Filter(Negate(is.null), setNames(lapply(needed, function(k) psd_blocks(comps[[k]], tr, te)), needed))
    for (mod in names(models)) {
      used <- intersect(models[[mod]], names(blocks))
      fit <- reml_mixed(y[tr], blocks[used])
      ev <- evaluate(y[te], fit$pred, ph$Env[te])
      rows[[length(rows) + 1]] <- data.frame(dataset = nm, partition = fi, model = mod,
        kernels_used = paste(used, collapse = "|"), converged = fit$converged, t(ev))
      log_msg(sprintf("  %s partition %d %-11s COR_within = %.4f", nm, fi, mod, ev[["COR_within"]]))
    }
  }
  write_out(do.call(rbind, rows), "P1_T4_S1_benchmark_5050.csv")
}
res <- do.call(rbind, rows)
write_out(paired_differences(res), "P1_T4_S1_benchmark_5050_paired_differences.csv")
cat("\nTable 4 (COR_within, mean of partitions):\n")
print(round(tapply(res$COR_within, list(res$dataset, res$model), mean), 3))
cat("\nTable S1 (COR_pooled, mean of partitions):\n")
print(round(tapply(res$COR_pooled, list(res$dataset, res$model), mean), 3))
