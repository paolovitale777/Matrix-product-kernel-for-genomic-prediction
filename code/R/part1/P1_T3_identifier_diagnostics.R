# =============================================================================
# Part I - Table 3, columns "Hadamard (B)", "F, lines shuffled", "Identifier kernel (R6)",
# and Section 3.4 (model R5).
# Python equivalent: id_order_test.py
#
# Models (all with the environmental and genomic main effects):
#   B_Hadamard    E, G, H                                           reference
#   F_original    E, G, H, C0, P0 built in the order of the file
#   F_shuffled    same, with the lines shuffled within each environment (environment blocks
#                 kept; splits/shuffle_within_env.csv.xz)
#   R5_id_fixed   E, G, H + identifier number as a fixed covariate within each environment
#   R6_id_kernel  E, G, H + K_id and K_id o K_E,  K_id(j,k) = exp(-|rank_j - rank_k| / h),
#                 h = 5% of the lines (diagnostic, order independent)
# Output: P1_T3_identifier_diagnostics.csv
# Usage:  Rscript code/R/part1/P1_T3_identifier_diagnostics.R [--datasets ...] [--partitions ...]
# =============================================================================
source("code/R/00_setup.R")
datasets   <- arg_values("--datasets", setdiff(GE_DATASETS, "Indica"))
partitions <- as.integer(arg_values("--partitions", 1:5))
parts <- read_split("partitions_5050.csv.xz")
shuf  <- read_split("shuffle_within_env.csv.xz")

triangles_in_order <- function(d, order) {        # C0, P0 built in `order`, mapped back to file order
  ph <- d$pheno[order, ]; rownames(ph) <- NULL
  bc <- build_components(d, ph, include_original = TRUE)$comps
  pos <- order(order)
  list(C0 = bc$C0[pos, pos], P0 = bc$P0[pos, pos])
}

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  bc <- build_components(d, include_original = TRUE)
  if (bc$info$degenerate) { log_msg(nm, "degenerate, skipped"); next }
  y <- ph$GY; envs <- ph$Env; gids <- ph$GID
  env_levels <- sortc(envs)
  num <- id_number(gids)
  st <- triangles_in_order(d, shuf$record[shuf$dataset == nm])
  K_id <- identifier_kernel(gids, num)
  kernels <- c(bc$comps[c("E", "G", "H", "C0", "P0")],
               list(C0_shuf = st$C0, P0_shuf = st$P0, ID = K_id, ID_E = K_id * bc$comps$E))
  X_int <- matrix(1, length(y), 1)
  z <- (num - mean(num)) / sqrt(mean((num - mean(num))^2))
  X_id <- cbind(X_int, z * sapply(env_levels, function(e) as.numeric(envs == e)))
  models <- list(B_Hadamard = list(c("E", "G", "H"), X_int),
                 F_original = list(c("E", "G", "H", "C0", "P0"), X_int),
                 F_shuffled = list(c("E", "G", "H", "C0_shuf", "P0_shuf"), X_int),
                 R5_id_fixed = list(c("E", "G", "H"), X_id),
                 R6_id_kernel = list(c("E", "G", "H", "ID", "ID_E"), X_int))
  for (fi in partitions) {
    is_te <- test_mask_5050(parts, nm, fi, gids)
    tr <- which(!is_te); te <- which(is_te)
    blocks <- Filter(Negate(is.null), lapply(kernels, psd_blocks, tr = tr, te = te))
    for (mod in names(models)) {
      used <- intersect(models[[mod]][[1]], names(blocks)); X <- models[[mod]][[2]]
      fit <- reml_fixed(y[tr], X[tr, , drop = FALSE], blocks[used], X[te, , drop = FALSE])
      ev <- evaluate(y[te], fit$pred, envs[te])
      rows[[length(rows) + 1]] <- data.frame(dataset = nm, partition = fi, model = mod,
        kernels_used = paste(used, collapse = "|"), converged = fit$converged, t(ev))
      log_msg(sprintf("%s p%d %-12s COR_within = %.4f", nm, fi, mod, ev[["COR_within"]]))
    }
  }
  write_out(do.call(rbind, rows), "P1_T3_identifier_diagnostics.csv")
}
res <- do.call(rbind, rows)
write_out(paired_differences(res), "P1_T3_identifier_paired_differences.csv")
print(round(tapply(res$COR_within, list(res$dataset, res$model), mean), 3))
