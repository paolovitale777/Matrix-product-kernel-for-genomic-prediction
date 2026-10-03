# =============================================================================
# Part II - Table 5 (untested lines in tested environments, 50:50 partitions), Table 6
# (untested environments, leave one environment out) and the "Share" column; PM10 values
# for Supplementary Table S2.
# Python equivalent: part2/code/part2_id_split.py --datasets ...
#
# Models (Table 2 of Part II; all include K_E):
#   M1 = B_Hadamard    linear genomic kernel G and G o E
#   M2 = GK_Hadamard   Gaussian genomic kernel GK and GK o E
#   M3 = ID_common     M1 + identifier kernel shared by all environments (K_id)
#   M4 = ID_env        M1 + environment-specific identifier kernel (K_id o S, S = same environment)
#   M5 = GK_ID_env     M2 + environment-specific identifier kernel
#   R6_both            M1 + K_id + K_id o S (used only for the variance share)
#   GK_ID_common       M2 + K_id (not reported)
# Share = var(ID_E) / (var(ID) + var(ID_E)) in R6_both, mean over the five partitions.
# Output: P2_T5_T6_identifier_split.csv (scenario cv1 = Table 5, loeo = Table 6)
# Usage:  Rscript code/R/part2/P2_T5_T6_identifier_split.R [--datasets ...] [--scenarios cv1 loeo]
# =============================================================================
source("code/R/00_setup.R")
datasets  <- arg_values("--datasets", GE_DATASETS)
scenarios <- arg_values("--scenarios", c("cv1", "loeo"))
parts <- read_split("partitions_5050.csv.xz")
MODELS <- list(B_Hadamard = c("E", "G", "H"), R6_both = c("E", "G", "H", "ID", "ID_E"),
               ID_common = c("E", "G", "H", "ID"), ID_env = c("E", "G", "H", "ID_E"),
               GK_Hadamard = c("E", "GK", "HK"), GK_ID_common = c("E", "GK", "HK", "ID"),
               GK_ID_env = c("E", "GK", "HK", "ID_E"))

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  bc <- build_components(d, include_original = FALSE)
  y <- ph$GY; envs <- ph$Env; gids <- ph$GID
  K_id <- identifier_kernel(gids, id_number(gids))
  S <- outer(envs, envs, "==") * 1
  GK <- gaussian_kernel(bc$comps$G)
  kernels <- list(E = bc$comps$E, G = bc$comps$G, H = bc$comps$H, ID = K_id, ID_E = K_id * S,
                  GK = GK, HK = bc$comps$E * GK)
  X <- matrix(1, length(y), 1)
  for (scen in scenarios) {
    sp <- if (scen == "cv1") lapply(1:5, function(fi) { m <- test_mask_5050(parts, nm, fi, gids); list(paste0("partition_", fi), which(!m), which(m)) })
          else lapply(sortc(envs), function(e) list(paste0("leave_", e), which(envs != e), which(envs == e)))
    for (s in sp) {
      tr <- s[[2]]; te <- s[[3]]
      blocks <- Filter(Negate(is.null), lapply(kernels, psd_blocks, tr = tr, te = te))
      for (mod in names(MODELS)) {
        used <- intersect(MODELS[[mod]], names(blocks))
        fit <- reml_fixed(y[tr], X[tr, , drop = FALSE], blocks[used], X[te, , drop = FALSE])
        ev <- evaluate(y[te], fit$pred, envs[te])
        v <- setNames(fit$theta, c(used, "resid"))
        rows[[length(rows) + 1]] <- data.frame(dataset = nm, scenario = scen, split = s[[1]], model = mod,
          var_ID = if ("ID" %in% used) v[["ID"]] else NA, var_ID_E = if ("ID_E" %in% used) v[["ID_E"]] else NA,
          converged = fit$converged, t(ev))
      }
      log_msg(nm, scen, s[[1]], "done")
    }
    write_out(do.call(rbind, rows), "P2_T5_T6_identifier_split.csv")
  }
}
res <- do.call(rbind, rows)
for (scen in scenarios) {
  cat("\n", if (scen == "cv1") "Table 5 (untested lines in tested environments)" else "Table 6 (untested environments)", "\n")
  r <- res[res$scenario == scen & res$model %in% c("B_Hadamard", "GK_Hadamard", "ID_common", "ID_env", "GK_ID_env"), ]
  print(round(tapply(r$COR_within, list(r$dataset, r$model), mean), 3))
}
b <- res[res$scenario == "cv1" & res$model == "R6_both", ]
if (nrow(b)) { cat("\nShare (environment-specific identifier variance):\n"); print(round(tapply(b$var_ID_E / (b$var_ID + b$var_ID_E), b$dataset, mean), 2)) }
