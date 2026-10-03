# =============================================================================
# Part I - Supplementary Table S2: order-independent models built from the connectivity
# vector c of Proposition 1 (c_s = total genomic relationship of the line of record s with
# all lines in the data).
# Python equivalent: recover_gain.py --datasets EYT_1 Wheat_1 Japonica
# Models (all with E, G, H):
#   R1_conn_fixed   c as a fixed covariate within each environment
#   R2_env_G        one genomic kernel per environment
#   R3_conn_kernel  c c' and (c c') o K_E
#   R4_conn_weight  D K_G D and (D K_G D) o K_E,  D = diag(c / mean c)
# plus B_Hadamard and F_original as references.
# Output: P1_S2_connectivity_models.csv
# Usage:  Rscript code/R/part1/P1_S2_connectivity_models.R [--datasets EYT_1 Wheat_1 Japonica]
# =============================================================================
source("code/R/00_setup.R")
datasets   <- arg_values("--datasets", c("EYT_1", "Wheat_1", "Japonica"))
partitions <- as.integer(arg_values("--partitions", 1:5))
parts <- read_split("partitions_5050.csv.xz")

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  bc <- build_components(d, include_original = TRUE)
  if (bc$info$degenerate) next
  y <- ph$GY; envs <- ph$Env; gids <- ph$GID; env_levels <- sortc(envs)
  G <- genomic_matrix(d); lines <- sortc(gids)
  tot <- rowSums(G[lines, lines]); c_rec <- tot[match(gids, lines)]          # connectivity()
  cz <- (c_rec - mean(c_rec)) / sqrt(mean((c_rec - mean(c_rec))^2))
  KG <- bc$comps$G; KE <- bc$comps$E
  extra <- list()
  for (e in env_levels) { s <- as.numeric(envs == e); extra[[paste0("G_", e)]] <- KG * outer(s, s) }
  extra$CC <- outer(cz, cz); extra$CC_E <- extra$CC * KE
  w <- c_rec / mean(c_rec); extra$DGD <- KG * outer(w, w); extra$DGD_E <- extra$DGD * KE
  kernels <- c(bc$comps[c("E", "G", "H", "C0", "P0")], extra)
  X_int <- matrix(1, length(y), 1)
  X_conn <- cbind(X_int, cz * sapply(env_levels, function(e) as.numeric(envs == e)))
  models <- list(B_Hadamard = list(c("E", "G", "H"), X_int), F_original = list(c("E", "G", "H", "C0", "P0"), X_int),
                 R1_conn_fixed = list(c("E", "G", "H"), X_conn),
                 R2_env_G = list(c("E", "G", "H", paste0("G_", env_levels)), X_int),
                 R3_conn_kernel = list(c("E", "G", "H", "CC", "CC_E"), X_int),
                 R4_conn_weight = list(c("E", "G", "H", "DGD", "DGD_E"), X_int))
  for (fi in partitions) {
    is_te <- test_mask_5050(parts, nm, fi, gids); tr <- which(!is_te); te <- which(is_te)
    blocks <- Filter(Negate(is.null), lapply(kernels, psd_blocks, tr = tr, te = te))
    for (mod in names(models)) {
      used <- intersect(models[[mod]][[1]], names(blocks)); X <- models[[mod]][[2]]
      fit <- reml_fixed(y[tr], X[tr, , drop = FALSE], blocks[used], X[te, , drop = FALSE])
      ev <- evaluate(y[te], fit$pred, envs[te])
      rows[[length(rows) + 1]] <- data.frame(dataset = nm, partition = fi, model = mod, converged = fit$converged, t(ev))
    }
    log_msg(nm, "partition", fi, "done")
  }
  write_out(do.call(rbind, rows), "P1_S2_connectivity_models.csv")
}
res <- do.call(rbind, rows)
print(round(tapply(res$COR_within, list(res$dataset, res$model), mean), 3))
