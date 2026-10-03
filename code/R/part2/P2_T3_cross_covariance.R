# =============================================================================
# Part II - Table 3 (Section 3.1): the complete genomic-pedigree product K A entered as a
# valid, order-independent cross-covariance (Equation 1), the spectral kernels, and the
# identifier kernel, on the data of Cuevas et al. (2025).
# Python equivalent: part2/code/part2_cross_ka_v2.py <file.Rdata>
#
# Models (single environment, five 80:20 splits per environment, splits/cuevas_splits.csv.xz):
#   KA        K + A
#   CP_file   K + A + C + P (triangles of K A in file order; reported in Part I, Table 6)
#   CROSS     K + A + gamma (K A + (K A)'),   gamma of free sign
#   SPECTRAL  K + A + K^1/2 A K^1/2 + A^1/2 K A^1/2
#   ID        K + A + identifier kernel
#   ID_CROSS  K + A + identifier kernel + gamma (K A + (K A)')
# The cross term is scaled by its spectral norm on the training block and every model is fitted
# with reml_robust (several starting values, best REML likelihood). gamma is tested by a
# likelihood-ratio test with 1 df on the training data: 2 (logL_CROSS - logL_KA) > 3.841.
# Output: P2_T3_cross_covariance.csv
# Usage:  Rscript code/R/part2/P2_T3_cross_covariance.R [--datasets Wheat_599 ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", c("Wheat_599", "EYT_16_17", "EYT_22_23", "EYT_23_24"))
cspl <- read_split("cuevas_splits.csv.xz")

rows <- list()
for (nm in datasets) {
  if (!file.exists(file.path(DATA_DIR, "Cuevas2025", paste0(nm, ".rds")))) { log_msg(nm, "not available, skipped"); next }
  x <- load_dataset(nm)
  K <- unname(x$K); A <- unname(x$A); lines <- x$pheno$Line
  KA <- K %*% A
  tri <- triangles(KA)
  Kh <- sqrtm_psd(K); Ah <- sqrtm_psd(A)
  kernels <- list(K = K, A = A, C = tri$C, P = tri$P, S1 = Kh %*% A %*% Kh, S2 = Ah %*% K %*% Ah,
                  ID = identifier_kernel(lines, id_digits(lines), h_min = 0))
  CROSS <- KA + t(KA)
  models <- list(KA = list(c("K", "A"), FALSE), CP_file = list(c("K", "A", "C", "P"), FALSE),
                 CROSS = list(c("K", "A"), TRUE), SPECTRAL = list(c("K", "A", "S1", "S2"), FALSE),
                 ID = list(c("K", "A", "ID"), FALSE), ID_CROSS = list(c("K", "A", "ID"), TRUE))
  for (env in x$env_names) {
    yv <- x$pheno[[env]]; ok <- which(!is.na(yv))
    if (length(ok) < 50) next
    y <- yv[ok]
    Ks <- lapply(kernels, function(M) M[ok, ok]); Xs <- CROSS[ok, ok]
    for (fi in 1:5) {
      te <- cspl$index[cspl$dataset == nm & cspl$env == env & cspl$split == fi]
      tr <- setdiff(seq_along(y), te)
      blocks <- lapply(Ks, psd_blocks, tr = tr, te = te)
      Stt <- Xs[tr, tr]; sc <- spectral_scale(Stt)
      sx <- list(tt = Stt / sc, ut = Xs[te, tr, drop = FALSE] / sc)
      for (mod in names(models)) {
        pn <- models[[mod]][[1]]; used <- pn[!vapply(blocks[pn], is.null, TRUE)]
        fit <- reml_robust(y[tr], blocks[used], if (models[[mod]][[2]]) list(sx) else list())
        rows[[length(rows) + 1]] <- data.frame(dataset = nm, env = env, split = fi, model = mod,
          COR = cor(fit$pred, y[te]), gamma = if (length(fit$gamma)) fit$gamma[1] else NA,
          logL = -fit$nll, converged = fit$converged, starts_ok = fit$starts_ok)
      }
      log_msg(nm, env, "split", fi, "done")
    }
  }
  write_out(do.call(rbind, rows), "P2_T3_cross_covariance.csv")
}
res <- do.call(rbind, rows)
env_means <- aggregate(COR ~ dataset + env + model, res, mean)
print(round(tapply(env_means$COR, list(env_means$dataset, env_means$model), mean), 3))
lrt <- function(a, b) {                              # training sets with significant gamma (P < 0.05)
  m <- merge(res[res$model == a, ], res[res$model == b, ], by = c("dataset", "env", "split"))
  tapply(2 * (m$logL.x - m$logL.y) > qchisq(0.95, 1), m$dataset, function(v) sprintf("%d of %d", sum(v), length(v)))
}
cat("\ngamma significant (CROSS vs KA):\n"); print(lrt("CROSS", "KA"))
cat("gamma significant (ID_CROSS vs ID):\n"); print(lrt("ID_CROSS", "ID"))
cat("converged:", sprintf("%.1f%%", 100 * mean(res$converged[res$model %in% c("CROSS", "ID_CROSS")])), "\n")
