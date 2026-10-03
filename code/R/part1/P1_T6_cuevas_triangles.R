# =============================================================================
# Part I - Table 6 (Section 3.7) and Part II Section 3.1: the genomic-pedigree product of
# Cuevas et al. (2025) on their own data (Gaussian kernel K, pedigree A, lines in file order).
# Python equivalent: cuevas_check_files.py <file.Rdata>
#
# Within each environment, five random 80:20 splits of the lines (splits/cuevas_splits.csv.xz)
# are predicted with single-environment models:
#   KA            K + A                                    baseline
#   CP_file       K + A + C + P, triangles of K A in the order of the file
#   CP_shuffled   K + A + C + P, triangles built after shuffling the line order
#                 (10 shuffles, splits/cuevas_shuffles.csv.xz)
#   CP_only       C + P (file order)
#   ID            K + A + identifier-proximity kernel
# C = D + U + U' and P = D + L + L' of M = K A (Cuevas et al. 2025, Equations 5-6).
# Output: P1_T6_cuevas_triangles.csv. Table 6 = mean COR over splits, then over environments;
#   "File order > all shuffles" = splits in which CP_file beats every CP_shuffled.
# Usage:  Rscript code/R/part1/P1_T6_cuevas_triangles.R [--datasets Wheat_599 EYT_16_17 ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", c("Wheat_599", "EYT_16_17", "EYT_22_23", "EYT_23_24"))
cspl <- read_split("cuevas_splits.csv.xz")
cshf <- read_split("cuevas_shuffles.csv.xz")

rows <- list()
for (nm in datasets) {
  if (!file.exists(file.path(DATA_DIR, "Cuevas2025", paste0(nm, ".rds")))) { log_msg(nm, "not in data/Cuevas2025, skipped"); next }
  x <- load_dataset(nm)
  K <- unname(x$K); A <- unname(x$A); lines <- x$pheno$Line
  tri <- triangles(K %*% A)
  kernels <- list(K = K, A = A, C = tri$C, P = tri$P)
  models <- list(KA = c("K", "A"), CP_file = c("K", "A", "C", "P"), CP_only = c("C", "P"))
  for (k in 1:10) {                                   # triangles built in a shuffled line order
    o <- cshf$line[cshf$dataset == nm & cshf$shuffle == k]
    t_o <- triangles(K[o, o] %*% A[o, o]); pos <- order(o)
    kernels[[paste0("Cs", k)]] <- t_o$C[pos, pos]; kernels[[paste0("Ps", k)]] <- t_o$P[pos, pos]
    models[[sprintf("CP_shuffled_%02d", k)]] <- c("K", "A", paste0("Cs", k), paste0("Ps", k))
  }
  kernels$ID <- identifier_kernel(lines, id_digits(lines), h_min = 0)
  models$ID <- c("K", "A", "ID")
  for (env in x$env_names) {
    yv <- x$pheno[[env]]; ok <- which(!is.na(yv))
    if (length(ok) < 50) { log_msg(nm, env, "no phenotypes, skipped"); next }
    y <- yv[ok]
    Ks <- lapply(kernels, function(M) M[ok, ok])
    for (fi in 1:5) {
      te <- cspl$index[cspl$dataset == nm & cspl$env == env & cspl$split == fi]
      tr <- setdiff(seq_along(y), te)
      blocks <- lapply(Ks, psd_blocks, tr = tr, te = te)
      X1 <- matrix(1, length(y), 1)
      for (mod in names(models)) {
        used <- models[[mod]][!vapply(blocks[models[[mod]]], is.null, TRUE)]
        fit <- reml_fixed(y[tr], X1[tr, , drop = FALSE], blocks[used], X1[te, , drop = FALSE])
        rows[[length(rows) + 1]] <- data.frame(dataset = nm, env = env, split = fi, model = mod,
                                               COR = cor(fit$pred, y[te]))
      }
      log_msg(nm, env, "split", fi, "done")
    }
  }
  res <- do.call(rbind, rows)
  res$group <- ifelse(startsWith(res$model, "CP_shuffled"), "CP_shuffled", res$model)
  write_out(res, "P1_T6_cuevas_triangles.csv")
}
res <- do.call(rbind, rows); res$group <- ifelse(startsWith(res$model, "CP_shuffled"), "CP_shuffled", res$model)
env_means <- aggregate(COR ~ dataset + env + group, res, mean)
print(round(tapply(env_means$COR, list(env_means$dataset, env_means$group), mean), 3))
best_shuf <- aggregate(COR ~ dataset + env + split, res[res$group == "CP_shuffled", ], max)
cpf <- merge(res[res$model == "CP_file", ], best_shuf, by = c("dataset", "env", "split"), suffixes = c("", "_maxshuf"))
cat("\nSplits in which the file order beats every shuffle:\n")
print(tapply(cpf$COR > cpf$COR_maxshuf, cpf$dataset, function(v) sprintf("%d of %d", sum(v), length(v))))
