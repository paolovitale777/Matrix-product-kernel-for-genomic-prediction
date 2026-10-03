# =============================================================================
# Part I - Section 3.2: covariance validity of the original and proposed components.
# Python equivalent: analysis_pipeline.py psd
# For each kernel: number of negative / positive eigenvalues, spectral negativity
# (sum |negative eigenvalues| / sum |eigenvalues|) and relative change caused by projection
# onto the PSD cone, sqrt(sum lambda_-^2 / sum lambda^2) (Supplementary S1, Equations S8-S9).
# The proposed components have rank <= 2 x (number of environments); Python obtains their
# eigenvalues with a randomised range finder, here they are computed the same way with R's
# random generator (the eigenvalues are the same up to rounding).
# Output: P1_S32_psd_diagnostics.csv
# Usage:  Rscript code/R/part1/P1_S32_psd_diagnostics.R [--datasets ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", GE_DATASETS)

spectrum <- function(ev, n, label, nm) {
  tol <- max(1, max(abs(ev), 0)) * n * EPS
  neg <- ev < -tol; pos <- ev > tol
  data.frame(dataset = nm, kernel = label, dimension = n, minimum_eigenvalue = min(ev, 0),
             maximum_eigenvalue = max(ev, 0), negative_eigenvalues = sum(neg), positive_eigenvalues = sum(pos),
             zero_eigenvalues = n - sum(neg) - sum(pos),
             spectral_negativity = if (sum(abs(ev))) sum(abs(ev[neg])) / sum(abs(ev)) else 0,
             relative_psd_projection_change = if (sum(ev^2)) sqrt(sum(ev[neg]^2) / sum(ev^2)) else 0)
}
low_rank_eigenvalues <- function(A, rank_bound, seed) {
  n <- nrow(A); if (!any(A != 0)) return(numeric(0))
  set.seed(seed); l <- min(n, max(2, rank_bound + 4))
  qrY <- qr(A %*% matrix(rnorm(n * l), n, l))
  Rm <- qr.R(qrY); keep <- abs(diag(Rm)) > max(1, norm(Rm, "2")) * n * EPS
  Q <- qr.Q(qrY)[, keep, drop = FALSE]
  if (!ncol(Q)) return(numeric(0))
  S <- crossprod(Q, A %*% Q)
  eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
}

rows <- list()
for (i in seq_along(datasets)) {
  nm <- datasets[i]
  d <- load_ge(nm)
  bc <- build_components(d, include_original = TRUE)
  n <- nrow(bc$comps$G)
  if (bc$info$degenerate) { log_msg(nm, "degenerate, skipped"); next }
  for (k in c("C0", "P0")) {
    ev <- eigen(bc$comps[[k]], symmetric = TRUE, only.values = TRUE)$values
    rows[[length(rows) + 1]] <- spectrum(ev, n, paste0("original_", k), nm)
  }
  for (k in intersect(c("CSP", "PSP", "TSP", "SGE"), names(bc$comps))) {
    ev <- low_rank_eigenvalues(bc$comps[[k]], 2 * length(bc$info$env_levels), 8800 + 10 * i)
    r <- spectrum(ev, n, paste0("proposed_", k), nm)
    r$mean_abs_diagonal <- mean(abs(diag(bc$comps[[k]])))
    rows[[length(rows) + 1]] <- r
  }
  log_msg(nm, "done")
  all_cols <- unique(unlist(lapply(rows, names)))
  write_out(do.call(rbind, lapply(rows, function(x) { x[setdiff(all_cols, names(x))] <- NA; x[all_cols] })), "P1_S32_psd_diagnostics.csv")
}
