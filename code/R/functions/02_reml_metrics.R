# =============================================================================
# 02_reml_metrics.R - PSD projection, REML fitting and accuracy metrics
#
# R translation of code/python/part1/code/benchmark_ulte.py (psd_blocks, reml_fit_predict,
# basic_metrics, evaluate, paired_differences), recover_gain.reml_fixed,
# benchmark_v4.reml_mixed and part2/code/part2_reml.reml_robust.
#
# Model:  y = X b + sum_k u_k + e,  u_k ~ N(0, theta_k K_k),  e ~ N(0, theta_e I),
#         optionally + free-sign cross terms gamma_s S_s in V (Part II, Equation 1).
# REML criterion (minimised):
#   0.5 * [ log|V| + log|X' V^-1 X| + r' V^-1 r + (n - p) log(2 pi) ],  r = y - X b_GLS.
# Variances are optimised on the log scale with L-BFGS-B and the same bounds and starting
# values as the Python code. Python uses finite-difference gradients; here the exact REML
# gradient is used (faster in R, same optimum):
#   d/d psi_j = 0.5 * [ tr(P M_j) - y' P M_j P y ],   P = V^-1 - V^-1 X (X'V^-1X)^-1 X' V^-1,
# with M_j = theta_j K_j for log-variances and M_j = S_j for cross terms.
# =============================================================================

# ---- training-anchored PSD projection (Python benchmark_ulte.psd_blocks) -----
# Eigen-decompose the training block, keep positive eigenvalues, extend to the test-by-train
# block through the same eigenvectors, and scale by the training mean diagonal. Test
# phenotypes are never used.
psd_blocks <- function(A, tr, te) {
  Att <- (A[tr, tr] + t(A[tr, tr])) / 2
  Aut <- A[te, tr, drop = FALSE]
  e <- eigen(Att, symmetric = TRUE)
  w <- e$values
  tol <- max(1, max(abs(w))) * length(tr) * EPS
  keep <- w > tol
  if (!any(keep)) return(NULL)
  V <- e$vectors[, keep, drop = FALSE]; w <- w[keep]
  Ktt <- V %*% (w * t(V))
  Kut <- (Aut %*% V) %*% t(V)
  scale <- mean(diag(Ktt))
  if (scale <= tol) return(NULL)
  list(tt = Ktt / scale, ut = Kut / scale, rank = length(w))
}

# ---- REML core -----------------------------------------------------------------
# blocks: list of list(tt, ut);  signed: list of list(tt, ut);  X, Xte: fixed-effect designs.
reml_core <- function(y, X, blocks, signed, par0, lower, upper, maxit, ftol, const = TRUE) {
  n <- length(y); m <- length(blocks); s <- length(signed); p <- ncol(X)
  Klist <- c(lapply(blocks, `[[`, "tt"), list(diag(n)))
  Slist <- lapply(signed, `[[`, "tt")
  cache <- new.env()
  evalV <- function(par) {
    key <- paste(format(par, digits = 17), collapse = ",")
    if (identical(cache$key, key)) return(cache$res)
    th <- exp(par[seq_len(m + 1)]); gam <- if (s) par[m + 1 + seq_len(s)] else numeric(0)
    V <- th[m + 1] * diag(n)
    for (k in seq_len(m)) V <- V + th[k] * Klist[[k]]
    for (j in seq_len(s)) V <- V + gam[j] * Slist[[j]]
    R <- tryCatch(chol(V), error = function(e) NULL)
    res <- NULL
    if (!is.null(R)) {
      Vi <- chol2inv(R)
      ViX <- Vi %*% X
      XtViX <- crossprod(X, ViX)
      ld <- determinant(XtViX, logarithm = TRUE)
      if (ld$sign > 0) {
        b <- solve(XtViX, crossprod(ViX, y))
        r <- y - X %*% b
        Pr <- Vi %*% r                                      # P y = V^-1 (y - X b)
        f <- 0.5 * (2 * sum(log(diag(R))) + as.numeric(ld$modulus) + sum(r * Pr) +
                    if (const) (n - p) * log(2 * pi) else 0)
        res <- list(f = f, th = th, gam = gam, Vi = Vi, ViX = ViX, XtViXi = solve(XtViX), Pr = Pr, b = b)
      }
    }
    cache$key <- key; cache$res <- res
    res
  }
  fn <- function(par) { r <- evalV(par); if (is.null(r)) 1e100 else r$f }
  gr <- function(par) {
    r <- evalV(par)
    if (is.null(r)) return(rep(0, length(par)))
    trP <- function(M) sum(r$Vi * M) - sum(diag(r$XtViXi %*% crossprod(r$ViX, M %*% r$ViX)))
    g <- numeric(length(par))
    for (k in seq_len(m + 1)) {
      M <- Klist[[k]]
      g[k] <- 0.5 * r$th[k] * (trP(M) - sum(r$Pr * (M %*% r$Pr)))
    }
    for (j in seq_len(s)) g[m + 1 + j] <- 0.5 * (trP(Slist[[j]]) - sum(r$Pr * (Slist[[j]] %*% r$Pr)))
    g
  }
  fit <- optim(par0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = maxit, factr = ftol / EPS, pgtol = 1e-5, lmm = 10))
  list(fit = fit, eval = evalV(fit$par), m = m, s = s)
}

reml_predict <- function(core, X, Xte, blocks, signed, y) {
  r <- core$eval
  alpha <- r$Vi %*% (y - X %*% r$b)
  pred <- as.numeric(Xte %*% r$b)
  for (k in seq_along(blocks)) pred <- pred + r$th[k] * as.numeric(blocks[[k]]$ut %*% alpha)
  for (j in seq_along(signed)) pred <- pred + r$gam[j] * as.numeric(signed[[j]]$ut %*% alpha)
  pred
}

# Python benchmark_ulte.reml_fit_predict (intercept only; maxiter 120, ftol 1e-7)
reml_fit_predict <- function(y, blocks, maxit = 120, ftol = 1e-7) {
  n <- length(y); X <- matrix(1, n, 1); Xte <- matrix(1, nrow(blocks[[1]]$ut), 1)
  vy <- max(mean((y - mean(y))^2), 1e-8); m <- length(blocks)
  lo <- log(vy * 1e-7); hi <- log(vy * 1e3)
  core <- reml_core(y, X, blocks, list(), rep(log(vy / (m + 1)), m + 1), rep(lo, m + 1), rep(hi, m + 1), maxit, ftol)
  list(pred = reml_predict(core, X, Xte, blocks, list(), y), theta = core$eval$th,
       converged = core$fit$convergence == 0, objective = core$fit$value)
}

# Python recover_gain.reml_fixed (general fixed effects; maxiter 200, ftol 1e-8; no constant)
reml_fixed <- function(y, X, blocks, Xte, maxit = 200, ftol = 1e-8) {
  vy <- max(mean((y - mean(y))^2), 1e-8); m <- length(blocks)
  lo <- log(vy * 1e-7); hi <- log(vy * 1e3)
  core <- reml_core(y, X, blocks, list(), rep(log(vy / (m + 1)), m + 1), rep(lo, m + 1), rep(hi, m + 1),
                    maxit, ftol, const = FALSE)
  list(pred = reml_predict(core, X, Xte, blocks, list(), y), theta = core$eval$th, beta = core$eval$b,
       converged = core$fit$convergence == 0)
}

# Python benchmark_v4.reml_mixed (PSD components + free-sign cross terms; maxiter 200, ftol 1e-8)
reml_mixed <- function(y, blocks, signed = list(), maxit = 200, ftol = 1e-8) {
  n <- length(y); X <- matrix(1, n, 1); Xte <- matrix(1, nrow(blocks[[1]]$ut), 1)
  vy <- max(mean((y - mean(y))^2), 1e-8); m <- length(blocks); s <- length(signed)
  lo <- log(vy * 1e-7); hi <- log(vy * 1e3)
  core <- reml_core(y, X, blocks, signed, c(rep(log(vy / (m + 1)), m + 1), rep(0, s)),
                    c(rep(lo, m + 1), rep(-10 * vy, s)), c(rep(hi, m + 1), rep(10 * vy, s)), maxit, ftol)
  list(pred = reml_predict(core, X, Xte, blocks, signed, y), theta = core$eval$th, gamma = core$eval$gam,
       converged = core$fit$convergence == 0, objective = core$fit$value)
}

# Python part2_reml.reml_robust: several starting values for gamma plus a warm start from the
# model without cross terms; keep the best REML likelihood; "converged" requires that the
# optimiser succeeded and moved away from its starting point.
spectral_scale <- function(Stt) max(abs(eigen((Stt + t(Stt)) / 2, symmetric = TRUE, only.values = TRUE)$values))

reml_robust <- function(y, blocks, signed = list(), maxit = 300, gamma_starts = c(0, 0.25, -0.25, 0.5, -0.5)) {
  n <- length(y); X <- matrix(1, n, 1); Xte <- matrix(1, nrow(blocks[[1]]$ut), 1)
  vy <- max(mean((y - mean(y))^2), 1e-8); m <- length(blocks); s <- length(signed)
  lo <- log(vy * 1e-7); hi <- log(vy * 1e3)
  lower <- c(rep(lo, m + 1), rep(-5 * vy, s)); upper <- c(rep(hi, m + 1), rep(5 * vy, s))
  th0 <- rep(log(vy / (m + 1)), m + 1)
  starts <- lapply(if (s) gamma_starts else 0, function(g) c(th0, rep(g * vy, s)))
  if (s) {
    base <- reml_core(y, X, blocks, list(), th0, lower[1:(m + 1)], upper[1:(m + 1)], maxit, 1e-10)
    starts[[length(starts) + 1]] <- c(base$fit$par, rep(0, s))
  }
  best <- NULL; best_ok <- FALSE; n_ok <- 0
  for (x0 in starts) {
    core <- reml_core(y, X, blocks, signed, x0, lower, upper, maxit, 1e-10)
    if (is.null(core$eval)) next
    ok <- core$fit$convergence == 0 && max(abs(core$fit$par - x0)) > 1e-6
    n_ok <- n_ok + ok
    if (is.null(best) || core$fit$value < best$fit$value - 1e-9) { best <- core; best_ok <- ok }
  }
  list(pred = reml_predict(best, X, Xte, blocks, signed, y), theta = best$eval$th, gamma = best$eval$gam,
       converged = best_ok, nll = best$fit$value, starts_ok = n_ok)
}

# Python benchmark_v4.signed_blocks: raw cross-term blocks scaled on the training records
signed_blocks <- function(S, KE, KG, tr, te) {
  scale <- sqrt(mean(diag(KE)[tr]) * mean(diag(KG)[tr]))
  Stt <- S[tr, tr] / scale
  list(tt = (Stt + t(Stt)) / 2, ut = S[te, tr, drop = FALSE] / scale)
}

# ---- metrics (Python benchmark_ulte.basic_metrics / evaluate) -------------------
MIN_RECORDS_PER_ENV <- 10
basic_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  sdp <- function(x) sqrt(mean((x - mean(x))^2))
  cor_ <- if (length(obs) > 2 && sdp(obs) > 0 && sdp(pred) > 0) cor(obs, pred) else NA_real_
  rmse <- if (length(obs)) sqrt(mean((obs - pred)^2)) else NA_real_
  nrmse <- if (length(obs) > 1 && sd(obs) > 0) rmse / sd(obs) else NA_real_
  k <- max(1, ceiling(0.10 * length(obs)))
  top_o <- order(-obs, method = "radix")[1:k]    # stable sort, as numpy mergesort
  top_p <- order(-pred, method = "radix")[1:k]
  c(COR = cor_, RMSE = rmse, NRMSE = nrmse, PM10 = 100 * length(intersect(top_o, top_p)) / k)
}

evaluate <- function(obs, pred, env) {
  pooled <- basic_metrics(obs, pred)
  per <- lapply(sortc(env), function(e) {
    m <- env == e
    if (sum(m) >= MIN_RECORDS_PER_ENV) c(n = sum(m), basic_metrics(obs[m], pred[m])) else NULL
  })
  per <- do.call(rbind, per[!vapply(per, is.null, TRUE)])
  out <- c(COR_pooled = pooled[["COR"]], RMSE_pooled = pooled[["RMSE"]],
           NRMSE_pooled = pooled[["NRMSE"]], PM10_pooled = pooled[["PM10"]])
  if (!is.null(per)) {
    out <- c(out, COR_within = mean(per[, "COR"]), RMSE_within = mean(per[, "RMSE"]),
             NRMSE_within = mean(per[, "NRMSE"]), PM10_within = mean(per[, "PM10"]),
             environments_scored = nrow(per))
  } else {
    out <- c(out, COR_within = NA, RMSE_within = NA, NRMSE_within = NA, PM10_within = NA, environments_scored = 0)
  }
  out
}

# Python benchmark_ulte.paired_differences: same-partition differences against a reference model
paired_differences <- function(metrics, reference = "B_Hadamard", metric_cols = c("COR_within", "PM10_within", "COR_pooled")) {
  rows <- list()
  for (ds in unique(metrics$dataset)) {
    g <- metrics[metrics$dataset == ds, ]
    ref <- g[g$model == reference, ]
    for (mod in setdiff(unique(g$model), reference)) {
      gm <- g[g$model == mod, ]
      common <- intersect(ref$partition, gm$partition)
      for (mc in metric_cols) {
        if (!mc %in% names(g)) next
        dif <- gm[[mc]][match(common, gm$partition)] - ref[[mc]][match(common, ref$partition)]
        k <- length(dif); sdv <- if (k > 1) sd(dif) else NA
        rows[[length(rows) + 1]] <- data.frame(dataset = ds, model = mod, metric = mc, reference = reference,
          partitions = k, mean_difference = mean(dif), sd_difference = sdv,
          se_difference = if (k > 1) sdv / sqrt(k) else NA, partitions_better = sum(dif > 0))
      }
    }
  }
  do.call(rbind, rows)
}

# ---- splits exported from the Python code (splits/ folder) ----------------------
read_split <- function(file) read.csv(xzfile(file.path(SPLIT_DIR, file)), stringsAsFactors = FALSE)
test_mask_5050 <- function(parts, dataset, partition, gids) {
  tg <- parts$GID[parts$dataset == dataset & parts$partition == partition]
  gids %in% tg
}
