# =============================================================================
# 01_data_kernels.R - data access, kernels and matrix-product components
#
# R translation of code/python/part1/code/analysis_pipeline.py (load_dataset,
# genomic_matrix, environmental_matrix, record_kernels, direction_mask,
# build_components) and of the kernels used in Part II.
# Each function names the Python function it reproduces.
# Notation follows Part I, Section 2:  K_G = Z_G G Z_G',  K_E = Z_E E Z_E',  M = K_E K_G.
# =============================================================================

EPS <- .Machine$double.eps
DEGENERACY_TOL <- 1e-10                        # ||K_E K_G|| / (||K_E|| ||K_G||) below this => product = 0
PARTIAL_DIRECTION <- list(EYT_2 = list(c("Bed2IR", "Bed5IR")),   # (earlier, later) environment pairs
                          EYT_3 = list(c("Bed2IR", "Bed5IR")))
H_FRACTION <- 0.05                             # bandwidth of the identifier kernel (fraction of lines)

sortc <- function(x) sort(unique(as.character(x)), method = "radix")   # C-locale order, as Python sorted()

# ---- environments: similarity order and declared direction ------------------
is_chronological <- function(name) name %in% c("Indica", "Japonica") || startsWith(name, "Wheat_")

environment_time <- function(label) {           # '2010' -> 2010, 'YT_13_14' -> 2013
  y <- as.integer(regmatches(label, regexpr("[0-9]+", label)))
  ifelse(y < 100, y + 2000L, y)
}

ordered_env_levels <- function(name, labels) {  # deterministic, never depends on record order
  lev <- sortc(labels)
  if (is_chronological(name)) lev <- lev[order(environment_time(lev), lev, method = "radix")]
  lev
}

direction_mask <- function(name, rec_env) {     # B[i, j] = 1 if env of record i precedes env of record j
  rec_env <- as.character(rec_env); n <- length(rec_env)
  if (is_chronological(name)) {
    t <- environment_time(rec_env)
    return(outer(t, t, "<") * 1)
  }
  B <- matrix(0, n, n)
  for (pr in PARTIAL_DIRECTION[[name]]) B[outer(rec_env == pr[1], rec_env == pr[2], "&")] <- 1
  B
}

# ---- data --------------------------------------------------------------------
# Python: analysis_pipeline.load_dataset (reads the compact package instead of the .RData files)
load_ge <- function(name, data_dir = DATA_DIR) {
  x <- load_dataset(name, data_dir)
  ph <- x$pheno
  gid_col <- if ("Line" %in% names(ph)) "Line" else "GID"
  names(ph)[names(ph) == gid_col] <- "GID"
  ph$GID <- as.character(ph$GID); ph$Env <- as.character(ph$Env)
  rownames(ph) <- NULL
  traits <- setdiff(names(ph), c("GID", "Env", "Hybrid_Name"))
  traits <- traits[vapply(ph[traits], is.numeric, TRUE)]
  G <- NULL; G_ids <- NULL; markers <- NULL; marker_ids <- NULL
  if (name %in% c("EYT_1", "EYT_2", "EYT_3") || startsWith(name, "Wheat_")) {
    G <- x$G
    G_ids <- if (!is.null(x$geno_line_order)) x$geno_line_order else rownames(G)
    if (setequal(G_ids, colnames(G))) G <- G[, match(G_ids, colnames(G)), drop = FALSE]
  }
  if (!is.null(x$markers)) { markers <- x$markers; marker_ids <- rownames(x$markers) }
  EC <- x$EC
  if (!is.null(EC)) { EC <- as.data.frame(EC); EC[[1]] <- as.character(EC[[1]]) }
  list(name = name, pheno = ph, traits = traits, G = G, G_ids = G_ids,
       markers = markers, marker_ids = marker_ids, EC = EC)
}

# Python: analysis_pipeline.genomic_matrix - supplied G, or centred markers / sum of variances;
# symmetrised and scaled to unit mean diagonal.
genomic_matrix <- function(d) {
  if (!is.null(d$G)) {
    G <- unname(d$G); ids <- as.character(d$G_ids)
  } else {
    X <- d$markers; ids <- as.character(d$marker_ids)
    X <- sweep(X, 2, colMeans(X))
    denom <- sum(colMeans(X^2))                  # sum of marker variances (ddof = 0)
    G <- tcrossprod(X) / denom
  }
  G <- (G + t(G)) / 2
  G <- G / mean(diag(G))
  dimnames(G) <- list(ids, ids)
  G
}

# Python: analysis_pipeline.environmental_matrix - standardised covariates, E = W W'/q,
# unit mean diagonal; identity when there are no covariates.
environmental_matrix <- function(d, env_levels) {
  if (is.null(d$EC)) return(diag(length(env_levels)))
  ec <- d$EC; labels <- as.character(ec[[1]])
  W <- as.matrix(ec[, -1, drop = FALSE]); storage.mode(W) <- "double"
  W <- sweep(W, 2, colMeans(W))
  sd0 <- sqrt(colMeans(W^2)); keep <- sd0 > 1e-12
  W <- sweep(W[, keep, drop = FALSE], 2, sd0[keep], "/")
  E <- tcrossprod(W) / ncol(W)
  o <- match(env_levels, labels)
  E <- E[o, o, drop = FALSE]
  E / mean(diag(E))
}

# Python: analysis_pipeline.record_kernels - record-level K_G and K_E for the records of `ph`.
record_kernels <- function(d, ph = d$pheno) {
  G <- genomic_matrix(d)
  gi <- match(ph$GID, rownames(G))
  KG <- unname(G[gi, gi])
  rec_env <- ph$Env
  env_levels <- ordered_env_levels(d$name, rec_env)
  ei <- match(rec_env, env_levels)
  E <- environmental_matrix(d, env_levels)
  KE <- E[ei, ei]
  list(KG = KG, KE = KE, rec_env = rec_env, env_levels = env_levels)
}

product_relative_norm <- function(M, KE, KG) norm(M, "F") / max(norm(KE, "F") * norm(KG, "F"), EPS)

triangles <- function(M) {                       # C0 = D + U + U',  P0 = D + L + L'  (positional)
  D <- diag(diag(M)); U <- M; U[lower.tri(U, diag = TRUE)] <- 0
  L <- M; L[upper.tri(L, diag = TRUE)] <- 0
  list(C = D + U + t(U), P = D + L + t(L))
}

# Python: analysis_pipeline.build_components - all raw (pre-projection) record-level kernels.
build_components <- function(d, ph = d$pheno, include_original = TRUE) {
  rk <- record_kernels(d, ph)
  KG <- rk$KG; KE <- rk$KE
  M <- KE %*% KG
  rel <- product_relative_norm(M, KE, KG)
  degenerate <- rel < DEGENERACY_TOL
  comps <- list(E = KE, G = KG, H = KE * KG)
  B <- direction_mask(d$name, rk$rec_env)
  directed <- any(B != 0)
  if (!degenerate) {
    if (include_original) { tr <- triangles(M); comps$C0 <- tr$C; comps$P0 <- tr$P }
    if (directed) {
      Hm <- 1 - B - t(B)
      comps$CSP <- B * M + t(B) * t(M)
      comps$PSP <- B * t(M) + t(B) * M
      comps$TSP <- 0.5 * Hm * (M + t(M))
    } else {
      comps$SGE <- 0.5 * (M + t(M))
    }
  }
  list(comps = comps, info = list(product_relative_norm = rel, degenerate = degenerate,
                                  directed = directed, directional_pairs = sum(B),
                                  env_levels = rk$env_levels))
}

# ---- line identifiers --------------------------------------------------------
id_number <- function(gids) {                    # Python id_order_test.id_number: first run of digits
  as.numeric(regmatches(gids, regexpr("[0-9]+", gids)))
}
id_digits <- function(gids) {                    # Python "".join(digits): all digits concatenated
  v <- gsub("[^0-9]", "", gids); v[v == ""] <- "0"; as.numeric(v)
}

# Identifier-proximity kernel K_id(j, k) = exp(-|rank_j - rank_k| / h), h = 5% of the lines,
# record level (Python id_order_test.run_dataset / part2_id_split.run_dataset).
identifier_kernel <- function(gids, num = id_number(gids), h_min = 1) {
  lines <- sortc(gids)
  line_num <- num[match(lines, gids)]
  rk <- rank(line_num)                           # average ranks for ties, as pandas
  r <- rk[match(gids, lines)]
  h <- max(h_min, H_FRACTION * length(lines))
  exp(-abs(outer(r, r, "-")) / h)
}

# Python part2_id_split.gaussian_kernel: exp(-d / median d), d_ij = K_ii + K_jj - 2 K_ij
gaussian_kernel <- function(KG) {
  dg <- diag(KG)
  D <- pmax(outer(dg, dg, "+") - 2 * KG, 0)
  exp(-D / median(D[D > 1e-10]))
}

sqrtm_psd <- function(A) {                       # symmetric PSD square root
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}
