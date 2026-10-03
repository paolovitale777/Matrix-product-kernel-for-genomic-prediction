# =============================================================================
# Part I - Table 1: datasets, lines in all environments, declared direction and numerical
# rank of the record-level matrix product M_GE = K_E K_G (Proposition 1).
# Python equivalent: analysis_pipeline.py qc  +  matrix_rank.py
# Output: code/R/outputs/P1_T1_datasets_rank.csv
# Usage:  Rscript code/R/part1/P1_T1_datasets_rank.R [--datasets ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", GE_DATASETS)
rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  rk <- record_kernels(d)
  # rank(M) = rank(E N G): M = Z_E (E N G) Z_G' and Z_E, Z_G have full column rank (Section 2.3),
  # so the rank is computed from the small environment x line matrix (fast and exact).
  G <- genomic_matrix(d); lines <- sortc(ph$GID)
  N <- table(factor(ph$Env, levels = rk$env_levels), factor(ph$GID, levels = lines))
  ENG <- environmental_matrix(d, rk$env_levels) %*% unclass(N) %*% G[lines, lines]
  s <- svd(ENG, nu = 0, nv = 0)$d
  rank_M <- if (s[1] < 1e-12 * max(abs(ENG), 1e-300) || max(abs(ENG)) < 1e-10 * nrow(ph)) 0 else sum(s > s[1] * 1e-10)
  M <- rk$KE %*% rk$KG
  per_line <- tapply(ph$Env, ph$GID, function(e) length(unique(e)))
  B <- direction_mask(nm, ph$Env)
  rows[[length(rows) + 1]] <- data.frame(
    dataset = nm, records = nrow(ph), lines = length(unique(ph$GID)),
    environments = length(unique(ph$Env)),
    environment_labels = paste(ordered_env_levels(nm, ph$Env), collapse = "|"),
    lines_in_all_env = sum(per_line == length(unique(ph$Env))),
    duplicate_line_env = sum(duplicated(ph[c("GID", "Env")])),
    missing_GY = sum(is.na(ph$GY)),
    genomic_source = if (!is.null(d$G)) "supplied G" else "centred markers",
    declared_direction = if (is_chronological(nm)) "calendar year" else if (any(B > 0)) "Bed2IR < Bed5IR" else "none",
    product_relative_norm = product_relative_norm(M, rk$KE, rk$KG),
    rank_M = rank_M, second_over_first_singular_value = if (length(s) > 1 && s[1] > 0) s[2] / s[1] else NA)
  log_msg(nm, "rank =", rank_M)
}
res <- do.call(rbind, rows)
write_out(res, "P1_T1_datasets_rank.csv")
print(res[, c("dataset", "environments", "lines", "records", "lines_in_all_env", "declared_direction", "genomic_source", "rank_M")])
