# =============================================================================
# Part I - Table 2 (Section 3.1): matrix-level permutation test.
# Python equivalent: analysis_pipeline.py order --permutations 20
#
# For each of 20 record permutations (splits/permutation_test.csv.xz):
#   original components C0, P0: the records are reordered, the triangles are rebuilt, mapped
#     back to the file order and compared with the reference triangles;
#   proposed components (C_SP, P_SP, T_SP or S_GE): ALL kernels and the direction mask are
#     rebuilt from the permuted records, mapped back and compared.
# delta = ||Q' A* Q - A||_F / ||A||_F  (Equation 21). An equivariant construction gives ~1e-16.
# Output: P1_T2_permutation_test.csv (per permutation) and P1_T2_permutation_summary.csv.
# Note: rebuilding the kernels for every permutation is slow for the EYT datasets
# (n = 3,064-4,820 records); use --datasets to run them one at a time.
# Usage:  Rscript code/R/part1/P1_T2_permutation_test.R [--datasets ...] [--permutations 20]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", GE_DATASETS)
n_perm <- as.integer(arg_values("--permutations", 20))
perms <- read_split("permutation_test.csv.xz")
rel_fro <- function(A, ref) norm(A - ref, "F") / max(norm(ref, "F"), EPS)

raw <- list(); summ <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno; n <- nrow(ph)
  rk <- record_kernels(d); M <- rk$KE %*% rk$KG
  rel <- product_relative_norm(M, rk$KE, rk$KG)
  row <- data.frame(dataset = nm, records = n, product_relative_norm = rel, degenerate = rel < DEGENERACY_TOL)
  if (rel < DEGENERACY_TOL) { log_msg(nm, "degenerate product, skipped"); summ[[nm]] <- row; next }
  ref_tri <- triangles(M)
  ref <- build_components(d, include_original = FALSE)
  keys <- intersect(c("CSP", "PSP", "TSP", "SGE"), names(ref$comps))
  B_ref <- direction_mask(nm, ph$Env)
  res <- matrix(NA, n_perm, 3 + length(keys), dimnames = list(NULL, c("C0_delta", "P0_delta", "mask_max_abs_change", paste0(keys, "_delta"))))
  for (r in seq_len(n_perm)) {
    # original components (Python original_permutation_deltas; computed here directly)
    p <- perms$record[perms$dataset == nm & perms$test == "original" & perms$permutation == r]
    inv <- order(p)
    tri_p <- triangles(M[p, p])
    res[r, "C0_delta"] <- rel_fro(tri_p$C[inv, inv], ref_tri$C)
    res[r, "P0_delta"] <- rel_fro(tri_p$P[inv, inv], ref_tri$P)
    # proposed components: rebuild everything from the permuted records
    p <- perms$record[perms$dataset == nm & perms$test == "proposed" & perms$permutation == r]
    inv <- order(p)
    ph_p <- ph[p, ]; rownames(ph_p) <- NULL
    res[r, "mask_max_abs_change"] <- max(abs(direction_mask(nm, ph_p$Env)[inv, inv] - B_ref))
    comp_p <- build_components(d, ph_p, include_original = FALSE)$comps
    for (k in keys) res[r, paste0(k, "_delta")] <- rel_fro(comp_p[[k]][inv, inv], ref$comps[[k]])
    log_msg(nm, "permutation", r)
  }
  raw[[nm]] <- data.frame(dataset = nm, permutation = seq_len(n_perm), res)
  row$C0_delta_median <- median(res[, "C0_delta"]); row$C0_delta_max <- max(res[, "C0_delta"])
  row$P0_delta_median <- median(res[, "P0_delta"]); row$P0_delta_max <- max(res[, "P0_delta"])
  row$proposed_delta_max <- max(res[, paste0(keys, "_delta")])
  row$mask_max_abs_change <- max(res[, "mask_max_abs_change"])
  if (ref$info$directed) row$reconstruction_delta <- rel_fro(0.5 * (ref$comps$CSP + ref$comps$PSP) + ref$comps$TSP, 0.5 * (M + t(M)))
  summ[[nm]] <- row
  write_out(do.call(rbind, raw), "P1_T2_permutation_test.csv")
}
all_cols <- unique(unlist(lapply(summ, names)))
summ <- do.call(rbind, lapply(summ, function(x) { x[setdiff(all_cols, names(x))] <- NA; x[all_cols] }))
write_out(summ, "P1_T2_permutation_summary.csv")
print(summ, digits = 3)
