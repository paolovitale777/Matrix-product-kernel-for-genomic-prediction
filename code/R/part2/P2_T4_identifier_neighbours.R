# =============================================================================
# Part II - Table 4 (Section 3.2): are identifier neighbours relatives, and is their phenotypic
# similarity shared across environments or environment-specific?  (model-free diagnostics)
# Python equivalents: part2_neighbour_env.py (Cuevas data), part2_neighbour_env_ge.py (G x E
# data), part2_id_relatedness.py (column "Related").
#
# Lines are ordered by the numeric part of their identifier. For standardised phenotypes:
#   Within (lag l):  correlation of lines l positions apart in the SAME environment
#   Ratio  (lag l):  the same across two environments, divided by the across-environment
#                    correlation of the same line ("Same line")
#   Related:         fraction of adjacent pairs among the 5% most related pairs (genomic G)
# Output: P2_T4_identifier_neighbours.csv and P2_T4_identifier_relatedness.csv
# Usage:  Rscript code/R/part2/P2_T4_identifier_neighbours.R
# =============================================================================
source("code/R/00_setup.R")
corr <- function(a, b) { ok <- !is.na(a) & !is.na(b); if (sum(ok) > 20) cor(a[ok], b[ok]) else NA }
stdz <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)      # pandas .std() (ddof = 1)
stdz0 <- function(v) (v - mean(v, na.rm = TRUE)) / sqrt(mean((v - mean(v, na.rm = TRUE))^2, na.rm = TRUE))  # numpy (ddof = 0)
neighbour_stats <- function(Y, lags) {
  envs <- names(Y); pairs <- combn(envs, 2, simplify = FALSE)
  same <- mean(sapply(pairs, function(p) corr(Y[[p[1]]], Y[[p[2]]])), na.rm = TRUE)
  do.call(rbind, lapply(lags, function(l) {
    nb <- function(a, b) corr(a[1:(length(a) - l)], b[(1 + l):length(b)])
    w <- mean(sapply(envs, function(e) nb(Y[[e]], Y[[e]])), na.rm = TRUE)
    cr <- mean(sapply(pairs, function(p) mean(c(nb(Y[[p[1]]], Y[[p[2]]]), nb(Y[[p[2]]], Y[[p[1]]])), na.rm = TRUE)), na.rm = TRUE)
    data.frame(lag = l, nb_within = w, nb_cross = cr, same_line_cross_env = same, ratio = cr / same)
  }))
}

rows <- list()
# --- Cuevas et al. (2025) data: one vector per environment, lines in file order
for (nm in c("Wheat_599", "EYT_16_17", "EYT_22_23", "EYT_23_24")) {
  if (!file.exists(file.path(DATA_DIR, "Cuevas2025", paste0(nm, ".rds")))) next
  x <- load_dataset(nm); o <- order(id_digits(x$pheno$Line), method = "radix")
  Y <- lapply(setNames(x$env_names, x$env_names), function(e) stdz0(x$pheno[[e]])[o])
  Y <- Y[vapply(Y, function(v) sum(!is.na(v)) > 50, TRUE)]
  if (length(Y) < 2) next
  rows[[nm]] <- data.frame(dataset = nm, neighbour_stats(Y, c(1, 2, 5, 10, 50)))
}
# --- G x E data: line x environment table of grain yield
for (nm in GE_DATASETS) {
  d <- load_ge(nm); ph <- d$pheno
  W <- tapply(ph$GY, list(ph$GID, ph$Env), mean)
  W <- W[order(id_digits(rownames(W)), method = "radix"), , drop = FALSE]
  Y <- lapply(setNames(colnames(W), colnames(W)), function(e) stdz(W[, e]))
  rows[[nm]] <- data.frame(dataset = nm, neighbour_stats(Y, c(1, 5, 50)))
}
nb <- do.call(rbind, rows); rownames(nb) <- NULL
write_out(nb, "P2_T4_identifier_neighbours.csv")

# --- Related: genomic relationship of adjacent identifiers (part2_id_relatedness.py)
rel <- list()
for (nm in GE_DATASETS) {
  d <- load_ge(nm); G <- genomic_matrix(d)
  lines <- sortc(d$pheno$GID); Gl <- G[lines, lines]
  o <- order(id_number(lines), method = "radix")
  adj <- Gl[cbind(o[-length(o)], o[-1])]
  allp <- Gl[upper.tri(Gl)]
  rel[[nm]] <- data.frame(dataset = nm, lines = length(lines), G_adjacent = mean(adj),
                          adjacent_in_top5pct = mean(adj >= quantile(allp, 0.95, type = 7)))
}
rel <- do.call(rbind, rel)
write_out(rel, "P2_T4_identifier_relatedness.csv")
t4 <- merge(rel[, c("dataset", "adjacent_in_top5pct")],
            reshape(nb[nb$lag %in% c(1, 5), c("dataset", "lag", "nb_within", "ratio", "same_line_cross_env")],
                    idvar = c("dataset", "same_line_cross_env"), timevar = "lag", direction = "wide"), all = TRUE)
print(t4, digits = 2)
