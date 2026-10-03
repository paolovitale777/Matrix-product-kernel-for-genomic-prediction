# =============================================================================
# Part I - Table 5 (Section 3.6): prediction of the latest calendar year from all earlier
# years, for all target records, lines already phenotyped in training and new lines.
# Python equivalent: prospective_gy.py (wheat)  and  prospective_gy.py --datasets Indica
#                    Japonica --ec previous (rice: covariates of the target year replaced by
#                    those of the previous year, because the weather of a future year is unknown)
# Models: A main effects, B Hadamard, C/D/E/F triangular variants, proposed (MD/MHD with a
#         declared direction).  Table 5 reports A, B, F (= F_original_HCP) and Proposed
#         (= MHD_proposed_HCPT).
# Output: P1_T5_prospective_year.csv
# Usage:  Rscript code/R/part1/P1_T5_prospective_year.R [--datasets ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", c("Indica", "Japonica", "Wheat_1", "Wheat_2", "Wheat_3", "Wheat_5", "Wheat_6"))
MIN_SUBSET <- 10

model_specs <- function(info) {                     # Python benchmark_ulte.model_specs
  m <- list(A_main = c("E", "G"), B_Hadamard = c("E", "G", "H"))
  if (info$degenerate) return(m)
  m <- c(m, list(C_original_C = c("E", "G", "C0"), D_original_P = c("E", "G", "P0"),
                 E_original_CP = c("E", "G", "C0", "P0"), F_original_HCP = c("E", "G", "H", "C0", "P0")))
  if (info$directed) c(m, list(MD_proposed_CPT = c("E", "G", "CSP", "PSP", "TSP"), MHD_proposed_HCPT = c("E", "G", "H", "CSP", "PSP", "TSP")))
  else c(m, list(MS_proposed_S = c("E", "G", "SGE"), MHS_proposed_HS = c("E", "G", "H", "SGE")))
}
ec_previous_year <- function(d, target_envs) {      # Python prospective_gy.ec_previous_year
  ec <- d$EC; lab <- as.character(ec[[1]])
  is_t <- lab %in% target_envs; yrs <- environment_time(lab)
  prev <- !is_t & yrs == max(yrs[!is_t])
  for (j in 2:ncol(ec)) ec[[j]][is_t] <- mean(as.numeric(ec[[j]][prev]))
  d$EC <- ec; d
}

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph <- d$pheno
  years <- environment_time(ph$Env)
  is_test <- years == max(years)
  if (!is.null(d$EC)) d <- ec_previous_year(d, sortc(ph$Env[is_test]))
  te <- which(is_test); tr <- which(!is_test)
  tested <- ph$GID[te] %in% ph$GID[tr]
  bc <- build_components(d)
  models <- model_specs(bc$info)
  needed <- unique(unlist(models))
  blocks <- Filter(Negate(is.null), setNames(lapply(needed, function(k) psd_blocks(bc$comps[[k]], tr, te)), needed))
  for (mod in names(models)) {
    fit <- reml_fit_predict(ph$GY[tr], blocks[intersect(models[[mod]], names(blocks))])
    for (sub in c("all", "tested_lines", "new_lines")) {
      m <- switch(sub, all = rep(TRUE, length(te)), tested_lines = tested, new_lines = !tested)
      if (sum(m) < MIN_SUBSET) next
      ev <- evaluate(ph$GY[te][m], fit$pred[m], ph$Env[te][m])
      rows[[length(rows) + 1]] <- data.frame(dataset = nm, target = paste(sortc(ph$Env[te]), collapse = "|"),
        subset = sub, n = sum(m), model = mod, ec = if (is.null(d$EC)) "none" else "previous_year",
        converged = fit$converged, t(ev))
    }
    log_msg(nm, mod, "done")
  }
  write_out(do.call(rbind, rows), "P1_T5_prospective_year.csv")
}
res <- do.call(rbind, rows)
r5 <- res[res$model %in% c("A_main", "B_Hadamard", "F_original_HCP", "MHD_proposed_HCPT"), ]
print(reshape(r5[, c("dataset", "subset", "n", "model", "COR_within")], idvar = c("dataset", "subset", "n"),
              timevar = "model", direction = "wide"), digits = 3)
