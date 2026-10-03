# =============================================================================
# Part II - Table 7 (breeding scenarios), Supplementary Table S1 (rice: three versions of the
# environmental covariates of the target year) and PM10 values for Supplementary Table S2.
# Python equivalent: part2/code/part2_scenarios.py
#   --scenarios forward  (Wheat_1-6; Indica, Japonica with --ec previous / observed / train_mean)
#   --scenarios gradient (EYT_2, EYT_3: Bed5IR predicted from Bed2IR and heat/drought trials)
#
# forward : the latest calendar year is predicted from all earlier years (all target records,
#           lines already tested in training, new lines).
# gradient: Bed5IR is predicted from the stress environments, for lines tested there
#           ("tested_lines") and for new lines (five 50:50 partitions; half of the lines are
#           removed from training entirely).
# Models M1-M5 as in P2_T5_T6_identifier_split.R.
# Output: P2_T7_S1_breeding_scenarios.csv (column ec: observed / previous / train_mean)
#   Table 7 uses ec = "previous" for the rice datasets and "observed" for the others.
# Usage:  Rscript code/R/part2/P2_T7_S1_breeding_scenarios.R [--datasets ...]
# =============================================================================
source("code/R/00_setup.R")
datasets <- arg_values("--datasets", c("Indica", "Japonica", "Wheat_1", "Wheat_2", "Wheat_3", "Wheat_5", "Wheat_6", "EYT_2", "EYT_3"))
parts <- read_split("partitions_5050.csv.xz")
MODELS <- list(B_Hadamard = c("E", "G", "H"), GK_Hadamard = c("E", "GK", "HK"), ID_common = c("E", "G", "H", "ID"),
               ID_env = c("E", "G", "H", "ID_E"), GK_ID_env = c("E", "GK", "HK", "ID_E"))
GRADIENT <- list(EYT_2 = list(train = c("Bed2IR", "EHT", "LHT"), target = "Bed5IR"),
                 EYT_3 = list(train = c("Bed2IR", "LHT"), target = "Bed5IR"))
MIN_SUBSET <- 10

ec_replace <- function(d, target_envs, mode) {       # Python part2_scenarios.ec_previous_year / ec_training_mean
  ec <- d$EC; lab <- as.character(ec[[1]]); is_t <- lab %in% target_envs
  src <- if (mode == "previous") { yrs <- environment_time(lab); !is_t & yrs == max(yrs[!is_t]) } else !is_t
  for (j in 2:ncol(ec)) ec[[j]][is_t] <- mean(as.numeric(ec[[j]][src]))
  d$EC <- ec; d
}

rows <- list()
for (nm in datasets) {
  d0 <- load_ge(nm); ph <- d0$pheno
  y <- ph$GY; envs <- ph$Env; gids <- ph$GID
  bc <- build_components(d0, include_original = FALSE)
  K_id <- identifier_kernel(gids, id_number(gids))
  S <- outer(envs, envs, "==") * 1
  GK <- gaussian_kernel(bc$comps$G)
  base <- list(E = bc$comps$E, G = bc$comps$G, H = bc$comps$H, ID = K_id, ID_E = K_id * S, GK = GK, HK = bc$comps$E * GK)
  X <- matrix(1, length(y), 1)
  if (nm %in% names(GRADIENT)) {
    g <- GRADIENT[[nm]]; in_tr <- envs %in% g$train; in_te <- envs %in% g$target
    sp <- c(list(list("tested_lines", which(in_tr), which(in_te))),
            lapply(1:5, function(fi) { nw <- test_mask_5050(parts, nm, fi, gids); list(paste0("new_lines_", fi), which(in_tr & !nw), which(in_te & nw)) }))
    scen <- "gradient"; ec_modes <- "observed"
  } else {
    years <- environment_time(envs); yr <- max(years)
    sp <- list(list(paste0("predict_", yr), which(years < yr), which(years == yr)))
    scen <- "forward"; ec_modes <- if (is.null(d0$EC)) "observed" else c("observed", "previous", "train_mean")
  }
  for (ecm in ec_modes) {
    for (s in sp) {
      tr <- s[[2]]; te <- s[[3]]
      if (length(te) < MIN_SUBSET) next
      kernels <- base
      if (ecm != "observed") {
        b2 <- build_components(ec_replace(d0, sortc(envs[te]), ecm), include_original = FALSE)$comps
        kernels$E <- b2$E; kernels$H <- b2$H; kernels$HK <- b2$E * GK
      }
      tested <- gids[te] %in% gids[tr]
      blocks <- Filter(Negate(is.null), lapply(kernels, psd_blocks, tr = tr, te = te))
      for (mod in names(MODELS)) {
        used <- intersect(MODELS[[mod]], names(blocks))
        fit <- reml_fixed(y[tr], X[tr, , drop = FALSE], blocks[used], X[te, , drop = FALSE])
        for (sub in c("all", "tested", "new")) {
          m <- switch(sub, all = rep(TRUE, length(te)), tested = tested, new = !tested)
          if (sum(m) < MIN_SUBSET) next
          ev <- evaluate(y[te][m], fit$pred[m], envs[te][m])
          rows[[length(rows) + 1]] <- data.frame(dataset = nm, scenario = scen, ec = ecm, split = s[[1]],
            subset = sub, n = sum(m), model = mod, converged = fit$converged, t(ev))
        }
      }
      log_msg(nm, scen, ecm, s[[1]], "done")
    }
  }
  write_out(do.call(rbind, rows), "P2_T7_S1_breeding_scenarios.csv")
}
res <- do.call(rbind, rows)
res$target <- ifelse(res$scenario == "gradient", ifelse(startsWith(res$split, "new"), "Bed5IR new lines", "Bed5IR tested lines"), "next year")
t7 <- res[(res$ec == "observed" & !res$dataset %in% c("Indica", "Japonica") | res$ec == "previous") &
          ((res$scenario == "forward" & res$subset == "all") | (res$scenario == "gradient" & res$subset == "all")), ]
cat("\nTable 7 (COR_within):\n")
print(round(tapply(t7$COR_within, list(paste(t7$dataset, t7$target), t7$model), mean), 3))
s1 <- res[res$dataset %in% c("Indica", "Japonica") & res$subset == "all", ]
if (nrow(s1)) { cat("\nTable S1 (rice, COR_within by covariate version):\n"); print(round(tapply(s1$COR_within, list(paste(s1$dataset, s1$ec), s1$model), mean), 3)) }
