# =============================================================================
# Part I - Table 3, columns "F, file order" and "F, random orders" (Section 3.3).
# Python equivalent: order_accuracy_sensitivity.py
#   --datasets EYT_1 Wheat_1 --orderings 10   and   --datasets <others> --orderings 3
#
# The phenotype table is re-ordered at random (splits/record_orders.csv.xz; ordering 0 = the
# file order), ALL kernels are rebuilt, and the Hadamard model (B), the triangular model F and
# the proposed model are refitted on the same 50:50 cultivar partitions. Partitions are
# defined by cultivar, so only the positional C0/P0 can change between orderings.
# Output: P1_T3_record_order_accuracy.csv (dataset x ordering x partition x model)
#   Table 3 "F, file order"    = mean COR_within of F_original_HCP at ordering 0;
#   Table 3 "F, random orders" = mean COR_within of F_original_HCP over orderings >= 1.
# Usage:  Rscript code/R/part1/P1_T3_record_order_accuracy.R [--datasets ...] [--orderings 10]
#         (default: 10 orderings for EYT_1 and Wheat_1, 3 for the other datasets, as in the paper)
# =============================================================================
source("code/R/00_setup.R")
datasets   <- arg_values("--datasets", setdiff(GE_DATASETS, "Indica"))
n_ord_arg  <- arg_values("--orderings", NA)
partitions <- as.integer(arg_values("--partitions", 1:5))
parts  <- read_split("partitions_5050.csv.xz")
orders <- read_split("record_orders.csv.xz")

rows <- list()
for (nm in datasets) {
  d <- load_ge(nm); ph0 <- d$pheno
  n_ord <- if (!is.na(n_ord_arg[1])) as.integer(n_ord_arg) else if (nm %in% c("EYT_1", "Wheat_1")) 10 else 3
  for (o in 0:n_ord) {
    ph <- if (o == 0) ph0 else { p <- orders$record[orders$dataset == nm & orders$ordering == o]; x <- ph0[p, ]; rownames(x) <- NULL; x }
    bc <- build_components(d, ph)
    if (bc$info$degenerate) { log_msg(nm, "degenerate, skipped"); break }
    use <- list(B_Hadamard = c("E", "G", "H"), F_original_HCP = c("E", "G", "H", "C0", "P0"),
                proposed_full = if (bc$info$directed) c("E", "G", "H", "CSP", "PSP", "TSP") else c("E", "G", "H", "SGE"))
    for (fi in partitions) {
      is_te <- test_mask_5050(parts, nm, fi, ph$GID)
      tr <- which(!is_te); te <- which(is_te)
      needed <- unique(unlist(use))
      blocks <- Filter(Negate(is.null), setNames(lapply(needed, function(k) psd_blocks(bc$comps[[k]], tr, te)), needed))
      for (mod in names(use)) {
        fit <- reml_fit_predict(ph$GY[tr], blocks[intersect(use[[mod]], names(blocks))])
        ev <- evaluate(ph$GY[te], fit$pred, ph$Env[te])
        rows[[length(rows) + 1]] <- data.frame(dataset = nm, ordering = o, partition = fi, model = mod,
                                               converged = fit$converged, t(ev))
      }
      log_msg(nm, "ordering", o, "partition", fi, "done")
    }
    write_out(do.call(rbind, rows), "P1_T3_record_order_accuracy.csv")
  }
}
res <- do.call(rbind, rows)
f <- res[res$model == "F_original_HCP", ]
cat("\nTable 3: F in file order vs random orders (mean COR_within)\n")
print(round(cbind(file_order = tapply(f$COR_within[f$ordering == 0], f$dataset[f$ordering == 0], mean),
                  random_orders = tapply(f$COR_within[f$ordering > 0], f$dataset[f$ordering > 0], mean)), 3))
