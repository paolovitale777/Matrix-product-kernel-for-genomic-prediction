# add_cuevas_file.R - add a dataset of Cuevas et al. (2025) downloaded from the CIMMYT repository.
#
# EYT_23_24 is not redistributed in this package. To reproduce its results, download
# EYT_23_24.Rdata from the CIMMYT Research Data & Software Repository Network
# (https://hdl.handle.net/11529/10549160; access on request) and run, from the repository root:
#     Rscript data/add_cuevas_file.R path/to/EYT_23_24.Rdata
# The file is converted to data/Cuevas2025/EYT_23_24.rds (same format as the other datasets);
# the R scripts (P1_T6_cuevas_triangles.R, P2_T3_cross_covariance.R, P2_T4_identifier_neighbours.R)
# then include it automatically, and rebuild_original_files.R writes it for the Python code.
# Note: splits/cuevas_splits.csv.xz and cuevas_shuffles.csv.xz contain the splits of the three
# datasets in the package; for EYT_23_24 run code/python/export_splits.py again after adding it.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !file.exists(args[1])) stop("usage: Rscript data/add_cuevas_file.R path/to/<name>.Rdata")
f <- args[1]
nm <- sub("\\.[Rr][Dd]ata$", "", basename(f))
e <- new.env(); objs <- load(f, envir = e)
stopifnot(all(c("A", "G", "K") %in% objs))
envs <- setdiff(objs, c("A", "G", "K"))
lines <- rownames(e$A)
ph <- data.frame(Line = lines, sapply(envs, function(v) as.numeric(e[[v]])), check.names = FALSE)
if (all(is.na(ph[, -1]))) stop("the file contains no phenotypes - it is not a usable dataset")
d <- outer(diag(e$G), diag(e$G), "+") - 2 * e$G
if (max(abs(exp(-d / median(d)) - e$K)) > 1e-10) warning("K is not exp(-d/median(d)) of G; K is recomputed from G when loaded")
lower <- function(X) list(dimnames = dimnames(X), n = nrow(X), values = X[lower.tri(X, diag = TRUE)])
obj <- list(name = nm, crop = "wheat", pheno = ph, A_lower = lower(e$A), G_lower = lower(e$G), env_names = envs,
            source = "Cuevas et al. (2025); CIMMYT Research Data repository, hdl:11529/10549160 (added by the user)",
            notes = "Wide phenotype table (one column per environment). K = exp(-d/median(d)), d from G.",
            object_order = objs, pheno_vector_names = names(e[[envs[1]]]),
            dimnames_K = dimnames(e$K), dimnames_G = dimnames(e$G))
dir.create(file.path("data", "Cuevas2025"), showWarnings = FALSE, recursive = TRUE)
saveRDS(obj, file.path("data", "Cuevas2025", paste0(nm, ".rds")), compress = "xz")
write.csv(ph, file.path("data", "Cuevas2025", paste0(nm, "_pheno.csv")), row.names = FALSE)
cat("written data/Cuevas2025/", nm, ".rds (", length(lines), " lines, ", length(envs), " environments)\n", sep = "")
