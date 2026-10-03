# rebuild_original_files.R - recreate the original .RData files from the compact package.
#
# The published Python code (code/python/part1 and code/python/part2) reads the data in their
# original format. Run, from the repository root:
#     Rscript data/rebuild_original_files.R
# This writes New_GE_Framework/<name>.RData (ten G x E datasets) and cuevas_data/<name>.Rdata
# (data of Cuevas et al. 2025) into code/python/part1 and code/python/part2, with the same
# objects, names, values and order as the files used for the papers (verified object by object;
# the Gaussian kernel K of the Cuevas data is recomputed from G and agrees to ~1e-16).
# The R code does not need these files: it reads data/ directly.

source("data/load_dataset.R")

write_ge <- function(name, dest) {
  d <- load_dataset(name)
  e <- new.env()
  if (name %in% c("EYT_1", "EYT_2", "EYT_3")) {
    e$Pheno <- d$pheno
    e$Geno  <- data.frame(Line = d$geno_line_order, d$G, check.names = FALSE)
    attr(e$Geno, "row.names") <- d$geno_row_names
    objs <- c("Pheno", "Geno")
  } else if (name %in% c("Indica", "Japonica")) {
    e$Markers <- data.frame(Line = rownames(d$markers), d$markers, check.names = FALSE)
    attr(e$Markers, "row.names") <- d$markers_row_names
    e$EC <- d$EC
    e$Pheno <- d$pheno
    objs <- c("Markers", "EC", "Pheno")
  } else {
    e$dat_ls <- list(Pheno = d$pheno, Markers = d$markers, Geno = d$G)
    objs <- "dat_ls"
  }
  for (dd in dest) save(list = objs, envir = e, file = file.path(dd, "New_GE_Framework", paste0(name, ".RData")))
}

write_cuevas <- function(name, dest) {
  d <- load_dataset(name)
  e <- new.env()
  e$A <- d$A; e$G <- d$G; e$K <- d$K
  for (v in d$env_names) e[[v]] <- setNames(d$pheno[[v]], d$pheno_vector_names)
  for (dd in dest) save(list = d$object_order, envir = e, file = file.path(dd, "cuevas_data", paste0(name, ".Rdata")))
}

dest <- file.path("code", "python", c("part1", "part2"))
for (dd in dest) for (sub in c("New_GE_Framework", "cuevas_data")) dir.create(file.path(dd, sub), recursive = TRUE, showWarnings = FALSE)
for (nm in DATASETS_GE) { write_ge(nm, dest); cat("written", nm, "\n") }
for (nm in DATASETS_CUEVAS) {
  if (file.exists(file.path("data", "Cuevas2025", paste0(nm, ".rds")))) { write_cuevas(nm, dest); cat("written", nm, "\n") }
  else cat("not included:", nm, "(download it from CIMMYT and add it with data/add_cuevas_file.R)\n")
}
