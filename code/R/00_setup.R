# =============================================================================
# 00_setup.R - paths and functions shared by all R scripts.
# Every script starts with  source("code/R/00_setup.R")  and must be run from the
# repository root (the folder that contains data/, splits/ and code/).
# Only base R is used (R >= 4.1); no extra packages are needed.
# =============================================================================
if (!dir.exists("data") || !dir.exists("code")) stop("Run the R scripts from the repository root folder.")
DATA_DIR  <- "data"
SPLIT_DIR <- "splits"
OUT_DIR   <- file.path("code", "R", "outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path("data", "load_dataset.R"))
source(file.path("code", "R", "functions", "01_data_kernels.R"))
source(file.path("code", "R", "functions", "02_reml_metrics.R"))

GE_DATASETS <- c("EYT_1", "EYT_2", "EYT_3", "Indica", "Japonica",
                 "Wheat_1", "Wheat_2", "Wheat_3", "Wheat_5", "Wheat_6")

# Command-line arguments:  Rscript code/R/part1/xx.R --datasets Wheat_6 Japonica --partitions 1 2
arg_values <- function(flag, default) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (is.na(i)) return(default)
  j <- i + 1; v <- character(0)
  while (j <= length(a) && !startsWith(a[j], "--")) { v <- c(v, a[j]); j <- j + 1 }
  v
}
write_out <- function(df, file) {
  write.csv(df, file.path(OUT_DIR, file), row.names = FALSE)
  cat("WROTE", file.path(OUT_DIR, file), "\n")
}
log_msg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }
