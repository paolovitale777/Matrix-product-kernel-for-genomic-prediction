# Data

Grain-yield datasets used in Part I and Part II (13 of the 14 are included; see EYT_23_24 below). Each dataset is one compact `.rds` file
(R serialisation, xz-compressed) and a CSV file with its phenotypes. The files can be read in R
(`readRDS`, or `load_dataset()` in `load_dataset.R`) and in Python (`load_dataset()` in
`load_dataset.py`, which uses the `rdata` package).

**Record order and line order are exactly those of the original files.** This matters: Part I
shows that the published triangular kernels depend on the order of the records.

## Files

### `GxE/` — ten genotype × environment datasets (Montesinos-López et al., 2026)

| Dataset | Crop | Env. | Lines | Records | Genomic information | Other |
|---|---|---|---|---|---|---|
| EYT_1 | wheat | 4 | 766 | 3,064 | genomic relationship matrix G (lower triangle) | traits DTHD, DTMT, GY, Height |
| EYT_2 | wheat | 5 | 775 | 3,875 | G | as EYT_1 |
| EYT_3 | wheat | 5 | 964 | 4,820 | G | as EYT_1 |
| Indica | rice | 3 | 327 | 981 | 16,383 SNP markers | 54 environmental covariates per year (EC) |
| Japonica | rice | 5 | 320 | 1,151 | 16,383 SNP markers | EC |
| Wheat_1, 2, 3, 5, 6 | wheat | 2 | 1,275–1,403 | 1,400–1,933 | 5,741 markers | — |

Original source: https://github.com/osval78/New_GE_Framework (Wheat_4 is not available there and
was not analysed).

### `Cuevas2025/` — data of Cuevas et al. (2025)

| Dataset | Lines | Env. | Stored | Recomputed when loaded |
|---|---|---|---|---|
| Wheat_599 | 599 | 4 | pedigree A, genomic G | Gaussian kernel K |
| EYT_16_17 | 1,028 | 6 | pedigree A, 3,714 SNPs | G = M M'/m, K |
| EYT_22_23 | 1,004 | 6 | pedigree A, 3,672 SNPs | G = M M'/m, K |
| EYT_23_24 | 898 | 5 | not included — see below | |

Phenotypes are stored as a wide table (one column per environment). Original source: CIMMYT
Research Data & Software Repository Network, https://hdl.handle.net/11529/10549160
("Replication Data for: Enhancing Genomic Prediction Accuracy in Precision Breeding through
Combined Pedigree and Genomic Matrices"); redistributed with the agreement of the dataset authors.

**EYT_23_24 is not included.** Download `EYT_23_24.Rdata` from the CIMMYT repository (access on
request) and add it with

    Rscript data/add_cuevas_file.R path/to/EYT_23_24.Rdata

which converts it to `data/Cuevas2025/EYT_23_24.rds`; the R scripts then include it automatically
(for the Python code, run `rebuild_original_files.R` and `code/python/export_splits.py` again).

## How the files were made smaller without changing the data

Matrices that can be recomputed exactly are not stored; `load_dataset()` recomputes them:

* Wheat_1–6: G = M M'/m from the stored markers (identical to the supplied G).
* Wheat_1–6 markers: stored as integer codes 0/1/2 with the imputed value (column mean) of each
  marker kept separately; the original marker matrix is restored exactly.
* EYT_16_17, EYT_22_23: G = M M'/m from the stored SNPs (identical to the supplied G).
* Cuevas data: K = exp(−d / median d), d_jk = G_jj + G_kk − 2 G_jk (equal to the supplied K up to
  floating-point rounding, relative difference < 4 × 10⁻¹⁶).
* Symmetric matrices (G of EYT_1–3, A, G of Wheat_599) are stored as their lower triangle.

`add_cuevas_file.R` adds a dataset of Cuevas et al. (2025) downloaded from CIMMYT.
`rebuild_original_files.R` writes the original `.RData` files (same objects, names, values and
order; verified object by object) for the Python code in `code/python/`.

## Fields of the `.rds` files

`name`, `crop`, `pheno` (data frame), and, depending on the dataset: `G_lower` (lower triangle of
G), `geno_line_order`, `markers` or `markers_int` + `markers_fill`, `EC`, `A_lower`, `env_names`,
`source`, `notes`. `load_dataset()` returns the ready-to-use objects `pheno`, `G`, `K`, `A`,
`markers`, `EC`.
