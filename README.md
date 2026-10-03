# Record order and matrix-product kernels in genomic prediction — data and code

Data, code (Python and R) and results for two companion papers submitted to *Frontiers in Genetics*:

* **Part I** — *Accuracy gains of triangular matrix-product kernels in genomic prediction reflect how the data are sorted.*
* **Part II** — *What the file order knew: line identifiers, environment-specific signal and non-linear kernels in genomic prediction.*

Authors: José Crossa, Jaime Cuevas, Osval A. Montesinos-López, Abelardo Montesinos-López, Paulino Pérez-Rodríguez, Paolo Vitale, Rodomiro Ortiz.

## Contents (three parts, also released as three ZIP files on Zenodo)

| Folder | Content | Size |
|---|---|---|
| `data/GxE/` | the ten genotype × environment datasets (EYT_1–3, Indica, Japonica, Wheat_1–3, 5, 6) | ~17 MB |
| `data/Cuevas2025/` | three datasets of Cuevas et al. (2025): Wheat 599, EYT_16_17, EYT_22_23 (EYT_23_24: see below) | ~4 MB |
| `code/` + `splits/` | Python code (published analyses), R code (same analyses), exported splits, outputs | ~4 MB |

* `data/` — one compact `.rds` file per dataset (phenotypes, markers or genomic relationships, pedigree, environmental covariates) plus a CSV of its phenotypes; loaders for R (`load_dataset.R`) and Python (`load_dataset.py`); `rebuild_original_files.R`, which recreates the original `.RData` files exactly; see `data/README_data.md`.
* `splits/` — every random element of the analyses (50:50 partitions, random record orders, within-environment shuffles, permutations of the matrix-level test, 80:20 splits of the Cuevas data), exported once from the Python code so that R and Python analyse identical splits.
* `code/python/part1`, `code/python/part2` — the Python code used for the papers, with the outputs reported in them (`outputs/`).
* `code/R/` — R translation of every analysis reported in the papers (one script per table), a comparison script and `run_all.R`; see `code/R/README_R.md`.
* `docs/` — guide for readers who want to check where each result comes from.

## Quick start

**R** (R ≥ 4.1, base R only), from this folder:

```r
source("data/load_dataset.R")
d <- load_dataset("Wheat_1")      # d$pheno, d$G, d$markers
```
```
Rscript code/R/run_all.R --quick        # every analysis on a small dataset (~30 min), then comparison with Python
Rscript code/R/run_all.R                # all analyses (several hours)
```

**Python** (3.11–3.12; `pip install -r code/python/part1/requirements.txt`):

```
Rscript data/rebuild_original_files.R   # writes the original .RData files used by the Python code
python code/python/part1/code/selftest_synthetic.py   # 29 checks, no data needed, seconds
```
and then `code/python/part1/run_all.bat` (or `.sh`) and `code/python/part2/run_all.bat`.

## Where each result comes from

| Paper | Table | R script (`code/R/`) | Python script | Python output |
|---|---|---|---|---|
| I | Table 1 | `part1/P1_T1_datasets_rank.R` | `analysis_pipeline.py qc`, `matrix_rank.py` | `dataset_qc_manifest.csv` |
| I | Table 2 | `part1/P1_T2_permutation_test.R` | `analysis_pipeline.py order` | `order_diagnostics_*.csv` |
| I | Section 3.2 | `part1/P1_S32_psd_diagnostics.R` | `analysis_pipeline.py psd` | `psd_diagnostics.csv` |
| I | Table 3 (B, F shuffled, R6), 3.4 (R5) | `part1/P1_T3_identifier_diagnostics.R` | `id_order_test.py` | `id_order_metrics.csv` |
| I | Table 3 (F file order, random orders) | `part1/P1_T3_record_order_accuracy.R` | `order_accuracy_sensitivity.py` | `order_accuracy_sensitivity*.csv` |
| I | Table 4, Table S1 | `part1/P1_T4_S1_benchmark_5050.R` | `benchmark_v4.py --scenario ulte` | `v4_ulte_metrics.csv` |
| I | Table 5 | `part1/P1_T5_prospective_year.R` | `prospective_gy.py` | `prospective_gy_metrics*.csv` |
| I | Table 6 | `part1/P1_T6_cuevas_triangles.R` | `cuevas_check_files.py` | `cuevas_files_*.csv` |
| I | Table S2 | `part1/P1_S2_connectivity_models.R` | `recover_gain.py` | `recover_gain_metrics.csv` |
| II | Table 3 | `part2/P2_T3_cross_covariance.R` | `part2_cross_ka_v2.py` | `part2_cross_v2_*.csv` |
| II | Table 4 | `part2/P2_T4_identifier_neighbours.R` | `part2_neighbour_env*.py`, `part2_id_relatedness.py` | `part2_neighbour_*.csv`, `part2_id_relatedness.csv` |
| II | Tables 5–6 (and S2) | `part2/P2_T5_T6_identifier_split.R` | `part2_id_split.py` | `part2_id_split_*.csv` |
| II | Table 7, Table S1 (and S2) | `part2/P2_T7_S1_breeding_scenarios.R` | `part2_scenarios.py` | `part2_scenarios_*.csv` |

## Data sources and licences

* The ten G × E datasets were published with Montesinos-López et al. (2026) at https://github.com/osval78/New_GE_Framework.
* The datasets of Cuevas et al. (2025) are deposited in the CIMMYT Research Data & Software Repository Network (https://hdl.handle.net/11529/10549160). Three of them (Wheat 599, EYT_16_17, EYT_22_23) are redistributed here with the agreement of the dataset authors (J. Cuevas and co-authors). The fourth, **EYT_23_24**, is not included: download it from the CIMMYT repository (access on request) and add it with `Rscript data/add_cuevas_file.R path/to/EYT_23_24.Rdata`. Please cite Cuevas et al. (2025) and the CIMMYT repository when using these data.
* Code: MIT License (`LICENSE`). Data: please cite the original publications and this repository (Zenodo DOI in `CITATION.cff`).

## How to cite

Please cite the two papers (references will be added on publication) and this repository through its Zenodo DOI.
