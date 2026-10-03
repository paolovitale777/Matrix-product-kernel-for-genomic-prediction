# R code

R translation of every analysis reported in Part I and Part II. It uses only base R (R ≥ 4.1)
and reads the data directly from `data/` (no conversion needed) and the splits from `splits/`.
All scripts are run from the repository root, e.g.

```
Rscript code/R/part1/P1_T4_S1_benchmark_5050.R --datasets Wheat_6 --partitions 1 2
Rscript code/R/compare_with_python.R
```

## Files

```
code/R/
  00_setup.R                       paths, sources the functions, command-line options
  functions/01_data_kernels.R      data access, kernels, direction masks, triangular and proposed components
  functions/02_reml_metrics.R      training-anchored PSD projection, REML fitting, metrics, splits
  part1/P1_T1_datasets_rank.R      Part I Table 1
  part1/P1_T2_permutation_test.R   Part I Table 2
  part1/P1_S32_psd_diagnostics.R   Part I Section 3.2
  part1/P1_T3_identifier_diagnostics.R   Part I Table 3 (B, F lines shuffled, R6) and Section 3.4 (R5)
  part1/P1_T3_record_order_accuracy.R    Part I Table 3 (F file order, F random orders)
  part1/P1_T4_S1_benchmark_5050.R  Part I Table 4 and Supplementary Table S1
  part1/P1_T5_prospective_year.R   Part I Table 5
  part1/P1_T6_cuevas_triangles.R   Part I Table 6 (also Part II Section 3.1)
  part1/P1_S2_connectivity_models.R      Part I Supplementary Table S2
  part2/P2_T3_cross_covariance.R   Part II Table 3
  part2/P2_T4_identifier_neighbours.R    Part II Table 4
  part2/P2_T5_T6_identifier_split.R      Part II Tables 5 and 6 (PM10 for Table S2)
  part2/P2_T7_S1_breeding_scenarios.R    Part II Table 7 and Supplementary Table S1 (PM10 for Table S2)
  compare_with_python.R            R versus published Python values, cell by cell
  run_all.R                        runs everything (option --quick: one small dataset per script)
  outputs/                         results written by the scripts
```

## The functions and the equations they implement

| R function (file) | Python function | What it does (paper reference) |
|---|---|---|
| `genomic_matrix` (01) | `analysis_pipeline.genomic_matrix` | supplied G, or centred markers / Σ marker variances; symmetrised; unit mean diagonal |
| `environmental_matrix` (01) | `environmental_matrix` | standardised covariates, E = W W′/q (Part I Eq. 12), unit mean diagonal; identity without covariates |
| `record_kernels` (01) | `record_kernels` | K_G = Z_G G Z_G′, K_E = Z_E E Z_E′ (Part I Eq. 7) |
| `ordered_env_levels`, `environment_time` (01) | same | environments ordered by label (and calendar year), never by file position |
| `direction_mask` (01) | `direction_mask` | B: b_ij = 1 if the environment of i precedes that of j (Eq. 13); calendar year or Bed2IR ≺ Bed5IR |
| `triangles` (01) | `triangles` / `build_components` | C0 = D + U + U′, P0 = D + L + L′ (Eqs. 5–6, 9) |
| `build_components` (01) | `build_components` | E, G, H = K_E∘K_G (Eq. 1), C0, P0, C_SP, P_SP, T_SP or S_GE (Eqs. 14–15) |
| `identifier_kernel` (01) | `id_order_test`, `part2_id_split` | K_id(j,k) = exp(−\|rank_j − rank_k\|/h), h = 5% of the lines |
| `gaussian_kernel` (01) | `part2_id_split.gaussian_kernel` | exp(−d/median d), d_ij = K_ii + K_jj − 2K_ij (Part II Section 2.4) |
| `psd_blocks` (02) | `benchmark_ulte.psd_blocks` | training-anchored projection and extension (Part I Section 2.6.2) |
| `reml_core` (02) | REML objective of all scripts | REML criterion and its exact gradient; L-BFGS-B on log variances |
| `reml_fit_predict`, `reml_fixed`, `reml_mixed`, `reml_robust` (02) | same names | the four REML variants of the Python code (same bounds, starts, tolerances) |
| `evaluate`, `basic_metrics` (02) | same | COR, NRMSE, PM10 within environments (averaged) and pooled (Part I Section 2.8) |
| `paired_differences` (02) | same | same-partition differences with the reference model |

## Agreement with the Python results

R and Python analyse the same splits (exported in `splits/`). Kernels, projections and metrics
are identical; the REML optimiser differs (R: exact gradient; Python: finite differences), so the
variance components and predictions agree to the optimisation tolerance. Typical differences
in predictive correlation are below 0.001 (see `docs/` and `compare_with_python.R`).

Two elements use R's random generator and therefore differ only by Monte Carlo noise: the
randomised eigenvalue routine in `P1_S32_psd_diagnostics.R` (eigenvalues identical to rounding)
and nothing else; all other random elements are read from `splits/`.

## Run time

R's default (reference) BLAS is about 6× slower than NumPy's. With it, `run_all.R --quick` takes
about 30 minutes and the full run several hours; the EYT datasets (3,000–4,800 records) dominate.
Linking R to OpenBLAS or MKL reduces the time several-fold without changing the results.
