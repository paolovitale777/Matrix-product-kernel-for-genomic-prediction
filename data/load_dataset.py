"""load_dataset.py - read one dataset of the compact data package in Python.

    from data.load_dataset import load_dataset          # from the repository root
    d = load_dataset("Wheat_1")                          # or "EYT_1", "Indica", "Wheat_599", ...
    d["pheno"].head(); d["G"].shape

Returns a dict with the same content as load_dataset() in R (data/load_dataset.R):
    pheno    pandas DataFrame of phenotypes (long format: Env, line, traits; Cuevas data: wide,
             one column per environment)
    G        genomic relationship matrix (numpy, lines x lines), with G_ids
    K, A     Gaussian kernel exp(-d / median d) and pedigree matrix (Cuevas data only)
    markers  marker matrix (lines x markers) with marker_ids, when available
    EC       environmental covariates (rice datasets)
Record order and line order are those of the original files. Requires: numpy, pandas, rdata.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import rdata

DATA_DIR = Path(__file__).resolve().parent
DATASETS_GE = ["EYT_1", "EYT_2", "EYT_3", "Indica", "Japonica",
               "Wheat_1", "Wheat_2", "Wheat_3", "Wheat_5", "Wheat_6"]
DATASETS_CUEVAS = ["Wheat_599", "EYT_16_17", "EYT_22_23", "EYT_23_24"]


def _path(name: str, data_dir: Path) -> Path:
    for sub in ("GxE", "Cuevas2025", "."):
        p = data_dir / sub / f"{name}.rds"
        if p.exists():
            return p
    raise FileNotFoundError(f"dataset {name} not found under {data_dir}/GxE or {data_dir}/Cuevas2025")


def _unlower(L) -> tuple[np.ndarray, list[str] | None]:
    n = int(np.asarray(L["n"]).ravel()[0])
    X = np.zeros((n, n))
    rows, cols = np.tril_indices(n)
    order = np.lexsort((rows, cols))            # R stores the lower triangle column by column
    X[rows[order], cols[order]] = np.asarray(L["values"], float)
    X = X + np.tril(X, -1).T
    dn = L.get("dimnames") if hasattr(L, "get") else None
    ids = [str(v) for v in dn[0]] if dn is not None and dn[0] is not None else None
    return X, ids


def _matrix(x) -> tuple[np.ndarray, list[str] | None]:
    """xarray DataArray or ndarray from rdata -> (array, row names)."""
    if hasattr(x, "values") and hasattr(x, "dims"):
        ids = [str(v) for v in x.coords[x.dims[0]].values] if x.dims[0] in x.coords else None
        return np.asarray(x.values, float), ids
    return np.asarray(x, float), None


def load_dataset(name: str, data_dir: str | Path = DATA_DIR) -> dict:
    x = rdata.read_rds(str(_path(name, Path(data_dir))))
    out = {"name": name, "crop": str(np.asarray(x["crop"]).ravel()[0]), "pheno": x["pheno"]}
    if "markers_int" in x:                      # wheat: integer codes + imputed column means
        M, ids = _matrix(x["markers_int"])
        fill = np.asarray(x["markers_fill"], float)
        r, c = np.where(np.isnan(M))
        M[r, c] = fill[c]
        out["markers"], out["marker_ids"] = M, ids
    elif "markers" in x:
        out["markers"], out["marker_ids"] = _matrix(x["markers"])
    if "G_lower" in x:
        out["G"], out["G_ids"] = _unlower(x["G_lower"])
        if "geno_line_order" in x:
            out["G_ids"] = [str(v) for v in x["geno_line_order"]]
    elif "markers" in out and out["crop"] == "wheat":
        M = out["markers"]
        out["G"], out["G_ids"] = M @ M.T / M.shape[1], out["marker_ids"]
    if "A_lower" in x:
        out["A"], _ = _unlower(x["A_lower"])
        g = np.diag(out["G"])
        d = g[:, None] + g[None, :] - 2 * out["G"]
        out["K"] = np.exp(-d / np.median(d))
        out["env_names"] = [str(v) for v in x["env_names"]]
    if "EC" in x:
        out["EC"] = x["EC"]
    return out


if __name__ == "__main__":
    for nm in DATASETS_GE + DATASETS_CUEVAS:
        try:
            d = load_dataset(nm)
        except FileNotFoundError as e:
            print(e)
            continue
        print(nm, d["pheno"].shape, "G", None if "G" not in d else d["G"].shape,
              "markers", None if "markers" not in d else d["markers"].shape)
