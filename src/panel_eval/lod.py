"""Limit of detection (LoD) by probit regression.

This follows the logic clinical labs use under CLSI EP17: run a dilution series,
record the detection rate at each concentration, fit detection probability
against log10 input, and read off the concentration giving 95% detection.

We fit a binomial GLM with a probit link:

    P(detect | x) = Phi(b0 + b1 * x),   x = log10(copies/mL)

LoD95 solves Phi(b0 + b1 * x) = 0.95, i.e. x = (z_0.95 - b0) / b1, then back to
the concentration scale. The confidence interval on x uses the delta method on
the fitted coefficient covariance, which is standard for this design.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning


def hit_rates(df: pd.DataFrame) -> pd.DataFrame:
    """Detection rate per (pathogen, input level) among spiked samples."""
    spiked = df[df["present"] == 1].copy()
    grouped = (
        spiked.groupby(["pathogen", "input_log10"])
        .agg(n=("called", "size"), hits=("called", "sum"))
        .reset_index()
    )
    grouped["hit_rate"] = grouped["hits"] / grouped["n"]
    return grouped


def _fit_probit(x: np.ndarray, hits: np.ndarray, n: np.ndarray):
    X = sm.add_constant(x)
    endog = np.column_stack([hits, n - hits])
    model = sm.GLM(endog, X, family=sm.families.Binomial(link=sm.families.links.Probit()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PerfectSeparationWarning)
        return model.fit()


def lod_for_pathogen(sub: pd.DataFrame, target: float = 0.95) -> dict:
    """Fit the probit model for one pathogen and return LoD at `target` detection."""
    sub = sub.sort_values("input_log10")
    x = sub["input_log10"].to_numpy(dtype=float)
    hits = sub["hits"].to_numpy(dtype=float)
    n = sub["n"].to_numpy(dtype=float)

    # Need variation in the response to fit a slope.
    if len(x) < 3 or sub["hit_rate"].nunique() < 2:
        return dict(lod_log10=np.nan, lod_copies=np.nan,
                    ci_low=np.nan, ci_high=np.nan, slope=np.nan, converged=False)

    fit = _fit_probit(x, hits, n)
    b0, b1 = fit.params
    z = norm.ppf(target)

    if b1 <= 0:
        return dict(lod_log10=np.nan, lod_copies=np.nan,
                    ci_low=np.nan, ci_high=np.nan, slope=float(b1), converged=False)

    x_lod = (z - b0) / b1

    # delta method: g(b) = (z - b0) / b1
    cov = fit.cov_params()
    dg_db0 = -1.0 / b1
    dg_db1 = -(z - b0) / (b1 ** 2)
    grad = np.array([dg_db0, dg_db1])
    var_x = float(grad @ np.asarray(cov) @ grad)
    se_x = np.sqrt(max(var_x, 0.0))

    return dict(
        lod_log10=float(x_lod),
        lod_copies=float(10 ** x_lod),
        ci_low=float(10 ** (x_lod - 1.96 * se_x)),
        ci_high=float(10 ** (x_lod + 1.96 * se_x)),
        slope=float(b1),
        converged=True,
    )


def lod_table(df: pd.DataFrame, target: float = 0.95) -> pd.DataFrame:
    rates = hit_rates(df)
    records = []
    for name, sub in rates.groupby("pathogen"):
        res = lod_for_pathogen(sub, target=target)
        res["pathogen"] = name
        records.append(res)
    cols = ["pathogen", "lod_copies", "ci_low", "ci_high", "lod_log10", "slope", "converged"]
    return pd.DataFrame(records)[cols].sort_values("lod_copies").reset_index(drop=True)


def probit_curve(sub: pd.DataFrame, grid: np.ndarray | None = None):
    """Return (x grid, fitted P) for plotting a single pathogen's dose response."""
    sub = sub.sort_values("input_log10")
    x = sub["input_log10"].to_numpy(dtype=float)
    hits = sub["hits"].to_numpy(dtype=float)
    n = sub["n"].to_numpy(dtype=float)
    if sub["hit_rate"].nunique() < 2:
        return None
    fit = _fit_probit(x, hits, n)
    if grid is None:
        grid = np.linspace(x.min() - 0.5, x.max() + 0.5, 200)
    Xg = sm.add_constant(grid)
    return grid, fit.predict(Xg)
