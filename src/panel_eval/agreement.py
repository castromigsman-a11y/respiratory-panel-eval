"""Quantitative agreement between the pipeline's signal and the known input.

A qualitative panel is validated on detection, but a sequencing platform also
carries quantitative information: more input template gives more reads. This
module calibrates read count to input load, then characterises agreement the
way a method-comparison study would, with correlation and a Bland-Altman plot
of estimated versus known load on the log10 scale.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def _positives(df: pd.DataFrame) -> pd.DataFrame:
    """Spiked samples that were actually called, i.e. have usable signal."""
    pos = df[(df["present"] == 1) & (df["called"] == 1)].copy()
    pos["log_reads"] = np.log10(pos["reads"] + 1)
    return pos


def calibrate(df: pd.DataFrame) -> dict:
    """Ordinary least squares of log10(reads+1) on known log10 input.

    The inverse mapping turns an observed read count into an estimated load,
    which is what the agreement analysis compares against truth.
    """
    pos = _positives(df)
    x = pos["input_log10"].to_numpy()
    y = pos["log_reads"].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    return dict(slope=float(slope), intercept=float(intercept), n=int(len(pos)))


def estimate_load(df: pd.DataFrame, cal: dict) -> pd.DataFrame:
    pos = _positives(df)
    # invert y = slope * x + intercept  ->  x_hat = (y - intercept) / slope
    pos["est_log10"] = (pos["log_reads"] - cal["intercept"]) / cal["slope"]
    return pos


def correlation(df: pd.DataFrame) -> dict:
    pos = _positives(df)
    r_p, p_p = pearsonr(pos["input_log10"], pos["log_reads"])
    r_s, p_s = spearmanr(pos["input_log10"], pos["log_reads"])
    return dict(pearson_r=float(r_p), pearson_p=float(p_p),
                spearman_r=float(r_s), spearman_p=float(p_s))


def bland_altman(df: pd.DataFrame, cal: dict) -> dict:
    """Bias and limits of agreement between estimated and known log10 load."""
    pos = estimate_load(df, cal)
    known = pos["input_log10"].to_numpy()
    est = pos["est_log10"].to_numpy()
    mean_load = (known + est) / 2
    diff = est - known
    bias = float(np.mean(diff))
    sd = float(np.std(diff, ddof=1))
    return dict(
        mean_load=mean_load,
        diff=diff,
        bias=bias,
        sd=sd,
        loa_low=bias - 1.96 * sd,
        loa_high=bias + 1.96 * sd,
    )
