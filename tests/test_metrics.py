"""Sanity tests for the evaluation core.

Run with: pytest -q
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from panel_eval import agreement, lod, metrics  # noqa: E402
from panel_eval.simulate_counts import simulate  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "config" / "panel.yaml"


def test_confusion_counts():
    y = np.array([1, 1, 0, 0])
    p = np.array([1, 0, 0, 1])
    c = metrics.confusion(y, p)
    assert (c["tp"], c["fn"], c["tn"], c["fp"]) == (1, 1, 1, 1)


def test_perfect_prediction_gives_unit_metrics():
    y = np.array([1, 1, 0, 0, 1, 0])
    tbl = metrics.dta_metrics(y, y)
    for m in ("sensitivity", "specificity", "ppv", "npv", "accuracy"):
        est = tbl.loc[tbl["metric"] == m, "estimate"].iloc[0]
        assert abs(est - 1.0) < 1e-9


def test_sensitivity_matches_manual():
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    p = np.array([1, 1, 1, 0, 0, 0, 0, 0])  # 3/4 caught, no false positives
    tbl = metrics.dta_metrics(y, p)
    sens = tbl.loc[tbl["metric"] == "sensitivity", "estimate"].iloc[0]
    spec = tbl.loc[tbl["metric"] == "specificity", "estimate"].iloc[0]
    assert abs(sens - 0.75) < 1e-9
    assert abs(spec - 1.0) < 1e-9


def test_kappa_bounds():
    y = np.array([1, 0, 1, 0, 1, 0])
    assert abs(metrics.kappa(y, y) - 1.0) < 1e-9
    assert metrics.kappa(y, 1 - y) <= 0.0


def test_lod_ordering_is_monotone_in_c50():
    """A pathogen configured with a lower c50 must not have a higher LoD."""
    df = simulate(CONFIG, seed=43)
    tbl = lod.lod_table(df)
    tbl = tbl[tbl["converged"]]
    # SARS_CoV_2 has the lowest configured c50, so the smallest LoD.
    assert tbl.iloc[0]["pathogen"] == "SARS_CoV_2"
    assert tbl["lod_copies"].is_monotonic_increasing


def test_lod_confidence_interval_is_finite():
    df = simulate(CONFIG, seed=43)
    tbl = lod.lod_table(df)
    conv = tbl[tbl["converged"]]
    assert np.isfinite(conv[["lod_copies", "ci_low", "ci_high"]].to_numpy()).all()
    assert (conv["ci_low"] < conv["lod_copies"]).all()
    assert (conv["lod_copies"] < conv["ci_high"]).all()


def test_read_depth_tracks_input():
    df = simulate(CONFIG, seed=43)
    corr = agreement.correlation(df)
    assert corr["spearman_r"] > 0.5


def test_simulation_is_deterministic():
    a = simulate(CONFIG, seed=43)
    b = simulate(CONFIG, seed=43)
    pd.testing.assert_frame_equal(a, b)
