"""Diagnostic test accuracy (DTA) metrics.

Given binary ground truth and binary calls, compute the quantities a clinical
lab reports when validating a qualitative assay: sensitivity, specificity, PPV,
NPV, F1 and accuracy, each with a Wilson score interval. Also ROC/AUC using a
continuous score (read count) and Cohen's kappa for categorical agreement.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score, roc_auc_score, roc_curve
from statsmodels.stats.proportion import proportion_confint


def confusion(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, int]:
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    return dict(tp=tp, fp=fp, tn=tn, fn=fn)


def _ratio_ci(numer: int, denom: int, alpha: float = 0.05) -> tuple[float, float, float]:
    if denom == 0:
        return (np.nan, np.nan, np.nan)
    point = numer / denom
    low, high = proportion_confint(numer, denom, alpha=alpha, method="wilson")
    return point, low, high


def dta_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> pd.DataFrame:
    """Return a tidy table of DTA metrics with Wilson 95% intervals."""
    c = confusion(y_true, y_pred)
    tp, fp, tn, fn = c["tp"], c["fp"], c["tn"], c["fn"]

    definitions = {
        "sensitivity": (tp, tp + fn),
        "specificity": (tn, tn + fp),
        "ppv": (tp, tp + fp),
        "npv": (tn, tn + fn),
        "accuracy": (tp + tn, tp + fp + tn + fn),
    }

    records = []
    for name, (numer, denom) in definitions.items():
        point, low, high = _ratio_ci(numer, denom)
        records.append(dict(metric=name, estimate=point, ci_low=low, ci_high=high, n=denom))

    precision = tp / (tp + fp) if (tp + fp) else np.nan
    recall = tp / (tp + fn) if (tp + fn) else np.nan
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and (precision + recall) else np.nan)
    records.append(dict(metric="f1", estimate=f1, ci_low=np.nan, ci_high=np.nan,
                        n=tp + fp + fn))

    return pd.DataFrame.from_records(records)


def kappa(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(cohen_kappa_score(y_true, y_pred))


def threshold_sweep(scores: np.ndarray, y_true: np.ndarray,
                    thresholds: np.ndarray | None = None) -> pd.DataFrame:
    """Sensitivity and specificity as the read-count calling threshold moves."""
    scores = np.asarray(scores)
    y_true = np.asarray(y_true).astype(int)
    if thresholds is None:
        top = int(np.percentile(scores, 99)) + 2
        thresholds = np.arange(1, max(top, 3) + 1)
    out = []
    for t in thresholds:
        pred = (scores >= t).astype(int)
        c = confusion(y_true, pred)
        sens = c["tp"] / (c["tp"] + c["fn"]) if (c["tp"] + c["fn"]) else np.nan
        spec = c["tn"] / (c["tn"] + c["fp"]) if (c["tn"] + c["fp"]) else np.nan
        out.append(dict(threshold=int(t), sensitivity=sens, specificity=spec,
                        youden=(sens + spec - 1)))
    return pd.DataFrame(out)


def roc_per_pathogen(df: pd.DataFrame, score_col: str = "reads") -> pd.DataFrame:
    """AUC per pathogen, using the continuous score to separate present/absent."""
    records = []
    for name, sub in df.groupby("pathogen"):
        y = sub["present"].to_numpy()
        s = sub[score_col].to_numpy()
        if len(np.unique(y)) < 2:
            auc = np.nan
        else:
            auc = float(roc_auc_score(y, s))
        records.append(dict(pathogen=name, auc=auc, n=len(sub),
                            positives=int(y.sum())))
    return pd.DataFrame(records).sort_values("auc", ascending=False).reset_index(drop=True)


def roc_points(y_true: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fpr, tpr, _ = roc_curve(y_true, scores)
    return fpr, tpr
