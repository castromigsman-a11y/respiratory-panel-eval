"""A small, honestly-evaluated classifier for the calling decision.

The simple rule "reads >= threshold" ignores everything except depth. A
platform described as AI-guided would learn to combine features. Here a
logistic model uses read depth and genome breadth to decide presence, and is
assessed by stratified cross-validation so the reported AUC is not the training
score. It is deliberately a modest model: the goal is a fair comparison against
the threshold rule, not a leaderboard number.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


FEATURES = ["reads", "log_reads", "breadth"]


def _design(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    X = df.copy()
    X["log_reads"] = np.log10(X["reads"] + 1)
    return X[FEATURES], X["present"].to_numpy().astype(int)


def cross_validated_scores(df: pd.DataFrame, seed: int = 43, folds: int = 5) -> dict:
    X, y = _design(df)
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=1000, class_weight="balanced"))
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    auc = float(roc_auc_score(y, proba))
    frac_pos, mean_pred = calibration_curve(y, proba, n_bins=8, strategy="quantile")
    # fit once on all data to expose coefficients for the report
    model.fit(X, y)
    coefs = dict(zip(FEATURES, model.named_steps["logisticregression"].coef_[0]))
    return dict(
        proba=proba,
        y=y,
        auc=auc,
        calibration=(mean_pred, frac_pos),
        coefficients={k: float(v) for k, v in coefs.items()},
    )


def _call_metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    tp = int(np.sum((y == 1) & (pred == 1)))
    fp = int(np.sum((y == 0) & (pred == 1)))
    tn = int(np.sum((y == 0) & (pred == 0)))
    fn = int(np.sum((y == 1) & (pred == 0)))
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    return dict(sensitivity=sens, specificity=spec, ppv=ppv)


def compare_to_threshold(df: pd.DataFrame, cv: dict) -> pd.DataFrame:
    """Compare the depth-only threshold rule with the learned rule.

    Beyond AUC, the learned probabilities are thresholded to match the depth
    rule's sensitivity, so the specificity and PPV columns are a like-for-like
    comparison of calling quality at the same detection rate. The gain comes
    from rejecting high-depth, low-breadth pileups that the depth rule accepts.
    """
    y = df["present"].to_numpy().astype(int)
    depth_calls = df["called"].to_numpy().astype(int)
    depth_auc = float(roc_auc_score(y, df["reads"].to_numpy()))
    depth_m = _call_metrics(y, depth_calls)

    proba = cv["proba"]
    target_sens = depth_m["sensitivity"]
    # pick the probability cut that reproduces the depth rule's sensitivity
    order = np.argsort(-proba)
    cut = 0.5
    for idx in order:
        pred = (proba >= proba[idx]).astype(int)
        if _call_metrics(y, pred)["sensitivity"] >= target_sens:
            cut = proba[idx]
            break
    learned_calls = (proba >= cut).astype(int)
    learned_m = _call_metrics(y, learned_calls)

    return pd.DataFrame([
        dict(rule="reads >= threshold (depth only)", auc=depth_auc, **depth_m),
        dict(rule="logistic (depth + breadth), 5-fold CV", auc=cv["auc"], **learned_m),
    ])
