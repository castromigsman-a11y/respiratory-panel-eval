"""Run the full diagnostic evaluation on a count matrix.

Reads a count matrix (the committed example, or the output of the real
pipeline, which shares the same schema), computes every analysis, writes
figures to results/figures, tables to results/tables, and assembles
results/REPORT.md. Run with no arguments to use the committed example data.

    python scripts/run_evaluation.py
    python scripts/run_evaluation.py --counts path/to/counts.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from panel_eval import agreement, classifier, lod, metrics  # noqa: E402
from panel_eval.simulate_counts import write_example  # noqa: E402

FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
PALETTE = ["#2b6cb0", "#c05621", "#2f855a", "#6b46c1", "#b7791f", "#c53030"]


def _style():
    plt.rcParams.update({
        "figure.dpi": 130,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
    })


CLINICAL_LOAD_LOG10 = 4.0  # >= 1e4 copies/mL, representative of clinical samples


def dta_at_clinical_load(df: pd.DataFrame) -> pd.DataFrame:
    """DTA using only clinically relevant positives, plus all true negatives.

    The full dilution series deliberately probes far below any clinical load, so
    pooled sensitivity understates real performance. Restricting the positive
    side to >= 1e4 copies/mL gives the number a lab would report, while the
    negative side keeps every pathogen-absent observation.
    """
    pos = df[(df["present"] == 1) & (df["input_log10"] >= CLINICAL_LOAD_LOG10)]
    neg = df[df["present"] == 0]
    sub = pd.concat([pos, neg])
    return metrics.dta_metrics(sub["present"].to_numpy(), sub["called"].to_numpy())


def pooled_detection(df: pd.DataFrame) -> pd.DataFrame:
    """Detection rate by input level, pooled across all pathogens."""
    spiked = df[df["present"] == 1]
    tbl = (spiked.groupby("input_log10")
           .agg(n=("called", "size"), detected=("called", "sum"))
           .reset_index())
    tbl["detection_rate"] = tbl["detected"] / tbl["n"]
    tbl["copies_per_ml"] = (10 ** tbl["input_log10"]).round().astype(int)
    return tbl[["copies_per_ml", "input_log10", "n", "detected", "detection_rate"]]


def load_counts(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"count matrix not found at {path}, generating example data")
        write_example(ROOT / "config" / "panel.yaml", path)
    df = pd.read_csv(path)
    if "called" not in df.columns:
        # infer from a threshold read from the panel config if needed
        df["called"] = (df["reads"] >= 5).astype(int)
    return df


def fig_threshold_sweep(df: pd.DataFrame) -> str:
    sweep = metrics.threshold_sweep(df["reads"].to_numpy(), df["present"].to_numpy())
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(sweep["threshold"], sweep["sensitivity"], "-o", color=PALETTE[0], label="sensitivity")
    ax.plot(sweep["threshold"], sweep["specificity"], "-o", color=PALETTE[1], label="specificity")
    best = sweep.loc[sweep["youden"].idxmax(), "threshold"]
    ax.axvline(best, color="grey", ls="--", lw=1, label=f"max Youden (t={int(best)})")
    ax.set_xlabel("calling threshold (reads)")
    ax.set_ylabel("rate")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Detection trade-off vs calling threshold")
    out = FIG / "threshold_sweep.png"
    fig.tight_layout(); fig.savefig(out); plt.close(fig)
    sweep.to_csv(TAB / "threshold_sweep.csv", index=False)
    return out.name


def fig_roc(df: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(5.2, 5))
    auc_tbl = metrics.roc_per_pathogen(df)
    for i, name in enumerate(auc_tbl["pathogen"]):
        sub = df[df["pathogen"] == name]
        fpr, tpr = metrics.roc_points(sub["present"].to_numpy(), sub["reads"].to_numpy())
        auc = auc_tbl.loc[auc_tbl["pathogen"] == name, "auc"].iloc[0]
        ax.plot(fpr, tpr, color=PALETTE[i % len(PALETTE)], lw=1.6,
                label=f"{name} (AUC {auc:.3f})")
    ax.plot([0, 1], [0, 1], color="grey", ls=":", lw=1)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_title("Per-pathogen detection ROC")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    out = FIG / "roc_per_pathogen.png"
    fig.tight_layout(); fig.savefig(out); plt.close(fig)
    return out.name


def fig_lod(df: pd.DataFrame) -> tuple[str, pd.DataFrame]:
    rates = lod.hit_rates(df)
    tbl = lod.lod_table(df)
    order = tbl["pathogen"].tolist()
    ncol = 3
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.1 * nrow), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for i, name in enumerate(order):
        ax = axes[i]
        sub = rates[rates["pathogen"] == name]
        ax.scatter(sub["input_log10"], sub["hit_rate"], color=PALETTE[i % len(PALETTE)],
                   zorder=3, s=28)
        curve = lod.probit_curve(sub)
        if curve is not None:
            gx, gy = curve
            ax.plot(gx, gy, color=PALETTE[i % len(PALETTE)], lw=1.6)
        row = tbl[tbl["pathogen"] == name].iloc[0]
        if row["converged"]:
            ax.axhline(0.95, color="grey", ls=":", lw=1)
            ax.axvline(row["lod_log10"], color="grey", ls="--", lw=1)
            ax.set_title(f"{name}\nLoD95 = {row['lod_copies']:.0f} copies/mL", fontsize=9)
        else:
            ax.set_title(name, fontsize=9)
        ax.set_ylim(-0.03, 1.05)
        ax.set_xlabel("log10 input (copies/mL)")
        if i % ncol == 0:
            ax.set_ylabel("detection rate")
    for j in range(len(order), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Limit of detection by probit regression (CLSI EP17 style)", y=1.02)
    out = FIG / "lod_probit.png"
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    tbl.to_csv(TAB / "lod95.csv", index=False)
    return out.name, tbl


def fig_agreement(df: pd.DataFrame) -> tuple[str, dict, dict]:
    cal = agreement.calibrate(df)
    corr = agreement.correlation(df)
    ba = agreement.bland_altman(df, cal)
    pos = agreement.estimate_load(df, cal)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    ax1.scatter(pos["input_log10"], pos["log_reads"], s=16, alpha=0.5, color=PALETTE[0])
    xs = np.linspace(pos["input_log10"].min(), pos["input_log10"].max(), 100)
    ax1.plot(xs, cal["slope"] * xs + cal["intercept"], color=PALETTE[5], lw=1.6,
             label=f"fit (r={corr['pearson_r']:.2f})")
    ax1.set_xlabel("known log10 input (copies/mL)")
    ax1.set_ylabel("log10(reads + 1)")
    ax1.set_title("Read depth vs input load")
    ax1.legend(frameon=False, fontsize=8)

    ax2.scatter(ba["mean_load"], ba["diff"], s=16, alpha=0.5, color=PALETTE[2])
    ax2.axhline(ba["bias"], color=PALETTE[5], lw=1.4, label=f"bias {ba['bias']:.2f}")
    ax2.axhline(ba["loa_low"], color="grey", ls="--", lw=1)
    ax2.axhline(ba["loa_high"], color="grey", ls="--", lw=1,
                label=f"95% LoA [{ba['loa_low']:.2f}, {ba['loa_high']:.2f}]")
    ax2.set_xlabel("mean of estimated and known log10 load")
    ax2.set_ylabel("estimated - known (log10)")
    ax2.set_title("Bland-Altman: quantitative agreement")
    ax2.legend(frameon=False, fontsize=8)

    out = FIG / "agreement.png"
    fig.tight_layout(); fig.savefig(out); plt.close(fig)
    pd.DataFrame([{**corr, **{k: ba[k] for k in ("bias", "sd", "loa_low", "loa_high")},
                   "cal_slope": cal["slope"], "cal_intercept": cal["intercept"]}]
                 ).to_csv(TAB / "agreement.csv", index=False)
    return out.name, corr, ba


def fig_classifier(df: pd.DataFrame) -> tuple[str, pd.DataFrame, dict]:
    cv = classifier.cross_validated_scores(df)
    comp = classifier.compare_to_threshold(df, cv)
    mean_pred, frac_pos = cv["calibration"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.bar(comp["rule"], comp["auc"], color=[PALETTE[1], PALETTE[0]])
    ax1.set_ylim(0.5, 1.0)
    ax1.set_ylabel("AUC")
    ax1.set_title("Calling rule comparison")
    ax1.tick_params(axis="x", labelrotation=12, labelsize=7.5)
    for i, v in enumerate(comp["auc"]):
        ax1.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=8)

    ax2.plot([0, 1], [0, 1], color="grey", ls=":", lw=1)
    ax2.plot(mean_pred, frac_pos, "-o", color=PALETTE[3], lw=1.5)
    ax2.set_xlabel("predicted probability")
    ax2.set_ylabel("observed frequency")
    ax2.set_title("Model calibration (5-fold CV)")
    out = FIG / "classifier.png"
    fig.tight_layout(); fig.savefig(out); plt.close(fig)
    comp.to_csv(TAB / "rule_comparison.csv", index=False)
    return out.name, comp, cv


def _md_table(df: pd.DataFrame, floatfmt: str = "{:.3f}") -> str:
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return ""
            if float(v).is_integer() and abs(v) < 1e6:
                return str(int(v))
            return floatfmt.format(v)
        return str(v)
    header = "| " + " | ".join(df.columns) + " |"
    sep = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.to_numpy()]
    return "\n".join([header, sep, *rows])


def write_report(df, dta, dta_clin, detect_ladder, kappa_val, auc_tbl,
                 figs, lod_tbl, corr, ba, comp, cv):
    n_samples = df["sample_id"].nunique()
    n_path = df["pathogen"].nunique()
    lines = []
    lines.append("# Respiratory panel evaluation report\n")
    lines.append(
        f"Generated from a count matrix of {n_samples} samples across {n_path} "
        f"pathogens. All figures and tables below are produced by "
        "`scripts/run_evaluation.py`.\n"
    )

    lines.append("## 1. Detection accuracy\n")
    lines.append(
        "At clinically relevant load (>= 1e4 copies/mL) against all "
        "pathogen-absent observations:\n"
    )
    lines.append(_md_table(dta_clin.round(3)))
    lines.append(
        "\nPooled across the entire dilution series, including concentrations "
        "far below any clinical sample (this understates clinical sensitivity by "
        "design and is shown for completeness):\n"
    )
    lines.append(_md_table(dta.round(3)))
    lines.append(f"\nCohen's kappa (called vs truth, all pathogens pooled): **{kappa_val:.3f}**.\n")
    lines.append("\nDetection rate by input level, pooled across pathogens:\n")
    lines.append(_md_table(detect_ladder.round(3)))
    lines.append(f"\n![threshold sweep](figures/{figs['sweep']})\n")

    lines.append("## 2. Per-pathogen discrimination (ROC / AUC)\n")
    lines.append(_md_table(auc_tbl.round(3)))
    lines.append(f"\n![roc](figures/{figs['roc']})\n")

    lines.append("## 3. Limit of detection (LoD95)\n")
    lines.append(
        "Probit fit of detection rate against log10 input, per pathogen, with a "
        "delta-method 95% interval on the LoD.\n"
    )
    lines.append(_md_table(lod_tbl.round(3)))
    lines.append(f"\n![lod](figures/{figs['lod']})\n")

    lines.append("## 4. Quantitative agreement\n")
    lines.append(
        f"Read depth tracks input load with Pearson r = {corr['pearson_r']:.3f} "
        f"(Spearman {corr['spearman_r']:.3f}) on the log scale. Bland-Altman bias "
        f"{ba['bias']:.3f} log10, 95% limits of agreement "
        f"[{ba['loa_low']:.3f}, {ba['loa_high']:.3f}].\n"
    )
    lines.append(f"\n![agreement](figures/{figs['agreement']})\n")

    lines.append("## 5. Learned calling rule\n")
    lines.append(
        "A logistic model on read depth and genome breadth, scored by 5-fold "
        "cross-validation, against the depth-only threshold rule.\n"
    )
    lines.append(_md_table(comp.round(3)))
    coef_line = ", ".join(f"{k}={v:+.2f}" for k, v in cv["coefficients"].items())
    lines.append(f"\nStandardised coefficients: {coef_line}.\n")
    lines.append(f"\n![classifier](figures/{figs['classifier']})\n")

    (ROOT / "results" / "REPORT.md").write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", type=Path, default=ROOT / "data" / "example_counts.csv")
    args = ap.parse_args()

    _style()
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    df = load_counts(args.counts)

    dta = metrics.dta_metrics(df["present"].to_numpy(), df["called"].to_numpy())
    dta_clin = dta_at_clinical_load(df)
    detect_ladder = pooled_detection(df)
    kappa_val = metrics.kappa(df["present"].to_numpy(), df["called"].to_numpy())
    auc_tbl = metrics.roc_per_pathogen(df)
    dta.to_csv(TAB / "dta_metrics_pooled.csv", index=False)
    dta_clin.to_csv(TAB / "dta_metrics_clinical_load.csv", index=False)
    detect_ladder.to_csv(TAB / "detection_by_level.csv", index=False)
    auc_tbl.to_csv(TAB / "auc_per_pathogen.csv", index=False)

    sweep_fig = fig_threshold_sweep(df)
    roc_fig = fig_roc(df)
    lod_fig, lod_tbl = fig_lod(df)
    agr_fig, corr, ba = fig_agreement(df)
    clf_fig, comp, cv = fig_classifier(df)

    figs = dict(sweep=sweep_fig, roc=roc_fig, lod=lod_fig,
                agreement=agr_fig, classifier=clf_fig)
    write_report(df, dta, dta_clin, detect_ladder, kappa_val, auc_tbl,
                 figs, lod_tbl, corr, ba, comp, cv)

    print("evaluation complete")
    print("DTA at clinical load (>= 1e4 copies/mL):")
    print(dta_clin.to_string(index=False))
    print(f"\nkappa = {kappa_val:.3f}")
    print("\nLoD95 (copies/mL):")
    print(lod_tbl[["pathogen", "lod_copies", "ci_low", "ci_high"]].to_string(index=False))
    print(f"\nrule comparison:\n{comp.to_string(index=False)}")


if __name__ == "__main__":
    main()
