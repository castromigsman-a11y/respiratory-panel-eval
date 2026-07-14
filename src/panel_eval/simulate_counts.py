"""Generate a synthetic count matrix that mimics the output of a Nanopore
metagenomic respiratory panel, together with the ground truth.

The point of the simulator is not to model sequencing chemistry in detail. It
is to produce a count matrix whose statistical behaviour is realistic enough
that a genuine diagnostic evaluation (detection limits, DTA metrics, method
comparison) is meaningful and reproducible on any machine, without needing
reference genomes or basecalling.

Each sample spikes at most one pathogen at a known input concentration. The
pipeline reports, per sample, the number of reads assigned to every pathogen in
the panel plus a coarse "breadth" score (fraction of the target genome that
received any read), which the real minimap2 step also produces.

Signal model, per spiked pathogen p at input level L = log10(copies/mL):

    mean_reads(L) = T * 10 ** (slope * (L - c50_p))

so the expected read count crosses the calling threshold T near c50_p. Counts
are drawn from a negative binomial (Poisson-Gamma) to give the overdispersion
seen in real metagenomic data. Off-target and pathogen-free samples get a small
background mean, which produces occasional false positives and keeps
specificity below one, as in a real panel.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


@dataclass
class Pathogen:
    name: str
    c50_log10: float


def _neg_binomial(mean: np.ndarray, dispersion: float, rng: np.random.Generator) -> np.ndarray:
    """Draw negative-binomial counts via a Gamma-Poisson mixture.

    `dispersion` is the Gamma shape: small values give heavier tails. The mean
    is preserved exactly regardless of dispersion.
    """
    mean = np.maximum(np.asarray(mean, dtype=float), 1e-9)
    rate = rng.gamma(shape=dispersion, scale=mean / dispersion)
    return rng.poisson(rate)


def load_panel(config_path: str | Path) -> dict:
    with open(config_path) as handle:
        return yaml.safe_load(handle)


def simulate(
    config_path: str | Path,
    seed: int = 43,
    slope: float = 0.9,
    background_mean: float = 0.35,
    signal_dispersion: float = 2.0,
    background_dispersion: float = 0.4,
) -> pd.DataFrame:
    """Return a long-format count matrix with ground truth.

    Columns:
        sample_id      unique sample identifier
        pathogen       one row per pathogen per sample
        present        1 if this pathogen was spiked into this sample, else 0
        input_log10    log10 copies/mL if present, else NaN
        reads          reads assigned to this pathogen by the pipeline
        breadth        fraction of the target genome covered by >=1 read
    """
    cfg = load_panel(config_path)
    rng = np.random.default_rng(seed)

    pathogens = [Pathogen(p["name"], float(p["c50_log10"])) for p in cfg["pathogens"]]
    levels = list(cfg["input_levels_log10"])
    reps = int(cfg["replicates_per_level"])
    n_neg = int(cfg["negative_samples"])
    threshold = int(cfg["calling_threshold"])
    pileup_rate = float(cfg.get("pileup_rate", 0.05))

    rows: list[dict] = []

    def emit_sample(sample_id: str, spiked: str | None, input_log10: float | None) -> None:
        for p in pathogens:
            is_present = int(spiked == p.name)
            if is_present:
                # true signal: reads scale with input, breadth reflects genuine
                # coverage of the target genome and stays high once detected.
                mean = threshold * 10 ** (slope * (input_log10 - p.c50_log10))
                reads = int(_neg_binomial(np.array([mean]), signal_dispersion, rng)[0])
                breadth = float(np.clip(1 - np.exp(-reads / 8.0) + rng.normal(0, 0.07), 0.02, 0.99))
            elif rng.random() < pileup_rate:
                # spurious pileup: high depth from contamination or index hopping,
                # but reads concentrate in one locus so breadth stays low. This is
                # exactly the case depth-only calling gets wrong.
                reads = int(_neg_binomial(np.array([4 * threshold]), 1.5, rng)[0])
                breadth = float(np.clip(rng.beta(2.0, 12.0), 0.0, 0.45))
            else:
                # clean background: little or no signal, negligible breadth.
                reads = int(_neg_binomial(np.array([background_mean]), background_dispersion, rng)[0])
                breadth = float(np.clip(rng.beta(1.2, 60) + 0.01 * (reads > 0), 0, 0.2))
            rows.append(
                dict(
                    sample_id=sample_id,
                    pathogen=p.name,
                    present=is_present,
                    input_log10=(float(input_log10) if is_present else np.nan),
                    reads=reads,
                    breadth=round(breadth, 4),
                )
            )

    # dilution series: one spiked pathogen per sample
    for p in pathogens:
        for level in levels:
            for r in range(reps):
                sid = f"{p.name}_L{level}_r{r+1}"
                emit_sample(sid, p.name, level)

    # pathogen-free samples for the specificity side
    for i in range(n_neg):
        emit_sample(f"neg_{i+1:03d}", None, None)

    df = pd.DataFrame(rows)
    df["called"] = (df["reads"] >= threshold).astype(int)
    return df.reset_index(drop=True)


def write_example(config_path: str | Path, out_path: str | Path, seed: int = 43) -> pd.DataFrame:
    df = simulate(config_path, seed=seed)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return df


if __name__ == "__main__":
    here = Path(__file__).resolve().parents[2]
    frame = write_example(here / "config" / "panel.yaml", here / "data" / "example_counts.csv")
    print(f"wrote {len(frame)} rows, {frame['sample_id'].nunique()} samples")
