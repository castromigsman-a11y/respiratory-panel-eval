"""Simulate Nanopore reads for the dilution series with Badread.

For each pathogen and each input level in config/panel.yaml, generate reads
from the pathogen reference at a coverage that stands in for input load, mix
them with a fixed background of human reads, and write one FASTQ per sample.
A sample sheet records the ground truth (which pathogen, which level) so the
detection step can label its output.

Requires Badread and a human reference FASTA. Maps log10 input level to target
pathogen coverage; the background coverage is held constant so that increasing
input raises the pathogen fraction, mirroring a real dilution series.

    python pipeline/02_simulate_reads.py \
        --references references --human GRCh38_subset.fasta --out reads
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

# log10(copies/mL) -> target pathogen coverage (x). Tuned so the calling
# threshold lands near the configured c50 values after alignment.
COVERAGE_BY_LEVEL = {1.0: 0.02, 1.5: 0.05, 2.0: 0.1, 2.5: 0.3, 3.0: 1.0,
                     3.5: 3.0, 4.0: 8.0, 4.5: 20.0, 5.0: 50.0, 6.0: 200.0}
BACKGROUND_COVERAGE = "2x"
READ_LENGTH = "5000,4000"          # mean, stdev
ERROR_MODEL = "nanopore2023"


def badread(reference: Path, quantity: str, seed: int) -> str:
    cmd = ["badread", "simulate", "--reference", str(reference),
           "--quantity", quantity, "--length", READ_LENGTH,
           "--error_model", ERROR_MODEL, "--seed", str(seed)]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def main(references: Path, human: Path, out: Path) -> None:
    cfg = yaml.safe_load((ROOT / "config" / "panel.yaml").read_text())
    reps = int(cfg["replicates_per_level"])
    out.mkdir(parents=True, exist_ok=True)
    sheet = out / "samples.csv"

    with open(sheet, "w") as manifest:
        manifest.write("sample_id,pathogen,input_log10,fastq\n")
        seed = 0
        for p in cfg["pathogens"]:
            ref = references / f"{p['name']}.fasta"
            for level in cfg["input_levels_log10"]:
                cov = COVERAGE_BY_LEVEL[float(level)]
                for r in range(reps):
                    seed += 1
                    sid = f"{p['name']}_L{level}_r{r+1}"
                    fastq = out / f"{sid}.fastq"
                    reads = badread(ref, f"{cov}x", seed)
                    reads += badread(human, BACKGROUND_COVERAGE, seed + 100000)
                    fastq.write_text(reads)
                    manifest.write(f"{sid},{p['name']},{level},{fastq}\n")
                    print(f"simulated {sid}")
    print(f"wrote reads and sample sheet to {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--references", type=Path, default=ROOT / "references")
    ap.add_argument("--human", type=Path, required=True, help="human background FASTA")
    ap.add_argument("--out", type=Path, default=ROOT / "reads")
    args = ap.parse_args()
    main(args.references, args.human, args.out)
