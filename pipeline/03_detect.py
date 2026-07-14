"""Align reads and produce the count matrix consumed by the evaluation layer.

For each sample in the sample sheet, map reads to the concatenated panel
reference with minimap2 (map-ont preset), then for every pathogen count the
reads whose primary alignment lands on its contig and compute genome breadth
(fraction of the reference covered by at least one read). The output columns
match data/example_counts.csv exactly, so scripts/run_evaluation.py runs on
real pipeline output without modification.

Requires minimap2, samtools and pysam.

    python pipeline/03_detect.py --reads reads --references references --out data/counts.csv
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def align(fastq: Path, reference: Path, bam: Path) -> None:
    sam = bam.with_suffix(".sam")
    with open(sam, "w") as out:
        subprocess.run(["minimap2", "-ax", "map-ont", str(reference), str(fastq)],
                       check=True, stdout=out)
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    sam.unlink()


def counts_and_breadth(bam: Path, contig_lengths: dict[str, int]) -> dict[str, tuple[int, float]]:
    import pysam

    result = {name: [0, set()] for name in contig_lengths}
    with pysam.AlignmentFile(bam, "rb") as handle:
        for read in handle.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            name = read.reference_name
            if name in result:
                result[name][0] += 1
                # sample covered positions coarsely for breadth
                for pos in range(read.reference_start, read.reference_end, 50):
                    result[name][1].add(pos)
    out = {}
    for name, (n_reads, covered) in result.items():
        length = contig_lengths[name]
        breadth = min(1.0, (len(covered) * 50) / length) if length else 0.0
        out[name] = (n_reads, round(breadth, 4))
    return out


def main(reads: Path, references: Path, out: Path) -> None:
    import pysam

    cfg = yaml.safe_load((ROOT / "config" / "panel.yaml").read_text())
    threshold = int(cfg["calling_threshold"])
    panel_ref = references / "panel_references.fasta"
    contig_lengths = dict(zip(pysam.FastaFile(str(panel_ref)).references,
                              pysam.FastaFile(str(panel_ref)).lengths))

    sheet = pd.read_csv(reads / "samples.csv")
    rows = []
    bam_dir = reads / "bam"
    bam_dir.mkdir(exist_ok=True)
    for _, s in sheet.iterrows():
        bam = bam_dir / f"{s['sample_id']}.bam"
        align(Path(s["fastq"]), panel_ref, bam)
        measured = counts_and_breadth(bam, contig_lengths)
        for pathogen in contig_lengths:
            n_reads, breadth = measured[pathogen]
            present = int(pathogen == s["pathogen"])
            rows.append(dict(
                sample_id=s["sample_id"],
                pathogen=pathogen,
                present=present,
                input_log10=(float(s["input_log10"]) if present else float("nan")),
                reads=n_reads,
                breadth=breadth,
                called=int(n_reads >= threshold),
            ))
        print(f"processed {s['sample_id']}")

    df = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", type=Path, default=ROOT / "reads")
    ap.add_argument("--references", type=Path, default=ROOT / "references")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "counts.csv")
    args = ap.parse_args()
    main(args.reads, args.references, args.out)
