"""Download reference genomes for the panel.

Fetches one RefSeq record per pathogen listed in config/panel.yaml and writes a
concatenated reference used later for read simulation and alignment. The human
background reference is large and licensing varies, so it is not downloaded
here: pass its path with --human when simulating reads.

Requires network access to NCBI and Biopython. This script is part of the
reproducible pipeline and is not exercised by the evaluation layer, which works
directly from a count matrix.

    python pipeline/01_fetch_references.py --email you@example.org
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def fetch(email: str, out_dir: Path) -> None:
    from Bio import Entrez, SeqIO

    Entrez.email = email
    cfg = yaml.safe_load((ROOT / "config" / "panel.yaml").read_text())
    out_dir.mkdir(parents=True, exist_ok=True)

    combined = out_dir / "panel_references.fasta"
    with open(combined, "w") as bundle:
        for p in cfg["pathogens"]:
            acc = p["accession"]
            print(f"fetching {p['name']} ({acc})")
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            # rename the contig to the pathogen so downstream counting is simple
            record.id = p["name"]
            record.description = acc
            SeqIO.write(record, out_dir / f"{p['name']}.fasta", "fasta")
            SeqIO.write(record, bundle, "fasta")
            time.sleep(0.4)  # be polite to the NCBI API
    print(f"wrote {combined}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="required by the NCBI Entrez API")
    ap.add_argument("--out", type=Path, default=ROOT / "references")
    args = ap.parse_args()
    fetch(args.email, args.out)
