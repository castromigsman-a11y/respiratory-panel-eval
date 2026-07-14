# Pipeline layer (optional, reproduces the upstream data)

The evaluation in `scripts/run_evaluation.py` runs on a count matrix and needs
no bioinformatics tools. This folder reproduces the steps that *produce* such a
matrix from simulated Nanopore reads, so the whole chain is transparent.

It is optional. The committed results use the synthetic count matrix in
`data/example_counts.csv`. Run the pipeline only if you want to regenerate the
matrix from reads.

## Tools

Install with the provided conda environment:

```bash
conda env create -f ../environment.yml
conda activate respiratory-panel-eval
```

This brings in Badread, minimap2, samtools, pysam and Biopython.

## Steps

```bash
# 1. download one reference per pathogen (needs an email for the NCBI API)
python 01_fetch_references.py --email you@example.org

# 2. simulate the dilution series; supply a human background FASTA
#    (a single chromosome or a GRCh38 subset is enough for a demo)
python 02_simulate_reads.py --human path/to/human_subset.fasta

# 3. align and build the count matrix
python 03_detect.py --out ../data/counts.csv

# 4. evaluate the real matrix with the same script used for the example
python ../scripts/run_evaluation.py --counts ../data/counts.csv
```

## Output schema

`03_detect.py` writes exactly the columns the evaluation layer expects:

| column | meaning |
| --- | --- |
| sample_id | unique sample identifier |
| pathogen | one row per pathogen per sample |
| present | 1 if this pathogen was spiked, else 0 (ground truth) |
| input_log10 | log10 copies/mL if present, else empty |
| reads | primary alignments to the pathogen contig |
| breadth | fraction of the pathogen genome covered by >= 1 read |
| called | 1 if reads >= calling_threshold |

Because the schema matches, the synthetic and real matrices are
interchangeable inputs to the evaluation.
