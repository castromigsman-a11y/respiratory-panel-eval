# respiratory-panel-eval

A Python framework for the analytical evaluation of a Nanopore metagenomic
pathogen panel. It takes a count matrix (reads and genome breadth per pathogen
per sample, with known ground truth) and computes the standard measures used to
characterise a diagnostic assay: limit of detection, diagnostic test accuracy
with confidence intervals, quantitative agreement, and a comparison between a
depth-only calling rule and a logistic rule that also uses genome breadth.

The example panel covers six targets: influenza A, influenza B, RSV,
rhinovirus, norovirus and SARS-CoV-2. All committed results are computed on
simulated data.

## What it measures

- **Detection accuracy:** sensitivity, specificity, PPV, NPV, F1 and accuracy
  with Wilson score intervals, reported at clinically relevant load
  (≥ 10⁴ copies/mL) and pooled across the full dilution series. Cohen's kappa
  for called vs true status.
- **Per-pathogen discrimination:** ROC curves and AUC, using read depth as the
  score.
- **Limit of detection (LoD95):** probit regression of detection rate on log10
  input over the dilution series, following the approach of CLSI EP17, with a
  delta-method 95% interval.
- **Quantitative agreement:** correlation between log read depth and log input
  load, and a Bland-Altman analysis.
- **Calling rule comparison:** a logistic regression on read depth and genome
  breadth, scored by 5-fold cross-validation and compared with the
  `reads >= threshold` rule at matched sensitivity.

The full output, with all tables and figures, is in
[results/REPORT.md](results/REPORT.md).

## Results on the example data

The example matrix contains 780 samples: a 10-level dilution series with 12
replicates per pathogen per level, plus 60 pathogen-free samples. These numbers
describe the simulator's behaviour under the parameters in `config/panel.yaml`;
they are not performance estimates for any real assay.

| Measure | Value |
|---|---|
| Sensitivity at ≥ 10⁴ copies/mL | 0.979 (95% CI 0.955–0.990) |
| Specificity | 0.970 (0.964–0.975) |
| PPV at ≥ 10⁴ copies/mL | 0.705 (0.659–0.748) |
| Sensitivity pooled over all dilution levels | 0.533 (0.497–0.570) |
| Per-pathogen AUC (read depth) | 0.73 (norovirus) to 0.90 (SARS-CoV-2) |
| LoD95 | 1.3 × 10³ (SARS-CoV-2) to 1.9 × 10⁴ (rhinovirus) copies/mL |
| Read depth vs input, Pearson r (log scale) | 0.85 |
| AUC, depth-only vs depth + breadth rule | 0.81 vs 0.92 |
| PPV at matched sensitivity (0.53), depth-only vs depth + breadth | 0.765 vs 0.997 |

![LoD probit fits](results/figures/lod_probit.png)

![Calling rule comparison](results/figures/classifier.png)

## Quickstart

```bash
pip install -r requirements.txt
python scripts/run_evaluation.py     # writes results/figures, results/tables, results/REPORT.md
pytest -q                            # tests for the statistical functions
```

Run with no arguments, the script regenerates the example count matrix if it is
missing and evaluates it. To evaluate a different matrix, for example one
produced by the pipeline layer:

```bash
python scripts/run_evaluation.py --counts data/counts.csv
```

The same steps are available through `make setup`, `make data`,
`make evaluate` and `make test`.

## Input data

The evaluation runs on a long-format count matrix with one row per pathogen per
sample. It can be produced in two ways:

1. **Simulated (default, no external tools).**
   `src/panel_eval/simulate_counts.py` draws read counts from a negative
   binomial model whose mean rises with input concentration, adds background
   reads, and injects spurious high-depth, low-breadth pileups as false
   positives. The committed results use this matrix.
2. **From simulated reads (optional).** The scripts in `pipeline/` download one
   reference per pathogen, simulate Nanopore reads with Badread over the
   dilution series, align with minimap2, and count reads and genome breadth per
   pathogen. The output uses the same schema. See
   [pipeline/README.md](pipeline/README.md).

## Repository layout

```
config/panel.yaml            panel definition and simulation parameters
src/panel_eval/              evaluation library
  simulate_counts.py         synthetic count matrix with ground truth
  metrics.py                 DTA metrics, ROC/AUC, kappa
  lod.py                     probit LoD95 with delta-method CI
  agreement.py               calibration, correlation, Bland-Altman
  classifier.py              cross-validated logistic calling rule
scripts/run_evaluation.py    runs everything, writes figures/tables/report
pipeline/                    optional Nanopore pipeline (reads to counts)
data/example_counts.csv      committed synthetic count matrix
results/                     committed figures, tables and REPORT.md
tests/test_metrics.py        tests
```

Changing the parameters in `config/panel.yaml` and re-running shows how each
metric responds. Random seeds are fixed, so results are reproducible.

## Limitations

All reported results come from simulated count data, not from sequencing of
real specimens, so they demonstrate the evaluation methods rather than the
performance of any assay. Several results are largely set by the simulator's
inputs: each pathogen's LoD95 follows from the `c50_log10` value assigned to it
in `config/panel.yaml`, and the advantage of the breadth-aware calling rule
reflects the fact that simulated false-positive pileups are generated with low
genome breadth by construction; real contamination, index hopping or
cross-mapping between related genomes may not separate as cleanly. The
simulation does not model nucleic acid extraction, matrix effects, sample
degradation, co-infections (each sample contains at most one pathogen), strain
diversity, or batch and run effects. Confidence intervals treat each
pathogen-sample row as an independent observation, although the six rows from
one sample would be correlated in real data, so the intervals may be too
narrow. PPV and NPV depend on the proportion of positives in the design and do
not transfer to clinical prevalence. The ≥ 10⁴ copies/mL "clinically relevant"
cut-off is a fixed assumption, and the calling-rule comparison is reported at
the pooled sensitivity of 0.53. The optional read-level pipeline has not been
run on clinical samples. Norovirus is included as a panel target although it is
primarily an enteric rather than a respiratory pathogen.

## License

MIT. See [LICENSE](LICENSE).
