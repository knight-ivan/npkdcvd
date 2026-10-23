# NPKDC-vd: Generative Variable Attribution for Sparse Nonparametric Multi-Class Classification

Reproducibility code for NPKDC-vd ("Naming Each Class").  For each class, NPKDC-vd
identifies the variables whose class-specific marginal departs from a shared
background (a generative attribution profile), as distinct from the variables that
separate the class from the others.

The method core is [`code/rodeo_core.R`](code/rodeo_core.R): a cross-fitted,
studentized density-derivative (rodeo) statistic with a deterministic largest-gap
selection rule, an optional reference (PIT) transform for a non-flat shared
background, and an optional marginal pre-screen.  `code/test_rodeo_core_R2.R` checks it.

## Requirements

* R (>= 4.1) with `parallel` and `e1071`.
* Python 3 with `scikit-learn` and `numpy`, only to download the digit data.
* MATLAB, only for the original classification comparisons in `code/original_matlab/`.

## Reproduce

```bash
bash code/run_R2.sh          # every result of the revised paper, B = 1000
bash code/run_R2.sh smoke    # quick check, B = 2; writes results_smoke/ only
```

The real-data steps need public data sets:

```bash
curl -L -o anuran.zip "https://archive.ics.uci.edu/static/public/406/anuran+calls+mfccs.zip"
unzip anuran.zip                                   # Frogs_MFCCs.csv
python3 code/fetch_digit_data.py /path/to/digits   # UCI optdigits and MNIST
ANURAN_CSV=$PWD/Frogs_MFCCs.csv DIGIT_DIR=/path/to/digits bash code/run_R2.sh
```

Outputs are written to `results_R2/` (not stored in the repository).

## Contents of `code/`

| Files | Purpose |
|---|---|
| `run_R2.sh` | Entry point for all results |
| `rodeo_core.R`, `npkdcvd_core.R`, `reviewer_metrics.R` | Method core and shared helpers |
| `rodeo_study1.R`, `rodeo_study2.R`, `rodeo_study3.R`, `rodeo_study5.R` | Simulation Studies 1, 2, 3, and 5 |
| `robustness_experiments.R`, `robustness_e8_e9.R` | Robustness designs E1–E7, ablation E8, no-rodeo baseline E9 |
| `real_anuran.R`, `real_waveform.R`, `real_digit_attribution.R` | Real-data and benchmark attribution (pooled-PIT and raw scale) |
| `study1_jaccard_figure.R`, `frog_bar_fig.R` | Figures |
| `N*_simulation.R`, `N*_figure*.R`, `V1_jaccard.R`, `V2_digit_simulation.R` | Earlier non-rodeo analyses that the paper still reports |
| `original_matlab/` | Original MATLAB implementation behind Table 1 and Supplement Tables S.1–S.4 |
| `rodeo_study4.R`, `rodeo_validate.R`, `rodeo_validate_bc.R`, `mnist_ntr_sensitivity.R` | Development checks |

## Change log

* **Second revision (October 2026).** Corrected a coordinate-alignment error in
  `rodeo_core.R` (evaluation columns are now restricted to the active coordinates; a
  regression test was added) and regenerated every rodeo-based result.  All analyses
  run in R through `run_R2.sh`; E8 and E9 are saved outputs; the real-data analyses
  report the pooled-PIT selector with the raw scale as a sensitivity; the original
  MATLAB code was added.
* **First revision (August 2026).** Corrected density rodeo with cross-fitting,
  studentization, and the deterministic largest-gap rule.

The associated manuscript is under review.
