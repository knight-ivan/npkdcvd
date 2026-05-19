# npkdcvd

**Reproducibility bundle for:**

> Foo, H.-M., Chen, W.-P. N., and Chang, Y.-C. I. (2026).
> Class-Specific Variable Detection for Sparse Nonparametric Multiclass Classification.
> *Journal of Computational and Graphical Statistics* (submitted).

---

## Overview

This repository contains all R code, simulation scripts, and pre-generated figures needed to reproduce the tables and figures in the main paper and Supplementary Material.

The core method is **NPKDC-vd** (Nonparametric Posterior Kernel Density Classifier with class-specific variable detection): a sparse nonparametric multiclass classifier that applies Rodeo-style bandwidth shrinkage independently to each class-conditional density, producing a per-class bandwidth vector that serves as a direct attribution signal.

---

## Repository structure

```
npkdcvd/
├── reproduce_all_tables.R     # Master script — reproduces all paper tables
├── reproduce_all_figures.R    # Master script — reproduces all paper figures
├── R/                         # Core source code
│   ├── 00_setup.R             # Dependency checks and shared utilities
│   ├── 01_data_and_core_npkdc.R   # Data generation, NPKDC-vd, KDE classifier, metrics
│   ├── 02_benchmark_methods.R     # OvR-LASSO, OvE-LASSO, RBF-SVM, XGBoost/TreeSHAP
│   ├── 03_run_N4_modern_comparison_B1000.R  # Study 4 full runner (B=1000)
│   ├── 04_make_N4_figures.R   # Study 4 figure generation
│   ├── 05_smoke_test.R        # Quick end-to-end test (~30 seconds)
│   └── 07_make_phase3_visualizations.R  # Supplementary bandwidth box-plots
├── scripts/                   # Per-study simulation and figure scripts
│   ├── N1_simulation.R / N1_figure.R   # Synthetic Study 1
│   ├── N2_simulation.R / N2_figure.R   # Synthetic Study 3 (detection probability)
│   ├── N3_simulation.R / N3_figure.R   # Synthetic Study 2 (high-dimensional)
│   ├── N6_simulation.R / N6_figure_v2.R  # Synthetic Study 5 (generative vs discriminative)
│   ├── V1_jaccard.R           # Real Data 1: Anuran species classification
│   └── V2_digit_simulation.R  # Real Data 3: Handwritten digit attribution
└── figures/                   # Pre-generated EPS figures (committed for reference)
```

---

## Quick start

### 1. Smoke test (~30 seconds)
Verifies that the core NPKDC-vd functions run correctly:
```r
Rscript R/05_smoke_test.R
```

### 2. Reproduce all tables
```r
Rscript reproduce_all_tables.R
```
Full runs use B = 500–1000 replications and can take several hours.
For a quick check, open the script and set `B_OVERRIDE <- 5`.

### 3. Reproduce all figures
```r
Rscript reproduce_all_figures.R
```
Run after `reproduce_all_tables.R` (figures read from `results/`).

---

## Required R packages

The NPKDC-vd core runs in **base R only** (no external packages).

The full benchmark comparison additionally requires:
```r
install.packages(c("glmnet", "e1071", "xgboost", "parallel"))
```

| Package  | Used for |
|----------|----------|
| `glmnet` | OvR-LASSO and OvE-LASSO benchmarks |
| `e1071`  | Radial-kernel SVM benchmark |
| `xgboost`| XGBoost classifier and TreeSHAP attribution |
| `parallel`| Parallel class estimation (optional, speeds up full runs) |

If any benchmark package is missing, the corresponding method is skipped and NPKDC-vd results are still produced.

---

## Session information

All results in the paper were produced with R 4.3.x on macOS.
Full `sessionInfo()` output is archived in `results/sessionInfo.txt` after running the master scripts.

---

## Data sources

| Dataset | Source |
|---------|--------|
| Anuran calls (7 species, 22 MFCC features) | UCI Machine Learning Repository |
| UCI Waveform (5000 obs, 40 features) | UCI Machine Learning Repository |
| UCI Optical Digits (8×8, 64 features) | UCI Machine Learning Repository |
| MNIST (28×28, balanced 100/class) | LeCun et al. (2010) |

The real datasets are downloaded automatically by the simulation scripts on first run.

---

## Citation

```bibtex
@article{foo2026npkdcvd,
  author  = {Foo, Hui-Mean and Chen, Wan-Ping Nicole and Chang, Yuan-chin Ivan},
  title   = {Class-Specific Variable Detection for Sparse Nonparametric Multiclass Classification},
  journal = {Journal of Computational and Graphical Statistics},
  year    = {2026},
  note    = {Submitted}
}
```

---

## License

MIT License. See `LICENSE` for details.
