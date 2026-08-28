# NPKDC-vd: Generative Variable Attribution for Sparse Nonparametric Multi-Class Classification

Reproducibility code for **NPKDC-vd** ("Naming Each Class"), a method that names, for
each class, the variables on which that class's *own* conditional density is
**concentrated** — its generative identity — rather than the variables that merely
*separate* it from the other classes.

This bundle implements the **corrected density-rodeo**: a standardized
density-derivative statistic `Z_j = ∂f̂/∂h_j`, cross-fitting, a squared
(non-cancelling) aggregate `T_j`, and a deterministic largest-gap selection rule, with
an optional marginal pre-screen for high dimension. The method core is
[`code/rodeo_core.R`](code/rodeo_core.R).

## Requirements
- **R** with `parallel` (fork-based; macOS/Linux). `e1071` for the SVM comparison.
- **Python 3** with `scikit-learn` + `numpy` — only to fetch the digit datasets.

## Reproduce

Run from `code/`. `RODEO_B` = number of splits (use `1000` for the final tables);
`RODEO_CORES` = core count.

```bash
cd code

# --- Simulation studies ---
RODEO_B=1000 Rscript rodeo_study1.R     # Study 1  overlapping sets / Jaccard
RODEO_B=1000 Rscript rodeo_validate.R   # Study 2  d = 200
RODEO_B=1000 Rscript rodeo_study3.R     # Study 3  convergence (Theorem 2)
RODEO_B=1000 Rscript rodeo_study4.R     # Study 4  class- vs per-observation stability
RODEO_B=1000 Rscript rodeo_study5.R     # Study 5  generative vs discriminative

# --- Robustness grid E1-E9 (+ ablation and a density-based competitor) ---
Rscript robustness_experiments.R all --B=500 --rule=rodeo
Rscript rodeo_robustness_compare.R

# --- Real data ---
# Digits (UCI optdigits 8x8 + MNIST) from OpenML, written to /tmp:
python3 fetch_digit_data.py /tmp
# Anuran Calls (MFCCs) from the UCI ML Repository:
#   curl -L -o anuran.zip "https://archive.ics.uci.edu/static/public/406/anuran+calls+mfccs.zip"
#   unzip anuran.zip   ->   Frogs_MFCCs.csv
ANURAN_CSV=/path/to/Frogs_MFCCs.csv RODEO_B=1000 Rscript real_anuran.R
RODEO_B=1000 Rscript real_waveform.R
Rscript V2_digit_simulation.R                    # digit classification
RODEO_B=1000 Rscript real_digit_attribution.R    # MNIST per-pixel attribution
Rscript frog_bar_fig.R                            # first-MFCC selection-frequency figure
```

## Why the handwritten-digit attribution is (correctly) uniform
`mnist_scoremap.R`, `uci8x8_attr.R`, and `mnist_ntr_sensitivity.R` document that the
digit attribution is a near-uniform per-pixel profile: **no small set of pixels
dominates** a digit's density, because a digit's identity is carried by the whole
*distributed* stroke configuration. The same rule that isolates a sparse dominant set
where one exists — the simulation studies, and the leading MFCC that anchors several
frog species — correctly reports its absence where none does.

## Layout
```
code/    all R scripts and fetch_digit_data.py; rodeo_core.R is the method core
```
Generated results (`*.rds`) and figures are produced by the scripts and are not stored
in the repository (see `.gitignore`). The associated manuscript is under review.
