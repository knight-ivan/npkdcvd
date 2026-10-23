# Original MATLAB implementation (classification and runtime comparisons)

These files belong to the original study.  They produced the classification and
runtime comparisons that the paper carries over unchanged: main-paper Table 1 and
Supplement Tables S.1–S.4 and Figure S.3.  The files are included as found
(2019–2020 versions); all attribution analyses of the paper use the R scripts in the
parent folder.

## Files

| File | Role |
|---|---|
| `RodeoPCNew.m` | Original NPKDC classifier: per-class rodeo bandwidths, training and prediction |
| `Rodeo_localNew.m`, `Rodeo_localNormal.m`, `Rodeo_GradNormal.m` | Rodeo bandwidth routines called by `RodeoPCNew.m` |
| `Ex1Com.m` | Study 1 design (cyclic relevant sets; relevant N(0.5, (0.02k)^2), irrelevant Uniform(0,1)); NPKDC-vd vs. LSPC accuracy and time, bandwidth Z-scores |
| `Ex2ComDistNew.m` | Class-combination design (2–5 of 5 classes, d = 10), a later variant of the design behind Supplement Tables S.3–S.4 |
| `ExAll.m` | Batch driver (MATLAB `batch`, 8 workers) with example settings |

Table 1 corresponds to `Ex1Com(10, 10, 1000, 30, 6, 150, 100, 1, 0)` (10 classes,
d = 30, 6 relevant variables, 150 training and 100 test points per class, 1000
replications).  Random seeds were not fixed in these programs, so reruns reproduce
the results only up to Monte Carlo error.

## Dependencies

* MATLAB with the Statistics and Machine Learning Toolbox (`normrnd`, `randsample`,
  `combnk`, `zscore`) and the Parallel Computing Toolbox (`parfor`, `batch`).
* The LSPC classifier (`LSPC_train.m`, `LSPC_test.m`) by M. Sugiyama, distributed by
  its author (originally at http://sugiyama-www.cs.titech.ac.jp/~sugi/software/LSPC/).
  It is third-party code and is not redistributed here.

## Archived outputs

`archived_outputs/` holds the tables exactly as these programs wrote them.

| File | Reported as |
|---|---|
| `Ex1_acc.tex` | Table 1, NPKDC-vd and LSPC columns |
| `Ex1SVM_acc.tex` | Table 1, SVM and SVMfs columns |
| `Ex1_prob.tex` | Supplement Table S.1 |
| `Ex1BW.tex` | Supplement Table S.2 |
| `Ex3DiffN0904.tex` | Supplement Table S.3 |
| `Ex3DiffNTime0904.tex` | Supplement Table S.4 |

Supplement Figure S.3 (`figures/Ex1.eps`) comes from the same Study 1 runs.

## Not preserved

* The SVM and SVMfs comparators in Table 1 used MATLAB's error-correcting output codes
  (ECOC) support-vector classifier with a Gaussian kernel and sequential feature
  selection; that driver script was not preserved.
* The 2018 script for Supplement Tables S.3–S.4 (five classes, d = 10, two relevant
  coordinates with standard deviations 0.1 and 0.2, eight Uniform(0,1) coordinates,
  training sizes 50–1000) was not preserved; `Ex2ComDistNew.m` implements the closely
  related later version of this class-combination design.

Some comments in the MATLAB files are written in Traditional Chinese (Big5 encoding);
this does not affect execution.
