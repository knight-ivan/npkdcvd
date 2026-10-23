#!/usr/bin/env bash
# Reproduce every attribution result of the R2 revision in R (>= 4.1; packages:
# parallel, e1071).  Table 1 and Supplement Tables S.1-S.4 come from original_matlab/.  Run from any directory:
#   bash code/run_R2.sh          # full run, B=1000 -> ../results_R2 and ../figures
#   bash code/run_R2.sh smoke    # quick check, B=2  -> ../results_smoke only
# The two real-data analyses need the public data sets (see README.md):
#   ANURAN_CSV=/path/Frogs_MFCCs.csv  DIGIT_DIR=/dir/with/digit/binaries
set -euo pipefail
cd "$(dirname "$0")"
command -v Rscript >/dev/null || { echo 'ERROR: Rscript is not installed or not on PATH.' >&2; exit 2; }
if [[ -z "${RODEO_CORES:-}" ]]; then
  RODEO_CORES="$(Rscript -e 'cat(max(1L, parallel::detectCores()-2L))')"
fi
export RODEO_CORES
# One BLAS thread per forked worker (avoids CPU oversubscription under mclapply).
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
if [[ "${1:-}" == "smoke" ]]; then
  export RODEO_B=2 RODEO_OUT="$(cd .. && pwd)/results_smoke"
  echo "SMOKE TEST (B=2) -> $RODEO_OUT; never use these numbers in the paper"
else
  : "${RODEO_B:=1000}"
  export RODEO_B
fi

Rscript test_rodeo_core_R2.R                 # regression checks of the corrected core
Rscript rodeo_study1.R
Rscript rodeo_study2.R
Rscript rodeo_study3.R
Rscript rodeo_study5.R                       # exact-PIT primary + two sensitivities
for setting in E1 E2 E3 E4 E5 E6 E7; do
  Rscript robustness_experiments.R "$setting" --B="$RODEO_B" --rule=rodeo
done
Rscript robustness_e8_e9.R                   # E8 ablation, E9 no-rodeo baseline
Rscript accuracy_gate_audit.R                # reads the archived ../data/N1, N6 outputs
Rscript real_waveform.R                      # generated benchmark, no download needed
if [[ -f "${ANURAN_CSV:-}" ]]; then
  Rscript real_anuran.R
else
  echo "SKIP real_anuran.R: set ANURAN_CSV to Frogs_MFCCs.csv (UCI Anuran Calls)"
fi
if [[ -f "${DIGIT_DIR:-/tmp}/mnist_X.bin" && -f "${DIGIT_DIR:-/tmp}/digits8_X.bin" ]]; then
  Rscript real_digit_attribution.R
else
  echo "SKIP real_digit_attribution.R: run fetch_digit_data.py and set DIGIT_DIR"
fi
if [[ -z "${RODEO_OUT:-}" ]]; then           # paper figures come from the full run only
  Rscript study1_jaccard_figure.R
  if [[ -f ../results_R2/real_anuran_results.rds ]]; then Rscript frog_bar_fig.R; fi
fi
printf '\nResults saved under %s\n' "${RODEO_OUT:-../results_R2}"
