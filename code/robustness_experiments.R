# robustness_experiments.R
# ---------------------------------------------------------------------------
# Theme-E robustness studies requested by the reviewers (R1.C5, R2.C6).
# SCAFFOLD: data-generating processes E1-E7 implemented; E8 (ablation) and
# E9 (density competitor) are stubbed. Sources the shared detector from
# npkdcvd_core.R so all studies use ONE implementation.
#
# Run one scenario:   Rscript robustness_experiments.R E4
# Run all:            Rscript robustness_experiments.R all
# Smoke test (fast):  Rscript robustness_experiments.R E4 --B=20
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(parallel)
}))
`%||%` <- function(a, b) if (is.null(a)) b else a
.get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(.get_script_dir(), "npkdcvd_core.R"))
source(file.path(.get_script_dir(), "rodeo_core.R"))   # Route-2 rodeo detectors

# ---- shared helpers --------------------------------------------------------
# One class's n x D design: relevant cols ~ N(mu_rel, sd_rel); others via irrel_fun.
make_class <- function(n, D, rel, mu_rel, sd_rel, irrel_fun) {
  X <- matrix(0, n, D)
  nr <- length(rel)
  X[, rel] <- matrix(rnorm(n * nr, mean = rep(mu_rel, each = n), sd = sd_rel),
                     n, nr)
  X[, -rel] <- irrel_fun(n, D - nr)
  X
}

disjoint_rel <- function(C, r_vec) {
  starts <- cumsum(c(1L, head(r_vec, -1)))
  lapply(seq_len(C), function(y) starts[y]:(starts[y] + r_vec[y] - 1L))
}

# irrelevant-column generators
irr_uniform <- function(n, m) matrix(runif(n * m), n, m)
irr_normal  <- function(n, m) matrix(rnorm(n * m), n, m)
irr_t3      <- function(n, m) matrix(rt(n * m, df = 3), n, m)
irr_ar1 <- function(n, m, rho = 0.5) {
  Z <- matrix(0, n, m); Z[, 1] <- rnorm(n)
  for (j in seq_len(m)[-1]) Z[, j] <- rho * Z[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
  Z
}

# ---- scenario registry -----------------------------------------------------
# Each scenario: list(C, D, r, n_grid, gen=function(n_vec) -> list(train,test,rel_sets))
# n_vec is length-C per-class training sizes (allows imbalance).
SCEN <- list()

# E1: shrinking bandwidth gap — relevant sd grows toward background as n_y rises,
#     so delta_y -> 0. Probes the finite-sample n_y >~ d log d / delta^2 boundary.
SCEN$E1 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(50, 100, 200, 400, 800),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    # Gap between relevant and background bandwidths SHRINKS as n grows (problem
    # gets harder with more data). Constant needs calibration to land the
    # P(exact) transition inside n_grid; direction is what matters for the study.
    sd_rel <- pmin(0.9, 1 - 1.2 * mean(nv)^(-0.4))
    tr <- lapply(seq_len(C), function(y) make_class(nv[y], D, rel[[y]], rep(0.5, 3), sd_rel, irr_normal))
    te <- lapply(seq_len(C), function(y) make_class(100,   D, rel[[y]], rep(0.5, 3), sd_rel, irr_normal))
    list(train = tr, test = te, rel_sets = rel)
  })

# E2: approximate sparsity — besides strong relevant vars, a block of WEAK vars
#     (intermediate sd) that are not in the truth set; tests over/under-selection.
SCEN$E2 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(100, 200, 400),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    weak <- 16:19  # weakly-structured, sd between relevant(0.1) and background(1)
    mkirr <- function(n, m) irr_normal(n, m)
    mk <- function(n, y) {
      X <- make_class(n, D, rel[[y]], rep(0.5, 3), 0.1, mkirr)
      X[, weak] <- matrix(rnorm(n * length(weak), 0, 0.5), n, length(weak))
      X
    }
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E3: class imbalance — unequal n_y across classes.
SCEN$E3 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(200),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    base <- nv[1]; nv <- round(base * c(0.4, 0.7, 1.0, 1.5, 2.0))
    mk <- function(n, y) make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_normal)
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E4: high dimension d > n_y.
SCEN$E4 <- list(C = 5, D = 300, r = rep(3, 5), n_grid = c(150),
  gen = function(nv) {
    C <- 5; D <- 300; rel <- disjoint_rel(C, rep(3, C))
    mk <- function(n, y) make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_normal)
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E5: correlation between a relevant and an irrelevant variable (per class).
SCEN$E5 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(150, 300),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C)); rho <- 0.7
    mk <- function(n, y) {
      X <- make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_normal)
      tgt <- setdiff(seq_len(D), unlist(rel))[1]        # an irrelevant col
      X[, tgt] <- rho * scale(X[, rel[[y]][1]]) + sqrt(1 - rho^2) * rnorm(n)
      X
    }
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E6: heavy-tailed irrelevant variables (t_3).
SCEN$E6 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(150, 300),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    mk <- function(n, y) make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_t3)
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E7: non-flat nuisance SHARED across classes + a "trap" column that is non-flat
#     but IDENTICAL across classes (must NOT be flagged) — direct test of the
#     concentration definition (R2.C2 / Theme A).
SCEN$E7 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(150, 300),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    trap <- 20L   # non-flat (bimodal) but same law in every class -> irrelevant
    mk <- function(n, y) {
      X <- make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_normal)
      X[, trap] <- ifelse(runif(n) < 0.5, rnorm(n, -2, 0.3), rnorm(n, 2, 0.3))
      X
    }
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)  # trap NOT in any rel set: success = never flagged
  })

# E8: ablation (KDE-no-detector; z-score vs gap rule) — TODO: wire in core rules.
# E9: density-based competitor (per-class marginal-spread ranking) — TODO.

# ---- runner ----------------------------------------------------------------
run_scenario <- function(name, B = 500, rule = c("rodeo", "rodeo_bc", "gap", "zscore"),
                         cores = max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores() - 2L))))) {
  rule <- match.arg(rule)
  sc <- SCEN[[name]]; if (is.null(sc)) stop("unknown scenario: ", name)
  # rodeo = Route-2 density rodeo (default); rodeo_bc = baseline-corrected (PIT);
  # gap/zscore = the marginal-sd concentration detector = density-based COMPETITOR
  # (E9) and the "KDE-no-rodeo" ablation baseline (E8).
  detect <- switch(rule,
    rodeo    = function(tr) rodeo_detect(tr,    keep = 40, mc.cores = 1L),
    rodeo_bc = function(tr) rodeo_detect_bc(tr, keep = 40, mc.cores = 1L),
    gap      = function(tr) detect_gap(tr),
    zscore   = function(tr) detect_zscore(tr))
  out <- list()
  for (n0 in sc$n_grid) {
    reps <- mclapply(seq_len(B), function(b) {
      set.seed(b * 1000 + n0)
      d  <- sc$gen(rep(n0, sc$C))
      at <- attr_metrics(detect(d$train), d$rel_sets)
      colMeans(at)
    }, mc.cores = cores)
    M <- do.call(rbind, reps)
    out[[as.character(n0)]] <- colMeans(M)
    cat(sprintf("[%s] n0=%d rule=%s  Prec=%.3f Rec=%.3f F1=%.3f P(exact)=%.3f\n",
                name, n0, rule, out[[as.character(n0)]]["prec"],
                out[[as.character(n0)]]["rec"], out[[as.character(n0)]]["f1"],
                out[[as.character(n0)]]["exact"]))
  }
  out
}

# ---- CLI -------------------------------------------------------------------
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) {
    which <- args[1]
    Bopt  <- as.integer(sub("--B=", "", grep("--B=", args, value = TRUE)))
    B     <- if (length(Bopt) && !is.na(Bopt)) Bopt else 500
    Ropt  <- sub("--rule=", "", grep("--rule=", args, value = TRUE))
    rule  <- if (length(Ropt)) Ropt[1] else "rodeo"   # Route-2 default
    names <- if (which == "all") setdiff(names(SCEN), c("E8", "E9")) else which
    res <- lapply(names, run_scenario, B = B, rule = rule)
    names(res) <- names
    saveRDS(res, file.path("..", "data",
            sprintf("robustness_%s_%s.rds", if (which == "all") "all" else which, rule)),
            compress = FALSE)
  }
}
