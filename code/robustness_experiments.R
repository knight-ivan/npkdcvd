# robustness_experiments.R
# ---------------------------------------------------------------------------
# Theme-E robustness studies requested by the reviewers (R1.C5, R2.C6).
# Data-generating processes E1-E7 are implemented here.  E8 (baseline-correction
# ablation) and E9 (no-rodeo density baseline) are produced by robustness_e8_e9.R,
# which calls run_scenario() below on the same generators and seeds.
# Sources the shared detector from npkdcvd_core.R and the rodeo from rodeo_core.R.
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
source(file.path(.get_script_dir(), "reviewer_metrics.R"))

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

# E1: weakening dispersion contrast (NOT asymptotically vanishing signal).
# Relevant mean remains 0.5 and sd is capped at 0.9, so a nonzero
# difference from the N(0,1) background persists at every sample size.
SCEN$E1 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(50, 100, 200, 400, 800),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    # Increasing relevant SD attenuates but does NOT eliminate the signal.
    # This is not evidence for a shrinking-to-zero asymptotic signal regime.
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
SCEN$E3 <- list(C = 5, D = 25, r = rep(3, 5), n_grid = c(100, 200, 300, 400),
  gen = function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C, rep(3, C))
    base <- nv[1]; nv <- round(base * c(0.4, 0.7, 1.0, 1.5, 2.0))
    mk <- function(n, y) make_class(n, D, rel[[y]], rep(0.5, 3), 0.15, irr_normal)
    list(train = lapply(seq_len(C), function(y) mk(nv[y], y)),
         test  = lapply(seq_len(C), function(y) mk(100,   y)),
         rel_sets = rel)
  })

# E4: high dimension d > n_y.
SCEN$E4 <- list(C = 5, D = 300, r = rep(3, 5), n_grid = c(75, 100, 150, 200),
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

# Optional E1V: a genuinely VANISHING class-specific signal (new R2 analysis).
# Here Delta_n = 2 n^(-1/4), mu_rel(n)=Delta_n and sd_rel(n)=1-Delta_n,
# so both location and dispersion contrasts vanish as n grows. At the
# smallest n=50, Delta_n=.752 < 1. This is a NEW DGP and must not be confused
# with the historical, bounded-contrast E1 estimates in the manuscript.
SCEN$E1V <- list(C=5, D=25, r=rep(3,5), n_grid=c(50,100,200,400,800),
  gen=function(nv) {
    C <- 5; D <- 25; rel <- disjoint_rel(C,rep(3,C))
    mk <- function(n,y) {
      delta <- 2 * nv[y]^(-0.25)  # fixed train/test DGP at each grid point
      make_class(n,D,rel[[y]],rep(delta,3),1-delta,irr_normal)
    }
    list(train=lapply(seq_len(C),function(y) mk(nv[y],y)),
         test=lapply(seq_len(C),function(y) mk(100,y)), rel_sets=rel)
  })

# E8/E9 are comparison summaries rather than additional DGPs; see
# robustness_e8_e9.R for baseline OFF/ON and the no-rodeo bandwidth baseline.

# ---- updated R2 runner -----------------------------------------------------
# Report raw class-by-replication outcomes AND independent-replicate MC intervals.
# Data columns include the class-specific actual n_y (E3 unequal sizes), plus
# sure-screening and conditional selection success for screened algorithms.
# Return value: named list by nominal n0 with vectors c(prec,rec,f1,exact).
run_scenario <- function(name, B = 500,
                         rule = c("rodeo", "rodeo_bc", "gap", "zscore"),
                         cores = reviewer_cores(), save_audit = TRUE) {
  rule <- match.arg(rule)
  sc <- SCEN[[name]]
  if (is.null(sc)) stop("unknown scenario: ", name)
  if (!is.finite(B) || B < 1L) stop("B must be positive")
  out <- list(); audit <- list()
  for (n0 in sc$n_grid) {
    reps <- parallel::mclapply(seq_len(B), function(b) {
      set.seed(b * 1000L + n0 + 100000L * match(name,names(SCEN)))
      d <- sc$gen(rep(n0, sc$C))
      # E3 changes the actual sample size by class; retain these values.
      actual_n <- vapply(d$train, nrow, integer(1))
      if (rule == "rodeo") {
        obj <- rodeo_detect_detail(d$train, keep=40, mc.cores=1L)
        det <- obj$detected; screens <- obj$screen
      } else if (rule == "rodeo_bc") {
        pool <- do.call(rbind,d$train)
        Xpit <- lapply(d$train,function(X) pit_transform(X,pool))
        obj <- rodeo_detect_detail(Xpit, keep=40, mc.cores=1L)
        det <- obj$detected; screens <- obj$screen
      } else {
        det <- if (rule == "gap") detect_gap(d$train) else detect_zscore(d$train)
        screens <- replicate(sc$C, seq_len(sc$D), simplify=FALSE)
      }
      rows <- reviewer_metric_rows(det,d$rel_sets,b,name,n0)
      rows$rule <- rule
      rows$actual_n_y <- actual_n
      rows$screen_all <- vapply(seq_len(sc$C),function(y)
        as.numeric(all(d$rel_sets[[y]] %in% screens[[y]])),numeric(1))
      rows$screened_dimension <- vapply(screens,length,integer(1))
      rows$exact_given_screen <- ifelse(rows$screen_all == 1,rows$exact,NA_real_)
      rows
    }, mc.cores = if (.Platform$OS.type == "windows") 1L else cores,
       mc.set.seed = FALSE)
    R <- do.call(rbind,reps)
    audit[[as.character(n0)]] <- R
    m <- colMeans(R[,c("precision","recall","f1","exact"),drop=FALSE])
    out[[as.character(n0)]] <- setNames(unname(m),c("prec","rec","f1","exact"))
    cat(sprintf("[%s] n0=%d rule=%s B=%d F1=%.3f exact=%.3f screening=%.3f conditional exact=%.3f\n",
                name,n0,rule,B,m["f1"],m["exact"],mean(R$screen_all),
                if (any(R$screen_all == 1)) mean(R$exact[R$screen_all==1]) else NA_real_))
  }
  raw <- do.call(rbind,audit)
  rownames(raw) <- NULL
  if (save_audit) {
    path <- reviewer_output_dir()
    stem <- sprintf("robustness_%s_%s",name,rule)
    write.csv(raw,file.path(path,paste0(stem,"_replications.csv")),row.names=FALSE)
    write.csv(reviewer_summary(raw),file.path(path,paste0(stem,"_summary.csv")),row.names=FALSE)
    screen <- aggregate(cbind(screen_all,exact) ~ setting+n_y, data=raw, FUN=mean)
    # Conditional exact is undefined if no class survives screening.
    cond <- aggregate(exact_given_screen ~ setting+n_y,data=raw,
                      FUN=function(x) if (all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE),na.action=na.pass)
    screen <- merge(screen, cond, by=c("setting","n_y"), all.x=TRUE, sort=TRUE)
    names(screen)[names(screen)=="exact_given_screen"] <- "conditional_exact"
    write.csv(screen,file.path(path,paste0(stem,"_screening.csv")),row.names=FALSE)
    saveRDS(list(raw=raw,summary=reviewer_summary(raw),screen=screen,seed_rule="b*1000+n0+100000*scenario_index",B=B),
            file.path(path,paste0(stem,"_details.rds")))
  }
  out
}

# ---- explicit CLI: avoid side effects when this script is sourced -----------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)) {
    target <- args[1L]
    Bopt <- grep("^--B=",args,value=TRUE)
    B <- if (length(Bopt)) as.integer(sub("^--B=","",Bopt[1L])) else reviewer_B(500L)
    Ropt <- grep("^--rule=",args,value=TRUE)
    rule <- if (length(Ropt)) sub("^--rule=","",Ropt[1L]) else "rodeo"
    targets <- if (target=="all") setdiff(names(SCEN),"E1V") else target
    res <- lapply(targets,run_scenario,B=B,rule=rule)
    names(res) <- targets
    saveRDS(res,file.path(reviewer_output_dir(),sprintf("robustness_%s_%s_aggregate.rds",target,rule)))
  }
}
