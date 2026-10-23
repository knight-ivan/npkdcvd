# R2 Study 5: FINAL adaptive selector (no oracle top-r), with the simulation's
# known shared-background transform as the PRIMARY target-aligned analysis.
# Raw-scale and pooled-empirical-PIT analyses are retained as sensitivities.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))
source(file.path(script_dir,"npkdcvd_core.R"))
source(file.path(script_dir,"rodeo_core.R"))
source(file.path(script_dir,"reviewer_metrics.R"))
C <- 3L; D <- 15L; N_Y <- 200L
truth <- list(c(1L,2L),3L,c(4L,5L))
generate <- function() lapply(seq_len(C),function(y) {
  X <- matrix(rnorm(N_Y*D),N_Y,D)
  if (y==1L) {
    X[,1L] <- rnorm(N_Y,.5,.2); X[,2L] <- rnorm(N_Y,.5,.2)
  } else if (y==2L) X[,3L] <- rnorm(N_Y,3,.2)
  else {X[,4L] <- rnorm(N_Y,.5,.2);X[,5L] <- rnorm(N_Y,.5,.2)}
  X
})
B <- reviewer_B()
reps <- reviewer_mc(seq_len(B),function(b) {
  set.seed(1313L*b+77L)
  X <- generate()

  # PRIMARY: q_j=N(0,1) is known by construction for every shared-background
  # coordinate in Study 5, so apply the exact PIT Q_j=Phi before the rodeo.
  Xexact <- lapply(X,pnorm)
  det_exact <- rodeo_detect(Xexact,keep=D,mc.cores=1L)

  # Sensitivity A: run on the raw scale.
  det_raw <- rodeo_detect(X,keep=D,mc.cores=1L)

  # Sensitivity B: pragmatic pooled empirical PIT when q is not supplied.
  det_emp <- rodeo_detect_bc(X,keep=D,mc.cores=1L)

  one <- function(det,setting) {
    r <- reviewer_metric_rows(det,truth,b,setting,N_Y)
    r$x3_class1 <- det[1L,3L]
    r
  }
  rbind(one(det_exact,"known_background_PIT"),
        one(det_raw,"raw_scale_sensitivity"),
        one(det_emp,"pooled_empirical_PIT_sensitivity"))
})
raw <- do.call(rbind,reps); rownames(raw) <- NULL
summ <- reviewer_summary(raw)
resdir <- reviewer_output_dir()
write.csv(raw,file.path(resdir,"study5_rodeo_replications.csv"),row.names=FALSE)
write.csv(summ,file.path(resdir,"study5_rodeo_summary.csv"),row.names=FALSE)
x3 <- aggregate(x3_class1 ~ replication+setting, data=raw[raw$class==1L,,drop=FALSE],FUN=mean)
x3summary <- reviewer_summary(transform(x3,precision=x3_class1),by="setting",metrics="precision")
x3summary$metric <- "P(X3 selected for class 1)"
write.csv(x3summary,file.path(resdir,"study5_x3_summary.csv"),row.names=FALSE)
saveRDS(list(data=raw,summary=summ,x3=x3summary,B=B,
             primary="known_background_PIT"),file.path(resdir,"study5_rodeo_results.rds"))
print(summ); print(x3summary)
