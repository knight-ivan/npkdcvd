# R2 Study 3: exact-recovery uncertainty + failure-type diagnostics.
# Fold randomness is reproducible inside each explicitly seeded replication.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))
source(file.path(script_dir, "npkdcvd_core.R"))
source(file.path(script_dir, "rodeo_core.R"))
source(file.path(script_dir, "reviewer_metrics.R"))
C <- 5L; D <- 25L; r_y <- rep(3L,C)
starts <- c(1L,4L,7L,10L,13L)
truth <- lapply(seq_len(C),function(y) starts[y]:(starts[y]+2L))
N_GRID <- c(30L,50L,75L,100L,150L,200L,300L,500L)
generate <- function(n) lapply(seq_len(C),function(y) {
  X <- matrix(rnorm(n*D),n,D)
  X[,truth[[y]]] <- matrix(rnorm(n*3L,.5,.5),n,3L)
  X
})
B <- reviewer_B()
by_n <- lapply(N_GRID,function(n) {
  reps <- reviewer_mc(seq_len(B),function(b) {
    set.seed(1000L*b+n)
    d <- rodeo_detect(generate(n),keep=D,mc.cores=1L)
    x <- reviewer_metric_rows(d,truth,b,"study3",n)
    x$failure_type <- ifelse(x$fp==0 & x$fn==0,"exact",
                  ifelse(x$fp>0 & x$fn>0,"mixed",ifelse(x$fp>0,"false_positive","false_negative")))
    x
  })
  do.call(rbind,reps)
})
raw <- do.call(rbind,by_n); rownames(raw) <- NULL
sm <- reviewer_summary(raw)
resdir <- reviewer_output_dir()
write.csv(raw,file.path(resdir,"study3_rodeo_replications.csv"),row.names=FALSE)
write.csv(sm,file.path(resdir,"study3_rodeo_summary.csv"),row.names=FALSE)
failure_counts <- as.data.frame(table(raw$n_y,raw$failure_type),stringsAsFactors=FALSE)
names(failure_counts) <- c("n_y","failure_type","class_replications")
write.csv(failure_counts,file.path(resdir,"study3_rodeo_failures.csv"),row.names=FALSE)
saveRDS(list(raw=raw,summary=sm,failures=failure_counts,B=B),file.path(resdir,"study3_rodeo_results.rds"))
print(sm[sm$metric %in% c("f1","exact"),])
cat("Failure types at n_y=500 (counts are class-replicates; CI uses independent replications):\n")
print(subset(failure_counts,n_y==500))
