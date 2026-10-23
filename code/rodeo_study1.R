# R2 Study 1: FINAL cross-fitted rodeo (adaptive largest gap; no oracle r_y).
# Retains identical original DGP as V1_jaccard.R and the old rodeo_study1.R.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))
source(file.path(script_dir, "npkdcvd_core.R"))
source(file.path(script_dir, "rodeo_core.R"))
source(file.path(script_dir, "reviewer_metrics.R"))
C <- 10L; D <- 30L; N_Y <- 150L
truth <- lapply(seq_len(C), function(y) (((y-1L)+0:5) %% D)+1L)
jacc <- function(a,b) {
  u <- length(union(a,b))
  if (u == 0L) return(1)  # J(empty,empty) conventional convention
  length(intersect(a,b))/u
}
jacmat <- function(sets) outer(seq_len(C),seq_len(C),Vectorize(function(y,yp) jacc(sets[[y]],sets[[yp]])))
J0 <- jacmat(truth)
generate <- function() lapply(seq_len(C),function(y) {
  X <- matrix(runif(N_Y*D),N_Y,D)
  for (k in seq_along(truth[[y]])) X[, truth[[y]][k]] <- rnorm(N_Y,.5,.02*k)
  X
})
B <- reviewer_B(); out <- reviewer_mc(seq_len(B), function(b) {
  set.seed(777L*b+13L)
  det <- rodeo_detect(generate(),keep=D,mc.cores=1L)
  sets <- lapply(seq_len(C),function(y) which(det[y,]!=0L))
  list(J=jacmat(sets),rows=reviewer_metric_rows(det,truth,b,"study1",N_Y),
       all_exact=as.integer(all(vapply(seq_len(C),function(y) setequal(sets[[y]],truth[[y]]),logical(1)))))
})
Jbar <- Reduce("+",lapply(out,`[[`,"J"))/B
raw <- do.call(rbind,lapply(out,`[[`,"rows")); rownames(raw) <- NULL
summ <- reviewer_summary(raw)
per_rep <- data.frame(replication=seq_len(B),
  all_class_exact=vapply(out,`[[`,integer(1),"all_exact"),
  max_abs_error=vapply(out,function(x) max(abs(x$J-J0)),numeric(1)))
resdir <- reviewer_output_dir()
write.csv(raw,file.path(resdir,"study1_rodeo_replications.csv"),row.names=FALSE)
write.csv(summ,file.path(resdir,"study1_rodeo_summary.csv"),row.names=FALSE)
write.csv(per_rep,file.path(resdir,"study1_rodeo_overall.csv"),row.names=FALSE)
saveRDS(list(true_jac=J0,est_jac=Jbar,overall=per_rep,summary=summ,replications=raw,
             B=B,selector="cross-fitted rodeo; deterministic largest gap"),
        file.path(resdir,"study1_rodeo_results.rds"))
cat(sprintf("Study 1 formal rodeo: B=%d | all-class exact=%.4f | max|mean J-J0|=%.6g\n",
 B,mean(per_rep$all_class_exact),max(abs(Jbar-J0))))
# The paper figure is drawn from study1_rodeo_results.rds by study1_jaccard_figure.R.
