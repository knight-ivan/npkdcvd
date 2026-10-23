# Standalone regression checks: run from any working directory via Rscript.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])))
source(file.path(script_dir,"rodeo_core.R"))
set.seed(37)
X <- matrix(rnorm(90*7),90,7)
X[,c(1L,5L)] <- matrix(rnorm(90*2,.5,.2),90,2)
fit <- X[1:50,,drop=FALSE]; ev <- X[51:90,,drop=FALSE]
sel <- c(2L,4L,7L)
# Explicit dimension guard detects exactly the old fit/eval-subset mismatch.
stopifnot(inherits(try(rodeo_tstat(ev,fit[,sel,drop=FALSE],rep(1,3)),silent=TRUE),"try-error"))
T <- rodeo_tstat(ev[,sel,drop=FALSE],fit[,sel,drop=FALSE],rep(1,3))
stopifnot(is.matrix(T),identical(dim(T),c(40L,3L)),all(is.finite(T)))
# The whole rodeo must be equivariant to a permutation of the feature columns;
# the original bug changes original coordinate identities after freezing.
perm <- c(4L,7L,1L,6L,2L,5L,3L)
set.seed(43)
a <- rodeo_select_class(X,max_iter=8L,lambda_mult=1.2)
set.seed(43)
b <- rodeo_select_class(X[,perm,drop=FALSE],max_iter=8L,lambda_mult=1.2)
stopifnot(isTRUE(all.equal(a$score,b$score[order(perm)],tolerance=1e-7)),
          identical(a$selected,b$selected[order(perm)]))
# Deterministic randomness given seed and serial execution.
set.seed(81);p <- rodeo_detect_detail(list(X,X),keep=5,mc.cores=1L)
set.seed(81);q <- rodeo_detect_detail(list(X,X),keep=5,mc.cores=1L)
stopifnot(identical(p$detected,q$detected),identical(p$screen,q$screen),
          all(vapply(p$screen,length,integer(1))==5L))
cat("PASS: dimension guard, aligned statistics, column permutation, reproducible screening.\n")
