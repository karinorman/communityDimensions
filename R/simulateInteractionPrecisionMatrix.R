## adaptation of Will's simulation code
library(Matrix) ## nearPD to get nearest positive definite matrix (so we can invert it)


J <- 20 ## species
N <- 1000 ## sites


sparsity <- c(.1, .25, .5, .75, .9, 1)
## how often we allow a non-zero entry in covariate matrix
## larger numbers mean less sparse (more non-zero entries)

## signal: how big the non-zero entry is

## generate a precision matrix
makePrecisionMat <- function(sparsity, signal, J) {
  precision.mat <- diag(J)
  precision.mat[upper.tri(precision.mat)] <- signal * (runif(J * (J - 1) / 2) < sparsity)
  precision.mat <- precision.mat + t(precision.mat)

  ## Check Pos-Def
  eig <- eigen(precision.mat)
  # print(sum(eig$values>0))

  if (sum(eig$values > 0) < J) {
    precision.mat <- nearPD(precision.mat) ## adjust if need be
    eig <- eigen(precision.mat$mat)
    # print(sum(eig$values>0))
    return(precision.mat$mat)
  } else {
    eig <- eigen(precision.mat)
    # print(sum(eig$values>0))
    return(precision.mat)
  }
}

## given a signal, change sparsity
testMat <- lapply(sparsity, makePrecisionMat, .05, J)
testMat <- lapply(sparsity, makePrecisionMat, .2, J)
testMat <- lapply(sparsity, makePrecisionMat, .5, J)

scenarios <- expand.grid(sparsity, c(0.05, 0.2, 0.5))

## change signal and sparsity at the same time
testMat <- mapply(makePrecisionMat, scenarios[, 1], scenarios[, 2], MoreArgs = list(J = J), SIMPLIFY = F)


## Note: when using these to simulate datasets, remember that this is the precision matrix not the covariance matrix (need to invert it)