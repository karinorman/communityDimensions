#### simulate a species covariance ####
library(mvtnorm)
library(Matrix)
library(MASS)

## sparsity: between 0 and 1 (not inclusive)
##           higher number means more sparse

## signal: magnitude of a non-zero entry (stick to stuff between 0 and 1)
## d: size of matrix
## this chunk modified from code originally from Will

## this matrix generation might not be what we want because a matrix can be
## dense but still low rank (well approximated by a small number of latent factors)

## from Will: any matrix can be approximated by a latent factor model with enough latent factors, the question is whether it can be approximated by a few latent factors (like the low rank assumption makes)
helperCreateMatrixNice <- function(sparsity, signal, d) {
  precision.mat <- diag(d)
  precision.mat[upper.tri(precision.mat)] <- signal * (runif(d * (d - 1) / 2) < sparsity)

  precision.mat <- precision.mat + t(precision.mat)

  ## Check Pos-Def
  eig <- eigen(precision.mat)

  if (sum(eig$values > 0) < d) {
    precision.mat <- nearPD(precision.mat) ## adjust if not pos-def
    # eig <- eigen(precision.mat$mat) # can print to see what's happening
    return(precision.mat$mat)
  } else {
    # eig <- eigen(precision.mat)
    return(precision.mat)
  }
}


## instead this function generates matrices with certain number of latent factors
## strength represents how large the covariances can be (larger means stronger relationships)
getMat <- function(numFactors, strength) {
  vectors <- lapply(1:numFactors, function(x, y) {
    runif(30, 0, y)
  }, strength)

  cov <- lapply(vectors, function(x) {
    x %*% t(x)
  })

  A <- matrix(0, nrow = 30, ncol = 30)

  for (i in 1:numFactors) {
    A <- A + cov[[i]]
  }
  return(A)
}
#### simulate reasonable latent factors given a plausible species covariance matrix ####


## WARNING: I need to run this part by Will and Perry to make sure it makes sense.



set.seed(13211)

numSpecies <- 30
numFactors <- 2
strength <- 1

testMat <- getMat(numFactors, strength)


#### given a matrix, simulate occurrence data ####

numCov <- 4
numSites <- 30
numSpecies <- 30

# mat: expects covariance matrix
simData <- function(mat, numCov, numSites, numSpecies) {
  X <- matrix(rnorm(numSites * numCov), numSites, numCov) ## value for every site


  X.coefs <- matrix(rnorm(numSites * numCov), numSites, numCov) ## each species gets a covariate for each of numCov things

  eta <- tcrossprod(as.matrix(X), X.coefs)

  sim_y <- matrix(NA, nrow = numSites, ncol = numSpecies)
  for (j in 1:numSpecies) {
    samples <- mvrnorm(1, eta[, j], mat) ## expects covariance matrix
    sim_y[, j] <- rbinom(numSites, size = 1, prob = pnorm(samples))
  }

  return(list(Y = sim_y, X = X, X.coefs = X.coefs))
}

test1 <- simData(testMat, numCov, numSites, numSpecies)
names(test1)
