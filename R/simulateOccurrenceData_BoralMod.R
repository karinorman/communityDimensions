#### simulate a species covariance ####
library(mvtnorm)
library(Matrix)
library(MASS)

## sparsity: between 0 and 1 (not inclusive)
##           higher number means more sparse

## signal: magnitude of a non-zero entry (stick to stuff between 0 and 1)
## d: size of matrix
## this chunk modified from code originally from Will
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

#### simulate reasonable latent factors given a plausible species covariance matrix ####


## WARNING: I need to run this part by Will and Perry to make sure it makes sense.

set.seed(13211)

numSpecies <- 30
numFactors <- 2

weakCorr <- helperCreateMatrixNice(0.1, 0.25, numSpecies)
## supposed to approximate weak correlations, still not sure if this is doing what I expect

## this chunk modified from code originally from Will
Sigma <- solve(weakCorr)
eig <- eigen(weakCorr)
d <- numSpecies
n <- 1000
Z <- matrix(rnorm(n * d), n)
X <- Z %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors) ## generate data from that covariance matrix

weakCorrF <- factanal(X, factors = numFactors) ## approximate the data with the given covariance matrix with specified number of factors

weakCorrF$loadings ## now treat these as the true latent factors

set.seed(21213)

strongCorr <- helperCreateNice(0.9, 1, numSpecies)
## supposed to approximate strong correlations, still not sure if this is doing what I expect


Sigma <- solve(strongCorr)
eig <- eigen(strongCorr)
Z <- matrix(rnorm(n * d), n)
X <- Z %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)


strongCorrF <- factanal(X, factors = numFactors)
strongCorrF$loadings

#### covariance implied by latent factors - now treated as ground truth ####

weakR <- weakCorrF$loadings %*% t(weakCorrF$loadings) ## covariance implied by latent factors

weakRadj <- nearPD(weakR) ## Not always invertible, so massage
## annoying thing is now you need weakRadj$mat to get the actual matrix

summary(c(weakR) - c(as.matrix(weakRadj$mat))) ## make sure not too far off
summary(c(as.matrix(weakRadj$mat)))

strongR <- strongCorrF$loadings %*% t(strongCorrF$loadings)
strongRadj <- nearPD(strongR)

summary(c(strongR) - c(as.matrix(strongRadj$mat)))
summary(c(as.matrix(strongRadj$mat)))

weakRadjCorr <- cov2cor(weakRadj) ## switch to correlation to remove anything weird that could be going on with the magnitude, might not be necessary
strongRadjCorr <- cov2cor(strongRadj)

weakRadjCorrPrec <- solve(weakRadjCorr) ## precision matrix version if you need it at some point
strongRadjCorrPrec <- solve(strongRadjCorr)


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

test1 <- simData(weakRadj$mat, numCov, numSites, numSpecies)
names(test1)
test2 <- simData(strongRadj$mat, numCov, numSites, numSpecies)