library(mvtnorm)
library(Matrix)
library(MASS)

## generates matrices with certain number of latent factors
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

set.seed(13211)

numFactors <- 2
strength <- 1

testMat <- getMat(numFactors, strength)


#### given a matrix, simulate occurrence data ####

numSites <- 30
numSpecies <- 30

# mat: expects covariance matrix
# no fixed effects
simData <- function(mat, signal, numSites, numSpecies) {
  X <- matrix(rnorm(numSites), numSites, 1) ## value for every site
  
  X.coef <- signal
  X.shift <- matrix(rnorm(numSites), numSites, 1) ## each species gets a shift from same effect, this will introduce nonlinearity in the covariance between species
  eta <- tcrossprod(as.matrix(X), X.coef)
  
  
  sim_y <- matrix(NA, nrow = numSites, ncol = numSpecies)
  for (i in 1:numSites) {
    samples <- mvrnorm(1, eta[i]+X.shift, mat) ## expects covariance matrix
    sim_y[i,] <- rbinom(numSites, size = 1, prob = pnorm(samples))
  }
  
  return(sim_y)
}

signal = 1

test1 <- simData(testMat,  signal, numSites, numSpecies)
test1


## thing that introduces more/less nonlinearity: how signal to shift magnitude differs

## use good species to site ratio from correctly specified

signal = c(0.1, 0.5, 1, 2, 5)
numFactors <- c(1, 2, 5, 10, 20)

