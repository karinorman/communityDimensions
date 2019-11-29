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
simData <- function(mat, numSites, numSpecies) {

  
  sim_y <- matrix(NA, nrow = numSites, ncol = numSpecies)
  for (i in 1:numSites) {
    samples <- mvrnorm(1, rep(0,numSpecies), mat) ## expects covariance matrix
    sim_y[i,] <- rbinom(numSites, size = 1, prob = pnorm(samples))
  }
  
  return(sim_y)
}

test1 <- simData(testMat,  numSites, numSpecies)
test1


## keep numSites and signal fixed

numSpecies <- c(5, 10, 15, 20, 25, 30)
numFactors <- c(1, 2, 5, 10, 20)

## then do a more realistic example

numSpecies = 50
numSites = 200
numFactors <- c(1, 2, 5, 10, 20)

