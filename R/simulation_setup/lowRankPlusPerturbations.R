require(Matrix)
library(mvtnorm)
library(Matrix) 
library(MASS)
require(boral)
library(dplyr)

## generates matrices with certain number of latent factors
## strength represents how large the covariances can be (larger means stronger relationships)
getMat <- function(numFactors, numSpecies, strength) {
  vectors <- lapply(1:numFactors, function(x, y) {
    runif(numSpecies, 0, y)
  }, strength)
  
  cov <- lapply(vectors, function(x) {
    x %*% t(x)
  })
  
  A <- matrix(0, nrow = numSpecies, ncol = numSpecies)
  
  for (i in 1:numFactors) {
    A <- A + cov[[i]]
  }
  return(A)
}

sparsity <- c(0.9, 0.5, 0.2)
signal <- 1

## sparsity: between 0 and 1 (not inclusive)
##           higher number means more sparse

## signal: magnitude of a non-zero entry (stick to stuff between 0 and 1)
## d: size of matrix (species)
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

## generate simulated given latent factor covariance and the perturbation matrix
simData <- function(mat, addM, mixture, numSpecies, numSites){
  #browser()
  sim_y=matrix(NA,nrow=numSites, ncol=numSpecies)
  for(i in 1:numSites){
    samples = mvrnorm(1,rep(0,numSpecies),mixture*mat+(1-mixture)*solve(as.matrix(addM))) 
    sim_y[i,] <- rbinom(numSpecies, size = 1, prob = pnorm(samples)) 
  }
  
  return(sim_y)
  
}

## change to the ratio that does best in correctly specified case
numSpecies = 30
numSites = 30

strength = 1 
signal = 1
numFactors = c(1, 2, 3, 5, 10, 20)
sparsity <- c(0.9, 0.5, 0.2)

mixture <- seq(0, 1, by=.1) ## interpolate between latent factor matrix and perturbation matrix

lfM <- getMat(numFactors[1],strength)
perturbM <- helperCreateMatrixNice(sparsity[1],signal, numSpecies)

test = simData(lfM, perturbM, mixture[5], numSpecies, numSites)