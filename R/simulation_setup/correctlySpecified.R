library(mvtnorm)
library(Matrix)
library(MASS)

## generates matrices with certain number of latent factors
## strength represents how large the covariances can be (larger means stronger relationships)
getMat <- function(numFactors, numSpecies, strength) {
  vectors <- lapply(1:numFactors, function(x, y) {
    runif(numSpecies, -y, y) ##  could go to rnorm() too
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

set.seed(13211)

numFactors <- 2
strength <- 1
numSpecies <- 30

testMat <- getMat(numFactors, numSpecies, strength)


#### given a matrix, simulate occurrence data ####

numSites <- 30
numSpecies <- 30

# mat: expects covariance matrix
# no fixed effects
simData <- function(mat, numSpecies,  numSites) {
#browser()
  
  sim_y <- matrix(NA, nrow = numSites, ncol = numSpecies)
  for (i in 1:numSites) {
    samples <- mvrnorm(1, rep(0,numSpecies), mat) ## expects covariance matrix
    sim_y[i,] <- rbinom(numSpecies, size = 1, prob = pnorm(samples))
  }
  
  return(sim_y)
}

test1 <- simData(testMat,  numSpecies, numSites)
test1


## keep numSites and strength fixed

numSpecies <- c(5, 10, 15, 20, 25, 30)
numFactors <- c(1, 2, 5, 10, 20)

## num sites = 30
scenarios = expand.grid(numSpecies=numSpecies, numFactors=numFactors) ## 30

setwd("~/Desktop/communityDimensions")
write.csv(scenarios,"R/simulation_setup/simulation_study_data/correctSpecificationScenarios.csv",row.names=F)

trueMats = mapply(getMat, scenarios[,2], scenarios[,1], strength, SIMPLIFY = F)

simulatedData = mapply(simData, trueMats, scenarios[,1], numSites, SIMPLIFY = F)

save(trueMats, file = "R/simulation_setup/simulation_study_data/matrices/testMats_correctlySpecified.RData")

save(simulatedData, file = "R/simulation_setup/simulation_study_data/observed_occurrence/testData_correctlySpecified.RData")

## then do a more realistic example

numSpecies = 50
numSites = 200
numFactors <- c(1, 2, 5, 10, 20)

scenarios = expand.grid(numSpecies=numSpecies, numFactors=numFactors) ## 5
write.csv(scenarios,"R/simulation_setup/simulation_study_data/correctSpecificationScenariosRealistic.csv",row.names=F)

trueMats = mapply(getMat, scenarios[,2], scenarios[,1], strength, SIMPLIFY = F)

simulatedData = mapply(simData, trueMats, scenarios[,1], numSites, SIMPLIFY = F)

save(trueMats, file = "R/simulation_setup/simulation_study_data/matrices/testMats_correctlySpecifiedRealistic.RData")

save(simulatedData, file = "R/simulation_setup/simulation_study_data/observed_occurrence/testData_correctlySpecifiedRealistic.RData")
