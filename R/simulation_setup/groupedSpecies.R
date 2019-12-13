### 100 species, 20 blocks, will need 700 sites to be in a good species/site ratio regime (might be slow to fit?)


## groups of species are related to one another but no other groups
## Will hypothesizes that a low rank model will have a hard time estimating this with only a few latent factors

require(Matrix)

bdiag(Diagonal(2), matrix(1:3, 3,4), diag(3:2))

block = matrix(c(1,0.9,0.9,0.9, 0.9, 
                 0.9, 1, 0.9, 0.9, 0.9,
                 0.9, 0.9, 1, 0.9, 0.9,
                 0.9, 0.9, 0.9, 1, 0.9,
                 0.9, 0.9, 0.9, 0.9, 1), 5,5)

sampleMat = bdiag(block, block, block, block, block,
      block, block, block, block, block,
      block, block, block, block, block,
      block, block, block, block, block
      ) 

simData <- function(mat, numSpecies,  numSites) {
  #browser()
  
  sim_y <- matrix(NA, nrow = numSites, ncol = numSpecies)
  for (i in 1:numSites) {
    samples <- mvrnorm(1, rep(0,numSpecies), mat) ## expects covariance matrix
    sim_y[i,] <- rbinom(numSpecies, size = 1, prob = pnorm(samples))
  }
  
  return(sim_y)
}

numSpecies = 100
numSites = 700

blockData = simData(sampleMat, numSpecies, numSites)

setwd("~/Desktop/communityDimensions")

save(sampleMat, file = "R/simulation_setup/simulation_study_data/matrices/testMats_blockDiag.RData")

save(blockData, file = "R/simulation_setup/simulation_study_data/observed_occurrence/testData_blockDiag.RData")
