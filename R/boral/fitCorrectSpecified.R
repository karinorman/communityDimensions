setwd("~/Desktop/communityDimensions/R/simulation_setup/test_data")

load(file="testData_correctlySpecified.RData")

numSpecies <- c(5, 10, 15, 20, 25, 30)
numFactors <- c(1, 2, 5, 10, 20)

scenarios = expand.grid(numSpecies, numFactors)

require(boral)

for(i in 1:18){ ## 5 latent factors
  ## use correct number of latent variables
      fit.boral <- boral(simulatedData[[i]], family = "binomial", lv.control=list(num.lv=scenarios[i,2]),  row.eff = "none", save.model = FALSE)
      save(fit.boral,file=paste("correctlySpecifiedResults/r",i,".RData", sep=""))
      #mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)

  
}

## boral can't handle more than 5 latent variables, so the rest are going to be misspecified

for(i in 19:nrow(scenarios)){ ## 5 latent factors
  ## use correct number of latent variables
  fit.boral <- boral(simulatedData[[i]], family = "binomial", lv.control=list(num.lv=5),  row.eff = "none", save.model = FALSE)
  save(fit.boral,file=paste("correctlySpecifiedResults/r",i,".RData", sep=""))
  #mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
  
  
}