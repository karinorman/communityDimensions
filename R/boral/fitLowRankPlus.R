setwd("~/Desktop/communityDimensions/R/simulation_setup/test_data")

load(file="testData_lowRankPlus.RData")

numSpecies = 15
numSites = 100

strength = 1 
signal = 1
numFactors = c(1, 2, 3, 5, 10)
sparsity <- c(0.9, 0.5, 0.2)


mixture <- seq(0, 1, by=.2) ## interpolate between latent factor matrix and perturbation matrix

scenarios = expand.grid(numFactors, sparsity, mixture)


require(boral)

for(i in 1:length(scenarios)){ ## 5 latent factors
  ## use correct number of latent variables
  
  if(scenarios[i,1]>5){
    lf = 5
  }else{
    lf = scenarios[i,1]
  }
  
  fit.boral <- boral(simulatedData[[i]], family = "binomial", lv.control=list(num.lv=lf),  row.eff = "none", save.model = FALSE)
  save(fit.boral,file=paste("lowRankPlus/r",i,".RData", sep=""))
  #mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
  
  
}
