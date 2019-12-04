setwd("~/Desktop/communityDimensions/R/simulation_setup/test_data")
source("../metrics.R")


load(file="testMats_correctlySpecified.RData")

numSpecies <- c(5, 10, 15, 20, 25, 30)
numFactors <- c(1, 2, 5, 10, 20)

scenarios = expand.grid(numSpecies, numFactors)

kl <- c()
frob <- c()
for(i in 1:18){
  

load(file=paste0("correctlySpecifiedResults/r",i,".RData"))


estimatedMat = fit.boral$lv.coefs.mean[,-1]%*%t(fit.boral$lv.coefs.mean[,-1])

## put lambda lambda together


## kl divergence

kl <- c(kl,klDivergence(trueMats[[i]]+diag(nrow(trueMats[[i]])), estimatedMat+diag(nrow(trueMats[[i]]))))

## frobenius

frob <- c(frob,frobeniusNorm(estimatedMat, trueMats[[i]]))

}

results = cbind.data.frame(scenarios[1:18,], kl, frob)
names(results)[1:2]= c("numSpecies", "numFactors")

require(ggplot2)

ggplot(results, aes(numSpecies, kl,col=as.factor(numFactors)))+geom_point()
ggplot(results, aes(numSpecies, frob,col=as.factor(numFactors)))+geom_point()


### misspecified ####

kl <- c()
frob <- c()
for(i in 19:nrow(scenarios)){
  
  
  load(file=paste0("correctlySpecifiedResults/r",i,".RData"))
  
  
  estimatedMat = fit.boral$lv.coefs.mean[,-1]%*%t(fit.boral$lv.coefs.mean[,-1])
  
  ## put lambda lambda together
  
  
  ## kl divergence
  
  kl <- c(kl,klDivergence(trueMats[[i]]+diag(nrow(trueMats[[i]])), estimatedMat+diag(nrow(trueMats[[i]]))))
  
  ## frobenius
  
  frob <- c(frob,frobeniusNorm(estimatedMat, trueMats[[i]]))
  
}

results = cbind.data.frame(scenarios[19:nrow(scenarios),], kl, frob)
names(results)[1:2]= c("numSpecies", "numFactors")

require(ggplot2)

ggplot(results, aes(numSpecies, kl,col=as.factor(numFactors)))+geom_point()
ggplot(results, aes(numSpecies, frob,col=as.factor(numFactors)))+geom_point()

## compare magnitudes to other ones

## this needs to go in markdown on github
