#### helper function ####
makePrecisionMat <- function(sparsity, signal, J) {
  precision.mat <- diag(J)
  precision.mat[upper.tri(precision.mat)] <- signal * (runif(J * (J - 1) / 2) < sparsity)
  precision.mat <- precision.mat + t(precision.mat)

  ## Check Pos-Def
  eig <- eigen(precision.mat)
  # print(sum(eig$values>0))

  if (sum(eig$values > 0) < J) {
    precision.mat <- nearPD(precision.mat) ## adjust if need be
    eig <- eigen(precision.mat$mat)
    # print(sum(eig$values>0))
    return(precision.mat$mat)
  } else {
    eig <- eigen(precision.mat)
    # print(sum(eig$values>0))
    return(precision.mat)
  }
}


J <- 5 ## species
N <- 50 ## sites

library(Matrix)
library(MASS)


#### with sparse interaction covariance ####
set.seed(4419)

cov <- as.matrix(runif(N, -5, 5), ncol = 1) ## value of covariate by site

beta0 <- as.matrix(c(-10, -2, 6, 8, 10), nrow = 1) ## intercept, drives prevalence
beta1 <- as.matrix(c(2, 1, 15, 4, 2), nrow = 1) ## slope


muLinear <- mapply(function(x, y) {
  x + y * cov
}, beta0, beta1)

speciesPrec <- makePrecisionMat(0.2, 0.5, J)

samp <- do.call("rbind", lapply(1:N, function(x) {
  mvrnorm(1, muLinear[x, ], solve(speciesPrec))
}))

sampN <- pnorm(samp)

newY <- matrix(unlist(lapply(c(sampN), function(x) {
  rbinom(1, 1, x)
})), nrow = nrow(sampN), byrow = F)

newX <- cov

mod <- boral(newY, newX, family = "binomial", lv.control = list(num.lv = 2), save.model = TRUE)


#### with induced latent factor due to misspecification ####
set.seed(4419)

cov = as.matrix(runif(N,-5,5),ncol=1) ## value of covariate by site
cov2= as.matrix(runif(N,-5,5),ncol=1)

beta0=as.matrix(c(-10, -2, 6, 8, 10),nrow=1) ## intercept, drives prevalence
beta1=as.matrix(c(2,1,15,4,2),nrow=1) ## slope
beta2=as.matrix(rep(5,J),nrow=1) ## covariate that we will leave out in the analysis
## This will be like having a latent factor

muLinear = mapply(function(x,y,z){x+y*cov+z*cov2},beta0,beta1,beta2)

speciesPrec=diag(J) ## no correlation between species except for the one induced by the missing covariate

samp = do.call("rbind",lapply(1:N,function(x){mvrnorm(1,muLinear[x,] ,solve(speciesPrec))}))

sampN=pnorm(samp)

newY =matrix(unlist(lapply(c(sampN),function(x){rbinom(1,1,x)})),nrow=nrow(sampN),byrow = F)

newX = cov

mod <- boral(newY,newX,family="binomial",lv.control=list(num.lv=2),save.model = TRUE)