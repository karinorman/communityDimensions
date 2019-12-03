### svd approximations

#### correct spec ####
load( file = "test_data/testMats_correctlySpecified.RData")

numSpecies <- c(5, 10, 15, 20, 25, 30)
numFactors <- c(1, 2, 5, 10, 20)

scenarios = expand.grid(numSpecies, numFactors)

test=svd(trueMats[[1]],nu = 1,nv = 1)

test$u %*% test$d[1] %*% t(test$v)

trueMats[[1]]


test=svd(trueMats[[7]],nu = 2,nv = 2)

test$u %*% (diag(2)*test$d[1:2]) %*% t(test$v)

trueMats[[7]]


## get a lower rank approximation to a higher rank matrix
svdApprox <- function(mat, numFactors){
#browser()
  test=svd(mat,nu = numFactors,nv = numFactors)
  
  test$u %*% (diag(numFactors)*test$d[1:numFactors]) %*% t(test$v)
}

#### low rank plus ####


load(file="test_data/testLFMats_lowRankPlus.RData")
load(file="test_data/testPerturbMats_lowRankPlus.RData")


numSpecies = 15
numSites = 100

strength = 1 
signal = 1
numFactors = c(1, 2, 3, 5, 10)
sparsity <- c(0.9, 0.5, 0.2)


mixture <- seq(0, 1, by=.2) ## interpolate between latent factor matrix and perturbation matrix

scenarios = expand.grid(numFactors = numFactors, sparsity = sparsity, mixture = mixture)

#lowRankApprox = mapply(svdApprox, trueMats, scenarios[,1], SIMPLIFY = F)

lowRankApprox = vector("list",length(trueMats))
for(i in 1:length(lfM)){
  mixture= scenarios[i,3]
  mat = lfM[[i]]
  addM = perturbM[[i]]
  lowRankApprox[[i]]=svdApprox(mixture*mat+(1-mixture)*solve(as.matrix(addM)), scenarios[i,1])
  print(i)
}


source("metrics.R")

kl <- c()
frob <- c()
for(i in 1:length(lfM)){
  mixture= scenarios[i,3]
  mat = lfM[[i]]
  addM = perturbM[[i]]
  kl <- c(kl, klDivergence(lowRankApprox[[i]]+diag(nrow(lfM[[i]])), mixture*mat+(1-mixture)*solve(as.matrix(addM)) +diag(nrow(lfM[[i]])) ))
  frob <- c(frob, frobeniusNorm(lowRankApprox[[i]], mixture*mat+(1-mixture)*solve(as.matrix(addM))))
}

results = cbind.data.frame(scenarios, kl, frob)


ggplot(subset(results, numFactors==1), aes(mixture, kl,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==1), aes(mixture, frob,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()


ggplot(subset(results, numFactors==2), aes(mixture, kl,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==2), aes(mixture, frob,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==5), aes(mixture, kl,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==5), aes(mixture, frob,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==10), aes(mixture, kl,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

ggplot(subset(results, numFactors==10), aes(mixture, frob,col=as.factor(numFactors)))+geom_point(cex=2)+geom_line()+facet_wrap(~sparsity)+theme_minimal()

