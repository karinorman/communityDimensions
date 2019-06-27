## Note: This simulation is based on Johnny's idea that we may not be able to identify correlation between species in the presence of too much zero inflation.


library(MASS)

J = 2
N = 1000


#### almost mutually exclusive ####

set.seed(53133)
x =  mvrnorm(N, rep(0,J), matrix(c(1,-0.99,-0.99,1),byrow=T,nrow=2))

xProb<- pnorm(x)

newY <- matrix(unlist(lapply(c(xProb), function(x) {
  rbinom(1, 1, x)
})), nrow = nrow(xProb), byrow = F)


## zero inflate this

ZI = seq(0,0.9,by=.1)


tryThis = lapply(ZI,function(z){
  matrix(unlist(lapply(c(newY), function(x) {
    x*(runif(1)>z)
  })), nrow = nrow(xProb), byrow = F)
})


estCorr=unlist(lapply(lapply(tryThis,cor),function(x){x[1,2]}))

plot(ZI,estCorr)
## why is the estimated correlation in the presence of no zero inflation not much better?

#### almost perfectly correlated ####

set.seed(53133)
x =  mvrnorm(N, rep(0,J), matrix(c(1,0.99,0.99,1),byrow=T,nrow=2))

xProb<- pnorm(x)

newY <- matrix(unlist(lapply(c(xProb), function(x) {
  rbinom(1, 1, x)
})), nrow = nrow(xProb), byrow = F)


## zero inflate this

ZI = seq(0,0.9,by=.1)


tryThis = lapply(ZI,function(z){
  matrix(unlist(lapply(c(newY), function(x) {
    x*(runif(1)>z)
  })), nrow = nrow(xProb), byrow = F)
})


estCorr=unlist(lapply(lapply(tryThis,cor),function(x){x[1,2]}))

plot(ZI,estCorr)
## why is the estimated correlation in the presence of no zero inflation not much better?


#### no correlation ####

set.seed(53133)
x =  mvrnorm(N, rep(0,J), matrix(c(1,0,0,1),byrow=T,nrow=2))

xProb<- pnorm(x)

newY <- matrix(unlist(lapply(c(xProb), function(x) {
  rbinom(1, 1, x)
})), nrow = nrow(xProb), byrow = F)


## zero inflate this

ZI = seq(0,0.9,by=.1)


tryThis = lapply(ZI,function(z){
  matrix(unlist(lapply(c(newY), function(x) {
    x*(runif(1)>z)
  })), nrow = nrow(xProb), byrow = F)
})


estCorr=unlist(lapply(lapply(tryThis,cor),function(x){x[1,2]}))

plot(ZI,estCorr)

## at least no real danger of getting a false correlation with zero inflated data