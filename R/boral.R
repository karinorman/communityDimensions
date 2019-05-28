detach("package:boral",unload=T) ## get rid of old boral if you have it

wd = "" ## where you have bryophytes stored
setwd(wd)
load(file="bryophytes.rda")

J <- 30 ## downsample to make things faster
N <- 150

set.seed(322223)
toKeep = sample(1:ncol(bryophytes$Y),J)
toKeepN = sample(1:nrow(bryophytes$Y),N)
newY = bryophytes$Y[,toKeep]
newY = newY[toKeepN,] ## occur/no occur
newX= bryophytes$X[toKeepN,] ## covariates, first column is an intercept so we remove this later on
## boral adds an intercept automatically

library(devtools)
install_github("sastoudt/boral",  ref = "ppc") ## my version that gives you posterior predictive information
require(boral)

#ptm <- proc.time()
mod <- boral(newY,newX[,-1],family="binomial",lv.control=list(num.lv=2),save.model = TRUE,row.eff = "fixed")
proc.time() - ptm  ## 919.882 15 minutes ish

## fixed effect per site
## num.lv: number of latent factors

## posterior predictive info only works for binomial (presence-absence) right now
## but if you need counts, let me know and I can add posterior predictive checks for that too

names(mod)

## I added yP, yStarP, ZP, ZStarP
