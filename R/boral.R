detach("package:boral",unload=T) ## get rid of old boral if you have it
library(here)

#### data ####
load(here::here("data", "bryophytes.rda"))
## source file originally from: https://github.com/guiblanchet/HMSC/blob/master/data/bryophytes.rda

names(bryophytes)
## Y: binary, occurrence
dim(bryophytes$Y) ## 204 (sites) x 60 (species)

## X: covariates
dim(bryophytes$X) ## 204 (sites) x 5 (number of covariates)
colnames(bryophytes$X)

## Tr: traits
dim(bryophytes$Tr) ## 5 (traits) x 60 (species)
row.names(bryophytes$Tr)

## Phylo: phylogenetic relationships between species
dim(bryophytes$Phylo) ## 60 (species) x 60 (species)

## Random: structure for random effects
dim(bryophytes$Random) ## 204 (sites) x 2
length(unique(bryophytes$Random[, 1])) ## 204 (trees)
length(unique(bryophytes$Random[, 2])) ## 28 (sites)


#### subsetting ####
J <- 30 ## downsample to make things faster
N <- 150

set.seed(322223)
toKeep = sample(1:ncol(bryophytes$Y),J)
toKeepN = sample(1:nrow(bryophytes$Y),N)
newY = bryophytes$Y[,toKeep]
newY = newY[toKeepN,] ## occur/no occur
newX= bryophytes$X[toKeepN,] ## covariates, first column is an intercept so we remove this later on
## boral adds an intercept automatically

#### model ####
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
