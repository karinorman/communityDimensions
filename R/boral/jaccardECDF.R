#### make ecdf plot ####
require(parallel)
require(ggplot2)
require(tidyr)
# data is fit.boral
# title is string for ggplot

makeJaccardECDF = function(data, title){
truth1 <- data$yP[nrow(data$yP), ]
parsetemp <- unlist(lapply(strsplit(colnames(data$yStarP), "\\["), function(x) {
  x[2]
}))

siteNum1 <- unlist(lapply(strsplit(parsetemp, ","), function(x) {
  x[1]
}))
speciesNum <- unlist(lapply(strsplit(parsetemp, ","), function(x) {
  x[2]
}))
speciesNum1 <- sub("\\]", "", speciesNum)


test <- cbind.data.frame(t(data$yStarP), speciesNum=speciesNum1, siteNum=siteNum1, truth=truth1)

test2 <- gather(test, -speciesNum, -siteNum, -truth, key = "iteration", value = "val")



jaccardCollapseOverSites <- function(s1, s2) {
  sub1 <- subset(trueDist, speciesNum == s1)
  sub2 <- subset(trueDist, speciesNum == s2)
  num <- which(sub1$truth == 1 & sub2$truth == 1)
  aa <- which(sub1$truth == 1)
  bb <- which(sub2$truth == 1)
  
  jaccard <- length(num) / (length(aa) + length(bb) - length(num))
  coOccur <- length(num)
  return(list(jaccard = jaccard, coOccur = coOccur))
}

trueDist <- subset(test2, iteration == 1)
trueDist$speciesNum <- as.numeric(as.character(trueDist$speciesNum))

pairs <- t(combn(unique(trueDist$speciesNum), 2))

trueJ <- mapply(jaccardCollapseOverSites, pairs[, 1], pairs[, 2], SIMPLIFY = F)

trueJJ <- lapply(trueJ, function(x) {
  x$jaccard
})
trueCO <- lapply(trueJ, function(x) {
  x$coOccur
})



jaccardCollapseOverSitesIteration <- function(s1, s2, Iteration) {
  sub1 <- subset(trueDist, speciesNum == s1 & iteration == Iteration)
  sub2 <- subset(trueDist, speciesNum == s2 & iteration == Iteration)
  num <- which(sub1$val == 1 & sub2$val == 1)
  aa <- which(sub1$val == 1)
  bb <- which(sub2$val == 1)
  jaccard <- length(num) / (length(aa) + length(bb) - length(num))
  coOccur <- length(num)
  return(list(jaccard = jaccard, coOccur = coOccur))
}


p1 <- rep(pairs[, 1], times = 100)
p2 <- rep(pairs[, 2], times = 100)

itN <- rep(901:1000, each = nrow(pairs))

test2$iteration <- as.numeric(as.character(test2$iteration))
trueDist <- subset(test2, iteration > 900)



#ptm <- proc.time()
itJ <- mcmapply(jaccardCollapseOverSitesIteration, p1, p2, itN, SIMPLIFY = F, mc.cores = 4)
#proc.time() - ptm ## 



itN <- rep(901:1000, each = nrow(pairs))
itJJ <- lapply(itJ, function(x) {
  x$jaccard
})
itCO <- lapply(itJ, function(x) {
  x$coOccur
})



byIterationJ <- cbind.data.frame(jaccard = unlist(itJJ), coOccur = unlist(itCO), iteration = itN)

g1=ggplot(byIterationJ,aes(x=jaccard,group=iteration))+stat_ecdf(geom="line",alpha=0.2,col="red")+stat_ecdf(data=data.frame(true = unlist(trueJJ)),aes(x=true, group=1),col="blue",geom="line")+theme_minimal(base_size = 18)+ylab("ECDF") +xlab("jaccard - similarity between species")+ggtitle(title)

return(g1)


}


## time this
## scenarios 6 has most species, will take longest

setwd("~/Desktop/communityDimensions/R/simulation_setup/test_data/correctlySpecifiedResults")

load("r6.RData")
ptm <- proc.time()
makeJaccardECDF(fit.boral, "")
proc.time() - ptm
## 2.5 minutes ish, not bad


load("r1.RData")
ptm <- proc.time()
makeJaccardECDF(fit.boral, "")
proc.time() - ptm ## 2 seconds

## this looks better than 6 which makes sense because species to site ratio

