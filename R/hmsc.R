#### new version of HMSC ####

# https://github.com/hmsc-r/HMSC
library(devtools)
install_url("https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz")
install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))

library(Hmsc)
load(here::here("data", "bryophytes.rda"))

J=30

set.seed(322223)
toKeep <- sample(1:ncol(bryophytes$Y), J)



newY <- bryophytes$Y[, toKeep]
newTr <- bryophytes$Tr[, toKeep]
newPhylo <- bryophytes$Phylo[toKeep, toKeep]


studyDesign <- data.frame(tree = bryophytes$Random[, 1], siteT = bryophytes$Random[, 2])
rL1 <- HmscRandomLevel(units = studyDesign$tree)
rL2 <- HmscRandomLevel(units = studyDesign$siteT)
mod <- Hmsc(Y = newY, X = bryophytes$X, XScale = TRUE, YScale = FALSE, Tr = t(newTr), TrScale = FALSE, C = newPhylo, studyDesign = studyDesign, ranLevels = list("tree" = rL1, "siteT" = rL2))

mod <- sampleMcmc(mod, samples = 2000, thin = 1, transient = 1000) ## transient is burnin

## Note: I briefly compared the results of this with the results of the code provided for the Ovaskainen et al '17 Ecology Letters (which used an original implementation). The predictions from the two models had visually similar densities. However, please double check that this model specification makes sense.

mpost = convertToCodaObject(mod)
names(mpost)

## these outputs are gnarly but the names match the hierarchy schematic in Figure 4 of Ovaskainen et al '17 Ecology Letters.

library(dplyr)

summary(mpost$Beta)$statistics %>% dim() ## 5 covariates, each species gets a set
summary(mpost$Beta)$statistics %>% head()

summary(mpost$Gamma)$statistics %>% dim() ## trait by covariate


summary(mpost$V)$statistics %>% dim() ## trait by covariate

summary(mpost$Sigma)$statistics %>% dim() ## per species

mpost$Eta[[1]][[1]] %>% dim() ## four factors, by site
mpost$Lambda[[1]][[1]] %>% dim() ## four factors, by species
mpost$Omega[[1]][[1]] %>% dim() ## species by species

