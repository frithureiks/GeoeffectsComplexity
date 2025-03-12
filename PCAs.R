library(rstan)
library(bayestestR)
library(sf)
library(spData)
library(plotly)
library(geosphere)
library(bayesplot)
library(factoextra)
library(blavaan)

DF = read.csv("ComplexityData.csv", na.strings = "")
DF = DF[which(!is.na(DF$longitude)),]
DF$longitude[DF$longitude < -12] = DF$longitude[DF$longitude < -12] + 360

calc_DF = DF[,c(6:19)]
calc_DF$ContrastiveStress = as.numeric(as.factor(calc_DF$ContrastiveStress))-1
for (col in 1:ncol(calc_DF)) {
  calc_DF[,col] = scale(calc_DF[,col], center = F)
}


model <- paste(c('latent =~', paste(colnames(calc_DF), collapse = '+')), collapse = '')
fit <- bsem(model, data = DF, save.lvs=TRUE, mcmcfile = T)
summary(fit)


