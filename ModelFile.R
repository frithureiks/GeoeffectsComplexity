library(rstan)
library(bayestestR)
library(sf)
library(spData)
library(plotly)
library(geosphere)
library(bayesplot)
options(mc.cores = parallel::detectCores())

DF = read.csv("ComplexityData.csv", na.strings = "")
DF = DF[which(!is.na(DF$longitude)),]
DF = DF %>% replace(is.na(.), "-1")
DF$ContrastiveStress[DF$ContrastiveStress == 'Yes'] = 1
DF$ContrastiveStress[DF$ContrastiveStress == 'No'] = 0
DF$longitude[DF$longitude < -12] = DF$longitude[DF$longitude < -12] + 360

DF$longitude = DF$longitude - 150
longitude = c()
latitude = c()
step = 5
for(i in seq(-180, 180, step)){
  for(e in seq(-56, 72, step)){
    longitude = c(longitude, i)
    latitude = c(latitude, e)
  }
}
tmp = data.frame(longitude, latitude)
pts <- st_as_sf(tmp, coords=1:2, crs=4326)
ii <- !is.na(as.numeric(st_intersects(pts, world)))
points2 = tmp[ii,]
points2$longitude = points2$longitude-150
points2$longitude[which(points2$longitude < -180)] = points2$longitude[which(points2$longitude < -180)]+360
latlong_df = rbind(DF[,c('longitude','latitude')],points2)
scalefactor = 1.5
latlong_df$longitude = latlong_df$longitude/scalefactor
latlong_df$latitude = latlong_df$latitude/scalefactor
geoDists = matrix(NA, nrow = nrow(latlong_df), ncol = nrow(latlong_df))
for(i in 1:nrow(latlong_df)){
  for (j in 1:nrow(latlong_df)){
    geoDists[i,j] = distHaversine(p1=latlong_df[i,c('longitude','latitude')], p2=latlong_df[j,c('longitude','latitude')],r=6378137*scalefactor)/1000000
  }
}


dat_phylo = list(
  GR = as.integer(DF$Grades),
  PL = as.integer(DF$Places),
  PI = as.integer(DF$Pitch),
  CL = as.integer(DF$Closure),
  GE = as.integer(DF$Geminate),
  AP = as.integer(DF$ArchiPhonemes),
  QU = as.integer(DF$Qualities),
  AN = as.integer(DF$Ancillary),
  TV = as.integer(DF$TotalVowels),
  DI = as.integer(DF$Diphthongs),
  TT = as.integer(DF$TotalTones),
  IS = as.integer(DF$InitialStress),
  FS = as.integer(DF$FinalStress),
  
  CS = as.numeric(DF$ContrastiveStress),
  

  PI_notmis = which(DF$Pitch!= -1),
  GE_notmis = which(DF$Geminate!= -1),
  AP_notmis = which(DF$ArchiPhonemes!= -1),
  DI_notmis = which(DF$Diphthongs!= -1),
  TT_notmis = which(DF$TotalTones!= -1),
  
  CS_notmis = which(DF$ContrastiveStress!= -1),
  

  PI_notmis_len = length(which(DF$Pitch!= -1)),
  GE_notmis_len = length(which(DF$Geminate!= -1)),
  AP_notmis_len = length(which(DF$ArchiPhonemes!= -1)),
  DI_notmis_len = length(which(DF$Diphthongs!= -1)),
  TT_notmis_len = length(which(DF$TotalTones!= -1)),
  
  CS_notmis_len = length(which(DF$ContrastiveStress!= -1)),
  
  PI_mis = which(DF$Pitch == -1),
  GE_mis = which(DF$Geminate == -1),
  AP_mis = which(DF$ArchiPhonemes == -1),
  DI_mis = which(DF$Diphthongs == -1),
  TT_mis = which(DF$TotalTones == -1),
  
  CS_mis = which(DF$ContrastiveStress == -1),
  
  PI_mis_len = length(which(DF$Pitch == -1)),
  GE_mis_len = length(which(DF$Geminate == -1)),
  AP_mis_len = length(which(DF$ArchiPhonemes == -1)),
  DI_mis_len = length(which(DF$Diphthongs == -1)),
  TT_mis_len = length(which(DF$TotalTones == -1)),
  
  CS_mis_len = length(which(DF$ContrastiveStress == -1)),
  
  
  
  means = c(mean(as.numeric(DF$Pitch[DF$Pitch != -1])), mean(as.numeric(DF$Geminate[DF$Geminate != -1])), mean(as.numeric(DF$ArchiPhonemes[DF$ArchiPhonemes != -1])),
      mean(as.numeric(DF$Diphthongs[DF$Diphthongs != -1])),
      mean(as.numeric(DF$TotalTones[DF$TotalTones != -1]))),
  
  sds = c(4*sd(as.numeric(DF$Pitch[DF$Pitch != -1])), 4*sd(as.numeric(DF$Geminate[DF$Geminate != -1])), 4*sd(as.numeric(DF$ArchiPhonemes[DF$ArchiPhonemes != -1])),
            4*sd(as.numeric(DF$Diphthongs[DF$Diphthongs != -1])),
            4*sd(as.numeric(DF$TotalTones[DF$TotalTones != -1]))),
  
  M = as.matrix(geoDists),
  F = as.numeric(as.factor(DF$Family)),
  NFam = max(as.numeric(as.factor(DF$Family))),

  N = nrow(DF),
  N2 = nrow(latlong_df),
  N_pred = nrow(points2),
  NParam = 14,
  ObsID = 1:nrow(DF),
  PredID = (1+nrow(DF)):nrow(latlong_df)
)

saveRDS(dat_phylo, 'dat_phylo_save.RData')

readRDS('dat_phylo_save.RData')

file <- file.path("mat_complex.stan")
fit <- stan(file, data = dat_phylo, iter=4000, chains = 4, refresh=100, cores = 4)
saveRDS(fit, "fit_FinalDataset.rds")
fit <- readRDS("fit_FinalDataset.rds")
print(fit)