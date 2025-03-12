library(ggplot2)
library(sf)
library(sp)
library(spData)
library(rstan)
library(plyr)
library(dplyr)
library(boot)
library(viridis)
library(matrixStats)
library(HDInterval)


rotate <- function(df, degree) {
  dfr <- df
  degree <- pi / (180/degree)
  dfr$longitude = round(df$longitude*cos(degree)-df$latitude*sin(degree), 3)
  dfr$latitude = round(df$longitude*sin(degree)+df$latitude*cos(degree), 3)
  return(dfr)
  
}

### Read In Polygon Coords
df = read.csv("PolygonCoordinates.csv")

### Center on Pacific and scale as we did the grid
df$Longitude = df$Longitude-150
df$Longitude[which(df$Longitude < -175)] = df$Longitude[which(df$Longitude < -175)]+360
scalefactor = 1.5
df$Longitude = df$Longitude/scalefactor
df$Latitude = df$Latitude/scalefactor


### Reading in data and generating grid, taken exactly from model script
DF = read.csv('ComplexityData.csv', stringsAsFactors = T)
DF = DF[which(!is.na(DF$longitude)),]
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
iii <- is.na(as.numeric(st_intersects(pts, world)))
WORLD <- map_data("world", wrap=c(-25,335), ylim=c(-50,75))
points2 = tmp[ii,]
points3 = tmp[iii,]
points2$longitude = points2$longitude-150
points3$longitude = points3$longitude-150
WORLD$long = WORLD$long-150
points2$longitude[which(points2$longitude < -170)] = points2$longitude[which(points2$longitude < -170)]+360
points3$longitude[which(points3$longitude < -170)] = points3$longitude[which(points3$longitude < -170)]+360
WORLD$long[which(WORLD$long < -180)] = WORLD$long[which(WORLD$long < -180)]+360
points3 = points3[points3$longitude != 30,]
latlong_df = rbind(DF[,c('longitude','latitude')],points2)
latlong_df2 = points3


WORLD$long = WORLD$long-25+180
df$Longitude = df$Longitude-25+180
DF$longitude = DF$longitude-25+180
points2$longitude = points2$longitude-25+180

### For labeling later
outcomes = list(
  "Grades",
  "Places",
  "Pitch",
  "Closure",
  "Geminate",
  "TotalConsonants",
  "Qualities",
  "Ancillary",
  "TotalVowels",
  "Diphthongs",
  "TotalTones",
  "MaxOnset",
  "MaxCoda",
  "ContrastiveStress"  
)

### Read In the Model results

fit <- readRDS("fit_FinalDataset.rds")
draws_lambda = extract(fit, pars = c('lambda'))$lambda

str(draws_lambda)


pred_m = matrix(NA, ncol = 14, nrow = nrow(DF))
for(var in 1:14){
  for (n in 1:nrow(DF)) {
    pred_m[n, var] = mean(draws_lambda[1,n,var])
  }
}
pred_m[,14] = boot::inv.logit(pred_m[,14])

DF$ContrastiveStress[DF$ContrastiveStress == levels(DF$ContrastiveStress)[1]] = NA

true_m = as.matrix(data.frame(
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
  CS = as.numeric(as.factor(as.character(DF$ContrastiveStress)))-1
))

true_m_stand = scale(true_m)

draws = extract(fit, pars = c('lambda', 'lambda_pred', 'a', 'a_mu', 'f'))

df$Polygon = gsub('Maritime SEA and Oceania','Southeast Asia and Australasia', df$Polygon)
### World Map
gg = ggplot() + 
  coord_cartesian(xlim = c(min(df$Longitude)+15, max(df$Longitude)-16), ylim = c(min(df$Latitude+6), max(df$Latitude)-6), clip = 'on') +
  geom_polygon(data=WORLD, aes(x=long, y=lat, group=group), colour='grey', fill=NA) +
  geom_path(data = df, aes(x = Longitude, y = Latitude, group = Polygon, color = Polygon), linewidth = 1) +
  geom_polygon(data = df, aes(x = Longitude, y = Latitude, group = Polygon, color = Polygon, fill = Polygon), alpha = .1) +
  geom_point(data = DF, aes(x = longitude, y = latitude), color = 'blue', size = 2, alpha = .7) +
  theme_classic()  +
  theme(panel.background = element_rect(fill = NA, color = "black"), legend.text=element_text(size=14))

pdf(file = 'polygons.pdf', width = 12, height = 6 )
print(gg)
dev.off()

png(file = 'polygons.png', width = 12, height = 6, res = 300, units = 'in' )
print(gg)
dev.off()

##Grid points

gg = ggplot() + 
  coord_cartesian(xlim = c(min(df$Longitude)+15, max(df$Longitude)-16), ylim = c(min(df$Latitude+6), max(df$Latitude)-6), clip = 'on') +
  geom_polygon(data=WORLD, aes(x=long, y=lat, group=group), colour='grey', fill=NA) +
  geom_point(data = points2, aes(x = longitude, y = latitude), shape = 'x', color = 'red', size = 3) +
  geom_point(data = DF, aes(x = longitude, y = latitude), color = 'blue', size = 2, alpha = .7) +
  theme_classic()  +
  theme(panel.background = element_rect(fill = NA, color = "black"))

pdf(file = 'gridpoints.pdf', width = 12, height = 6 )
print(gg)
dev.off()

png(file = 'gridpoints.png', width = 12, height = 6, res = 300, units = 'in' )
print(gg)
dev.off()

### Plotting Everything
columns= c('X', 'Y', 'Outcome') 
smoothDFs = list(
  "North America" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "The Americas" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "American Pacific Rim" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "Eurafrica" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "North Eurasia" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "East Asia" = data.frame(matrix(nrow = 0, ncol = length(columns))),
  "Southeast Asia and Australasia" = data.frame(matrix(nrow = 0, ncol = length(columns)))
)
for (i in 1:7){
  colnames(smoothDFs[[i]]) = columns
}

combDF_big = data.frame()
bigplot_df = data.frame()

for (n in 1:14){
  if(n != 10){
  ### Bind the Smooth Grid Points to the Outcomes
  post_df = data.frame(
    'mean' = numeric(length(draws$lambda_pred[1,,1]))
  )
  for (i in 1:length(draws$lambda_pred[1,,1])){
    post_df$mean[i] = mean(draws$lambda_pred[,i,n])
  }
  post_df$idx = seq(1, nrow(post_df), 1)
  combDF = cbind(points2, post_df)
  if (n==14){
    combDF$mean = inv.logit(combDF$mean)
  }
  
  combDF$outcome = outcomes[[n]]
  
  combDF_big = rbind(combDF_big, combDF)
  ### Construct A list of DFs, one for each region (done this way so points can belong to multiple regions)
  combDF$inPoly = 0
  regions = as.list(unique(df$Polygon))
  plotDFs = list()
  for (i in 1:length(regions)){
    tempdf = df[df$Polygon == regions[i],]
    combDF$inPoly = point.in.polygon(combDF$longitude, combDF$latitude, 
                                     tempdf$Longitude, tempdf$Latitude)
    
    regionDF = combDF[combDF$inPoly != 0,]
    regionDF$region = regions[[i]]
    plotDFs[[i]] = regionDF
    regionDF$Outcome = outcomes[[n]]
    regionDF = regionDF[c('longitude', 'latitude', 'mean', 'Outcome')]
  }
  
  ### Put em all into one big DF, data points in multiple regions will be in once per region 
  plotDF <- ldply(plotDFs, data.frame)
  

  ### Graph each one of the regions, if statements based on which view we are using
  for(x in 1:length(plotDFs)){
    if (plotDFs[[x]]$region[1] == 'North America' | plotDFs[[x]]$region[1] == 'North Eurasia'){
      tmp_df = plotDFs[[x]]
      tmp_df$Xs = tmp_df$longitude
      tmp_df$outcome = outcomes[[n]]
      bigplot_df = rbind(bigplot_df, tmp_df)
      
    }
    if (plotDFs[[x]]$region[1] == 'American Pacific Rim' | plotDFs[[x]]$region[1] == 'East Asia' | plotDFs[[x]]$region[1] == 'The Americas'){
      tmp_df = plotDFs[[x]]
      tmp_df$Xs = tmp_df$latitude
      tmp_df$outcome = outcomes[[n]]
      bigplot_df = rbind(bigplot_df, tmp_df)
    }
    if (plotDFs[[x]]$region[1] == 'Eurafrica'){
      tmp_df = plotDFs[[x]]
      tmp_df$Xs = tmp_df$latitude
      tmp_df$outcome = outcomes[[n]]
      bigplot_df = rbind(bigplot_df, tmp_df)
    }
    if (plotDFs[[x]]$region[1] == 'Southeast Asia and Australasia'){
      rotated_plotDF = rotate(plotDFs[[x]], 45)
      tmp_df = plotDFs[[x]]
      tmp_df$Xs = rotated_plotDF$longitude
      tmp_df$outcome = outcomes[[n]]
      tmp_df$region = "Southeast Asia and Australasia"
      bigplot_df = rbind(bigplot_df, tmp_df)
    }
  }
  }
  

  
}

regions = list("North America", "The Americas", "American Pacific Rim", "Eurafrica", "North Eurasia",
               "East Asia", "Southeast Asia and Australasia")

directions = c("West - East", "South - North", "South - North", "South - North", "West - East",
               "South - North", "Northwest - Southeast")


bigplot_summary = data.frame()
for (i in 1:length(unique(bigplot_df$region))){
  tmp_df = bigplot_df[bigplot_df$region == unique(bigplot_df$region)[i],]
  for (e in 1:length(unique(tmp_df$outcome))) {
    tmp_df2 = tmp_df[tmp_df$outcome == unique(tmp_df$outcome)[e],]
    for (x in unique(tmp_df2$Xs)) {
      tmp_df3 = tmp_df2[tmp_df2$Xs == x,]
      m = median(tmp_df3$mean)
      interval = hdi(tmp_df3$mean)
      l = interval[1]
      u = interval[2]
      bigplot_summary = rbind(bigplot_summary, data.frame(
        'region' = unique(bigplot_df$region)[i],
        'outcome' = unique(tmp_df$outcome)[e],
        'Xs' = x,
        'mean' = m,
        'lower' = l,
        'upper' = u
      ))
    }
  }
}

bigplot_summary$outcome = gsub('ArchiPhonemes', 'TotalConsonants', bigplot_summary$outcome)
bigplot_summary$region = gsub('Maritime SEA and Oceania','Southeast Asia and Australasia', bigplot_summary$region)
bigplot_summary$outcome = factor(bigplot_summary$outcome, levels = c("Grades", "Places", "Closure", "Pitch", "Geminate", "TotalConsonants","Qualities","Ancillary","TotalVowels","TotalTones", "ContrastiveStress", "MaxOnset","MaxCoda"))

regionPlots = list()
for (i in c(1:7)){
  if (i == 4){
    gg = ggplot() +
      geom_linerange(data = bigplot_summary[bigplot_summary$region == "Eurafrica",], aes(x = Xs, y = mean,ymin = lower, ymax = upper), colour="blue",linewidth=1.5, lineend = 'round', alpha = .4)+
      geom_line(data = bigplot_summary[bigplot_summary$region == "Eurafrica",], aes(Xs, mean), colour="black",linewidth=1, lineend = 'round')+
      theme_classic()  + facet_wrap(~outcome, scales = 'free_y') + ggtitle("Eurafrica, South - North")+ xlab('Geographical variable') + ylab('Posterior means') +
      theme(panel.background = element_rect(fill = NA, color = "black")) + geom_vline(data = bigplot_summary[bigplot_summary$region == "Eurafrica",], aes(xintercept = 0), color = 'red', alpha = 0.8)+
      geom_vline(data = bigplot_summary[bigplot_summary$region == "Eurafrica",], aes(xintercept = 25), color = 'red', alpha = 0.8)
  }else{
    
    gg = ggplot() +
      geom_linerange(data = bigplot_summary[bigplot_summary$region == regions[i],], mapping = aes(x = Xs, y = mean, ymin = lower, ymax = upper), colour="blue",linewidth=1.5, lineend = 'round', alpha = .4)+
      geom_line(data = bigplot_summary[bigplot_summary$region == regions[i],], mapping = aes(Xs, mean), colour="black",linewidth=1, lineend = 'round')+
      theme_classic()  + facet_wrap(~outcome, scales = 'free_y') + ggtitle(paste(c(regions[i], ', ', directions[i]), collapse = ''))+
      theme(panel.background = element_rect(fill = NA, color = "black"))
  }
  



regionPlots = list()
for (i in c(1:7)){
  if (i == 4){
    gg = ggplot() +
      geom_point(data = bigplot_df[bigplot_df$region == 'Eurafrica',], aes(Xs, mean), size = 4, alpha = .3, color = 'blue') +
      #stat_smooth(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), method = 'gam', formula = y ~ splines::ns(x, 5), se = FALSE, geom = 'line', linewidth = 1, n = 500, lineend = 'round')+
      stat_summary(data = bigplot_df[bigplot_df$region == 'Eurafrica',], aes(Xs, mean), fun.y=median, colour='black', geom='line',linewidth=1, lineend = 'round')+
      #stat_smooth(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), geom = 'line', color = 'black', alpha = .7, lineend = 'round', linewidth = 1)+
      theme_classic()  + facet_wrap(~outcome, scales = 'free_y') + ggtitle('Eurafrica, South - North')+ xlab('Geographical variable') + ylab('Posterior means') +
      theme(panel.background = element_rect(fill = NA, color = 'black')) + geom_vline(data = bigplot_df[bigplot_df$region == 'Eurafrica',], aes(xintercept = 0), color = 'red', alpha = 0.8)+
      geom_vline(data = bigplot_df[bigplot_df$region == 'Eurafrica',], aes(xintercept = 25), color = 'red', alpha = 0.8)
  }else{
  gg = ggplot() +
    geom_point(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), size = 4, alpha = .3, color = 'blue') +
    #stat_smooth(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), method = 'gam', formula = y ~ splines::ns(x, 5), se = FALSE, geom = 'line', linewidth = 1, n = 500, lineend = 'round')+
    stat_summary(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), fun.y=median, colour='black', geom='line',linewidth=1, lineend = 'round')+ xlab('Geographical variable') + ylab('Posterior means') +
    #stat_smooth(data = bigplot_df[bigplot_df$region == regions[i],], aes(Xs, mean), geom = 'line', color = 'black', alpha = .7, lineend = 'round', linewidth = 1)+
    theme_classic()  + facet_wrap(~outcome, scales = 'free_y') + ggtitle(paste(c(regions[i], ', ', directions[i]), collapse = ''))+
    theme(panel.background = element_rect(fill = NA, color = 'black'))
  }

regionPlots[[i]] = gg
pdf(file = paste(c( 'RegionPlot_', as.character(i), '.pdf'), collapse = ''), width = 12, height = 8 )
print(gg)
dev.off()
png(file = paste(c( 'RegionPlot_png_', as.character(i), '.png'), collapse = ''), width = 12, height = 8, res = 300, units = "in" )
print(gg)
dev.off()
}

pdf(file = 'Contour.pdf', width = 12, height = 8 )
for(e in 1:14){
  if(e != 10){
  gg = ggplot() +
    coord_cartesian(xlim = c(min(combDF_big$longitude), max(combDF_big$longitude)), ylim = c(min(combDF_big$latitude), max(combDF_big$latitude-5)), clip = 'on') +
    geom_raster(data = combDF_big[combDF_big$outcome == outcomes[[e]],], aes(longitude, latitude, fill = mean), alpha = 1, interpolate = TRUE) + 
    scale_fill_viridis() + 
    theme_bw() +
    ggtitle(paste(outcomes[[e]], 'Contour')) +
    geom_polygon(data=WORLD, aes(x=long, y=lat, group=group), colour='black', fill=NA) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position="bottom")
  print(gg)
  }
}

dev.off()

ggplot() + geom_point(data = rotated_plotDF, mapping = aes(x = longitude, y = latitude, color = mean))


diff_df = data.frame()
for (i in c(1:7)){
  
    subdf = bigplot_df[bigplot_df$region == regions[i],]
    subdf2 = aggregate( mean ~Xs + outcome + region ,subdf, median)
    for(e in 1:14){
      if (e != 10){
      subdf3 = subdf2[subdf2$outcome == outcomes[[e]],]
      subdf3 = subdf3[order(subdf3$Xs),]
      subdf3$mean = scale(subdf3$mean, center = F)
      subdf3$d = NA
      subdf3$d[2:nrow(subdf3)] = subdf3$mean[2:nrow(subdf3)]-subdf3$mean[1:(nrow(subdf3)-1)]
      subdf3$d = abs(subdf3$d)
      diff_df = rbind(diff_df, subdf3)
      }
    }
 }

regionDiffs = list()
for (i in c(1:7)){
  gg = ggplot() +
    geom_line(data = diff_df[diff_df$region == regions[i],], aes(Xs, d, colour=outcome),linewidth=1, lineend = 'round')+
    theme_classic()  + ggtitle(regions[i])+
    theme(panel.background = element_rect(fill = NA, color = "black"))
  regionDiffs[[i]] = gg
}
