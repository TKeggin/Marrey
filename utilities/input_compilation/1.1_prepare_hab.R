
# Compilation of GaSM input for use in the marine environment ####
# This creates an all_geo_hab object for gasm input containing sea surface temperature and depth.
# Thomas Keggin
# 09/10/2019

# 00. set the session ####

library("tidyverse")
library("raster")
library("readxl")

# set variables

resolution <- 6

# create blank rasters to aggregate/resample into
resTemplate            <- raster(nrow=180/resolution, ncol=360/resolution, crs=NA)
resTemplateOne         <- raster(nrow=180, ncol=360, crs=NA)
extent(resTemplate)    <- c(-180, 180, -90, 90)
extent(resTemplateOne) <- c(-180, 180, -90, 90)

# function: linear equation function
find.delta <- function(time){
  y <- m*time+b
  y
}

# function: raster to elevation (stolen from Oskar's compilation script)
convert.grey.to.elev <- function(grey_elev_raster){
  elevation <- grey_elev_raster
  values(elevation) <- (values(grey_elev_raster)-155)*40
  values(elevation)[values(grey_elev_raster) > 230] <- 3000 + ((values(grey_elev_raster)[values(grey_elev_raster) > 230]-230)*300)
  values(elevation)[values(grey_elev_raster) < 5] <- -6000 - ((5-values(grey_elev_raster)[values(grey_elev_raster) < 5])*1000)
  return(elevation)
}

# function: glacial scaling
glacial.scale <- function(input,step,factor){
  
 input[[step]]+(lgm*factor)
}

# 01. TEMP load in all koeppen raster files ####

koeppenNames <- list.files("./data/raw_koeppen_temperature_200-0Ma/", pattern = "*.asc")     # names of all the data files

koeppenList <- c()                                # set a blank list in which to put the imported data files

for(i in koeppenNames){                           # import all the data files and store them in dataList
  
  koeppen <- raster(paste("data/raw_koeppen_temperature_200-0Ma/",i, sep = ""))
  koeppenList <- c(koeppenList,koeppen)
  
}

# 02. TEMP replace koeppen band numbers with temperature values ####
temp <- c(26,22,16,5,-10,-20)                            # data from v7
bands <- c(1:6)                                          # Koeppen bands

koeppenTimesteps <- c(1:length(koeppenList))             # make a sequence; 1 to the number of number files, to use in the nested for loop

for(d in koeppenTimesteps){                              # replace Koeppen band numbers with temperature values in all imported rasters
  
  extent(koeppenList[[d]]) <- extent(-180, 180, -90, 90) # fix the extent
  
  for(b in bands){
    
    values(koeppenList[[d]])[values(koeppenList[[d]]) == b] = temp[b]
    
  }
}

# 03. TEMP smooth the interface between koeppen zones ####

koeppenListSmooth <- koeppenList

for(k in koeppenTimesteps){
  
  koeppenListSmooth[[k]] <- focal(koeppenList[[k]], w = matrix(1,3,3),
                                  fun = mean,
                                  pad = TRUE, padValue = -15)  # add this to arbitrarily replace edge NAs - crop the image
  
}

# focal na.rm = TRUE

#rm(koeppenList)

# 04. TEMP create temperature raster for each geo time step ####

deltaTemp.df <- read_excel("./data/Koeppen\ Zonal\ Temperaturesv7.xlsx", range = "A1:B42") # using v7 tropical delta temperature

# categorise geo timesteps into koeppen timestep bins

# load in elevation raster timestep times
elevationTimes <- read_excel("./data/elevation/elevation_times.xlsx")   # adapted from Cris Scotese's "Animation FrameAge_v8 .xls" file (headers changed)
elevationStart <- filter(elevationTimes, cum_frame == 1)               # save the present day raster
elevationTimes <- filter(elevationTimes, age != 0 & age <= 200)        # filter out present day and rasters beyond 200 million years
elevationTimes <- rbind(elevationStart, elevationTimes)

geoPeriod    <- 200 # million years
timeInterval <- elevationTimes$age[length(elevationTimes$age)] - elevationTimes$age[length(elevationTimes$age)-1] # i.e. ultimate timestep - penultimate timestep

geoTimes     <- seq(0, geoPeriod, timeInterval)   # sequence of times in ka from present
geoTimesteps <- seq(1,length(geoTimes))           # serial id for timesteps

koeppenTimes <- seq(0,200000,5000)          # sequence of times in ka from present

deltaTemp <- deltaTemp.df$`Temp Delta`      # vector of delta temp data
deltaBins <- .bincode(geoTimes, koeppenTimes, right = TRUE, include.lowest = TRUE) # assign geo timesteps into delta bins

#rm(deltaTemp.df)

# find the bounding koeppen rasters for each geo time step
lowBand  <- c()
highBand <- c()

for(bin in deltaBins){ 
  
  low  <- bin
  high <- bin+1
  
  lowBand  <- c(lowBand,  low)
  highBand <- c(highBand, high)
  
}

# create the rasters

geoTempList <- c()                     # make a nice list for the rasters to live in

for(step in geoTimesteps){
  
  kLow  <-  lowBand[step]              # select the geo step's low bound
  kHigh <- highBand[step]              # select the goe step's high bound
  
  r1 <- koeppenListSmooth[[kLow]]      # the low bounding raster
  r2 <- koeppenListSmooth[[kHigh]]     # the high bounding raster
  
  t1 <- koeppenTimes[[kLow]]           # the low bounding time
  t2 <- koeppenTimes[[kHigh]]          # the high bounding time
  tg <- geoTimes[step]                 # the time of the geo time step
  
  rd  <- r2-r1                         # the difference between the bounding rasters
  td1 <- t2-t1                         # the difference between the bounding times
  td2 <- tg-t1                         # the difference between the geo time and lower bounding time
  
  p <- td2/td1                         # the proportion of time the geo time step is between the bounding times
  
  vd <- rd*p                           # the raster difference between the low bounding raster and the geo time step raster
  
  rx <- vd+r1                          # calculate the geo time step raster from low bounding raster and the p raster difference
  
  geoTempList <- c(geoTempList, rx)    # put the new raster in the raster house
  
}

#rm(koeppenListSmooth, r1, r2, rd, rx, vd)

# 05. TEMP find the delta temp from present for each geo time step #### this could be done before 04. ####

# calculate linear eq. for each time bin.
deltaVar <- data.frame(c(),c())        # create empty dataframe in which to store equation variables

# calculate linear eq. for each time step
for(bin in unique(deltaBins)){ 
  
  # time
  x1 <- koeppenTimes[bin]
  x2 <- koeppenTimes[bin+1]
  
  # delta temp
  y1 <- deltaTemp[bin]
  y2 <- deltaTemp[bin+1]
  
  # slope
  m <- (y2-y1)/(x2-x1)
  
  # intercept
  b <- y1-m*x1
  
  newVar <- c(m,b)
  deltaVar <- rbind(deltaVar, newVar)
  
}
colnames(deltaVar) <- c("m","b")

#rm(deltaTemp)

# calculate delta temp for each geo time 
geoDeltaTemp <- c() # create a blank vector in which to store delta temp for each timestep

for(i in geoTimesteps){ # the serial number for all the timesteps
  
  
  time <- geoTimes[i]                     # time = the actual time in ka
  deltaBin <- deltaBins[i]                # deltaBin = the assigned bin for that time
  m <- deltaVar[deltaBin,1]               # m = the slope for that bin
  b <- deltaVar[deltaBin,2]               # b = the y intercept for that bin
  delta <- find.delta(time)               # delta = the delta temp calculated for that year
  geoDeltaTemp <- c(geoDeltaTemp,delta)   # add the result to the geoDeltaTemp vector
  
}

#rm(deltaVar)

# 06. TEMP update geo step temperature rasters with delta temp from present ####

for(s in geoTimesteps){
  
  geoTempList[[s]] <- aggregate(geoTempList[[s]], fact=resolution) # set the resolution
  geoTempList[[s]] <- geoTempList[[s]] + geoDeltaTemp[s]
  
}

# 07. TEMP correct temperature for glaciation events ####

lgm            <- raster("./data/glacial/sat.nc")                               # load glacial data
lgmRes         <- res(lgm)[1]                                               # find the extent of glacial data

# some filthy if statements to match the lgm resolution and extent to that of the temperature rasters
if(lgmRes == resolution){                                 # if the resolution matches do nothing
} else {if(lgmRes > resolution){                          
  lgm <- resample(lgm, resTemplate, method = 'bilinear')  # if the lgm res is higher, resample to target res
} else { if(lgmRes < resolution){                         
  lgm <- aggregate(resample(lgm, resTemplateOne, method = 'bilinear'), fact = resolution)}}                   # if the lgm res is lower, aggregate to half the target res (starts at 2)
}

scaleFactor <- read_csv("./data/glacial/frame_table_constant_rate.csv") %>% filter(intensity > 0) # load glacial scaling factors

for(i in seq(1:length(scaleFactor$timestep))){
  
  geoTempList[[scaleFactor$timestep[i]+1]] <-  glacial.scale(geoTempList,
                                                    scaleFactor$timestep[i]+1,
                                                    scaleFactor$intensity[i])
}


# 08. DEPTH load in all elevation raster files ####

geoDepthListRaw <- c()                               # set a blank list in which to put the imported data files
depthExt <- extent(-180, 180, -90, 90)               # set the extent to match temp rasters
depthTemplate <- raster(nrow=180, ncol=360, crs=NA)  # set a blank raster to resample into

depthNames <- list.files("./data/elevation/all/", pattern = "*.tif")     # names of all the data files
  
for(d in depthNames){                           # import all the data files and store them in geoDepthListRAW
  
  depth <- raster(paste("./data/elevation/all/",d, sep = ""))
  extent(depth) <- depthExt
  depth <- resample(depth, depthTemplate, method = "bilinear") # how does this resampling affect the accuracy of the data? # check before and after reef zones
  depth <- aggregate(depth, fact=resolution)                   # set the resolution to target resolution
  geoDepthListRaw <- c(geoDepthListRaw,depth)
  print(paste(era, d))
}

#rm(depth, depthExt, depthTemplate)

# 09. DEPTH interpolate depth to fit geo time steps ####

# The raw elevation data are in rasters (.tiff) that don't follow a consistent time step,
# so using the smallest time step as the constant and interpolating between rasters.

# filter imported rasters by elevationTimes
targetElevations <- c()
for(raster in c(elevationTimes$cum_frame)){
  
  rasteri <- geoDepthListRaw[[raster]]
  targetElevations <- c(targetElevations, rasteri)
}

elevationBins <- .bincode(geoTimes, elevationTimes$age, right = TRUE, include.lowest = TRUE) # assign geo timesteps into raster bins

# find find the bounding elevation rasters for each geo time step

lowBand  <- c()
highBand <- c()

for(bin in elevationBins){ 
  
  low  <- bin
  high <- bin+1
  
  lowBand  <- c(lowBand,  low)
  highBand <- c(highBand, high)
  
}

# create the rasters

geoDepthList <- c()                   # make a nice list for the rasters to live in

for(step in geoTimesteps){
  
  kLow  <-  lowBand[step]             # select the geo step's low bound
  kHigh <- highBand[step]             # select the goe step's high bound
  
  r1 <- targetElevations[[kLow]]      # the low bounding raster
  r2 <- targetElevations[[kHigh]]     # the high bounding raster
  
  t1 <- elevationTimes$age[[kLow]]    # the low bounding time
  t2 <- elevationTimes$age[[kHigh]]   # the high bounding time
  tg <- geoTimes[step]                # the time of the geo time step
  
  rd  <- r2-r1                        # the difference between the bounding rasters
  td1 <- t2-t1                        # the difference between the bounding times
  td2 <- tg-t1                        # the difference between the geo time and lower bounding time
  
  p <- td2/td1                        # the proportion of time the geo time step is between the bounding times
  
  vd <- rd*p                          # the raster difference between the low bounding raster and the geo time step raster
  
  rx <- vd+r1                         # calculate the geo time step raster from low bounding raster and the p raster difference
  
  geoDepthList <- c(geoDepthList, rx) # put the new raster in the raster house
  
}


# 10. DEPTH convert to real elevation values ####

for(t in geoTimesteps){
  
  geoDepthList[[t]] <- convert.grey.to.elev(geoDepthListRaw[[t]])
  
}

# 11. DEPTH filter out terrestrial habitat ####

for(s in geoTimesteps){
  
  values(geoDepthList[[s]])[values(geoDepthList[[s]]) >= 0] <- NA
  
}

# 12. output data to temporary file ####
# want to keep geoTempList and geoDepthList

save(geoTempList, geoDepthList, geoTimes, geoTimesteps, file = "./temp.Rdata")
