##################################################
## Description: Creates distance matrices
##
## Date: 2018-03-16 15:26:39
## Author: Oskar Hagen (oskar@hagen.bio)
##
## Update Lydian Boschman 2019-07-23
## Update Thomas Keggin   2019-10-18
## 
##################################################

# load libraries ####
lib <- c("raster","matrixStats","sp","gdistance","geosphere","parallel")
sapply(lib,library, character.only=TRUE)

# load data from previous step ####

load("./temp.Rdata")
geoTimes  <- seq(1,length(geoDepthList))
inputData <- geoDepthList

# set variables ####

OutputDir   <- paste0("../../Input/WorldMap200-0Ma_multiple_6d")
crossing_NA <- 0 # Set to 0 (conductance) making land impassible. See gdistance package documentation.

# check, or create, output directories ####

# create dirs if not existing
if (!dir.exists(paste0(OutputDir,"/all_geo_hab"))){
  dir.create(file.path(OutputDir, "all_geo_hab"))
  dir.create(file.path(OutputDir, "geo_dist_m"))
  dir.create(file.path(OutputDir, "geo_dist_m", "geo_dist_m_ti"))
}

# create plot dir if not existing
if (!file.exists(file.path(OutputDir, "/plot"))) {
  dir.create(file.path(OutputDir, "plot"))
}

# create distance matrices ####

t_start <- length(inputData)  # the starting raster index
t_end   <- 1                  # the final raster index (present day)

for (i in t_start:t_end){
  
  Step_1        <- geoDepthList[[i]]
  formattedyear <- geoTimes[i]
  
  #plot
  jpeg(file.path(OutputDir, "plot", paste0(formattedyear,".jpg") ), width = 680, height = 480)
  par(mar=c(0,0,0.2,0.5)+0.2, oma=c(0,0,0,0))
  plot(Step_1, legend.width=1,  legend.shrink=0.64, axes=FALSE, box=FALSE, xlab="", ylab="")
  title(paste("GaSM world @", formattedyear), line=-2.5, cex.main=3)
  dev.off()
  
  hab_180                  <- Step_1        # this is setting up the conductance (cost of dispersal) values for each cell in the raster
  hab_180[is.na(hab_180)]  <- crossing_NA   # this gives the NA valued cells a cost for crossing (land)
  hab_180[!is.na(hab_180)] <- 1             # this gives habitable cells a cost for crossing (1 = no change in cost)
  
  # create a transition object (conductance)
  trans_hab_180 <- transition(hab_180, transitionFunction=min, directions=8) # create matrix with least cost value between each pair of cells (symmetrical?)
  trans_hab_180 <- geoCorrection(trans_hab_180, type = "c", scl = T)         # correct for map distortion
  #go around the world
  points_180 <- rasterToPoints(Step_1)
  points_180_habitable <- points_180[, 1:2] # filter here for habitable cells
  
  # calculate the least-cost distance between points.
  geo_dist_m_ti <- costDistance(trans_hab_180,
                                points_180_habitable,
                                points_180_habitable)
  
  save(geo_dist_m_ti,file=file.path(OutputDir,"geo_dist_m", "geo_dist_m_ti" , paste0("geo_dist_m_ti_t_",i,".RData",sep="")) )
  
  cat("Done with", formattedyear, "\n")
  
}

# create all_geo_hab ####

# merge all temperature rasters into a single dataframe
masterTemp <- as.data.frame(geoTempList[[1]], xy=T)
masterTemp <- masterTemp[,-3]

for(raster in rev(geoTempList)){
  
  raster.df <- as.data.frame(raster, xy=T)
  masterTemp <- cbind(masterTemp,raster.df[,3])
  
}

colnames(masterTemp) <- c("x","y",rev(format(round(geoTimes, 2), nsmall = 2)))

# merge all depth rasters into a single dataframe
masterDepth <- as.data.frame(geoDepthList[[1]], xy=T)
masterDepth <- masterDepth[,-3]

for(raster in rev(geoDepthList)){
  
  raster.df <- as.data.frame(raster, xy=T)
  masterDepth <- cbind(masterDepth,raster.df[,3])
  
}

colnames(masterDepth) <- c("x","y",rev(format(round(geoTimes, 2), nsmall = 2)))

# filter out terrestrial habitat from temperature
masterTemp[is.na(masterDepth)] <- NA

# create and save all_geo_hab object
all_geo_hab <- list(temp = masterTemp, depth = masterDepth)

save(all_geo_hab, file = file.path(OutputDir,"all_geo_hab", paste0("all_geo_hab.RData",sep="")))
# saveRDS()


