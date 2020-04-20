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

setwd("C:/Users/thoma/OneDrive/Documents/Marrey/utilities/input_compilation")

# load libraries ####
lib <- c("raster","matrixStats","sp","gdistance","geosphere","parallel", "dplyr")
sapply(lib,library, character.only=TRUE)

# load data from previous step ####

load("./temp.Rdata")

# set variables ####

OutputDir   <- paste0("../../Input/6d")
crossing_NA <- 0     # Set to 0 (conductance) making land impassible. See gdistance package documentation.
depth_cut   <- -1000 # set the depth cut-off

# check, or create, output directories ####

# create dirs if they don't exist
if (!dir.exists(paste0(OutputDir,"/distances_full"))){
  dir.create(file.path(OutputDir, "distances_full"))
}

# create plot dir if it doesn't exist
if (!file.exists(file.path(OutputDir, "/plot"))) {
  dir.create(file.path(OutputDir, "plot"))
}

# create landscapes ####

# merge all temperature rasters into a single dataframe
masterTemp <- as.data.frame(geoTempList[[1]], xy=T)
masterTemp <- masterTemp[,-3]

for(raster in geoTempList){
  
  raster.df <- as.data.frame(raster, xy=T)
  masterTemp <- cbind(masterTemp,raster.df[,3])
  
}

colnames(masterTemp) <- c("x","y",format(round(geoTimes, 2), nsmall = 2))

# merge all depth rasters into a single dataframe
masterDepth <- as.data.frame(geoDepthList[[1]], xy=T)
masterDepth <- masterDepth[,-3]

for(raster in geoDepthList){
  
  raster.df <- as.data.frame(raster, xy=T)
  masterDepth <- cbind(masterDepth,raster.df[,3])
  
}

colnames(masterDepth) <- c("x","y",format(round(geoTimes, 2), nsmall = 2))

# filter out uninhabitable cells by depth cut off
masterDepth[masterDepth < depth_cut] <- NA
masterTemp[is.na(masterDepth)] <- NA

# explicitly assign rownames
rownames(masterTemp)  <- 1:dim(masterTemp)[1]
rownames(masterDepth) <- 1:dim(masterDepth)[1]

# create and save landscapes object
landscapes <- list(temp = masterTemp, depth = masterDepth)
saveRDS(landscapes, file = file.path(OutputDir, paste0("landscapes.rds",sep="")))


# create distance matrices ####

t_start <- length(geoDepthList)  # the starting raster index
t_end   <- 1                     # the final raster index (present day)

for (i in t_start:t_end){
  
  raster_i <- geoDepthList[[i]]
  crs(raster_i) <- "+proj=longlat +datum=WGS84"
  age      <- geoTimes[i]
  
  conductObj                     <- raster_i      # this is setting up the conductance (cost of dispersal) values for each cell in the raster
  conductObj[!is.na(conductObj)] <- 1             # this gives habitable cells a cost for crossing (1 = no change in cost)
  conductObj[is.na(conductObj)]  <- crossing_NA   # this gives the NA valued cells a cost for crossing (land)
  
  # create a transition object (based on conductance)
  transObj <- transition(conductObj, transitionFunction=min, directions=8) # create matrix with least cost values between each pair of cells (symmetrical?)
  transObj <- geoCorrection(transObj, type = "r", scl = F) * 1000          # correct for map distortion. The output values are in m, the "*1000" converts to km
  # filter by out cells by depth cut off
  df_i            <- as.data.frame(raster_i, xy=TRUE, na.rm = TRUE) # this will remove NA cells
  colnames(df_i)  <- c("x","y","depth")
  df_i_habitable  <- filter(df_i, depth >= depth_cut)[, 1:2]
  mat_i_habitable <- data.matrix(df_i_habitable)
  
  # calculate the least-cost distance between points using the transition object and target (habitable) cells
  dist_mat <- costDistance(transObj,
                                mat_i_habitable,
                                mat_i_habitable)
  
  # number rows and columns
  rownames(dist_mat) <- rownames(masterDepth)[which(!is.na(masterDepth[,i+2]))]
  colnames(dist_mat) <- rownames(masterDepth)[which(!is.na(masterDepth[,i+2]))]
  
  # save the distance matrix
  saveRDS(dist_mat,file=file.path(OutputDir,"distances_full", paste0("distances_full_",i-1,".rds",sep="")) )
  
  # filter out landscapes raster cells by the depth cut off
  # depth
  values(geoDepthList[[i]])[values(geoDepthList[[i]]) < depth_cut] <- NA
  # temp
  geoTempList[[i]][is.na(geoDepthList[[i]][])] <- NA
  
  # plot
  jpeg(file.path(OutputDir, "plot", paste0(round(age, digits = 2),".jpg") ), width = 680, height = 480)
  par(mar=c(0,0,0.2,0.5)+0.2, oma=c(0,0,0,0))
  plot(geoDepthList[[i]], legend.width=1,  legend.shrink=0.64, axes=FALSE, box=FALSE, xlab="", ylab="")
  title(paste("GaSM world @", round(age, digits = 2)), line=-2.5, cex.main=3)
  dev.off()
  
  cat("Done with", round(age, digits = 2), "\n")
  
}



