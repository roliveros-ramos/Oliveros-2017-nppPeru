# 
# create database for joint period (csv)

library(ncdf4)
library(kali)
library(lubridate)
library(data.table)
library(stlplus)
library(mgcv)

source("paperPP_functions.R")

# Input file names
seawifsFile = "input/vgpm.s.peru.1997-2009.nc4"
modisFile   = "input/vgpm.m.peru.2002-2017.nc4"
shelfFile   = "input/shelfBreak_peps.csv" # add to package 'peru'?
coastFile   = "input/costa_peps.csv" 

# Output file names
fullDataFile  = "input/fullData.RData"
modelDataFile = "input/modelData.RData"

# open nc files
ncSeawifs = nc_open(seawifsFile)
ncModis   = nc_open(modisFile)
# get variables
npps  = ncvar_get(ncSeawifs, "npp")    # seawifs
time1 = ncvar_get(ncSeawifs, "time")  
nppm  = ncvar_get(ncModis, "npp")      # modis
time2 = ncvar_get(ncModis, "time")
lat   = ncvar_get(ncSeawifs, "lat")    # coordinates
lon   = ncvar_get(ncSeawifs, "lon")
nc_close(ncSeawifs) # close files
nc_close(ncModis)
time = sort(union(time1, time2))

# auxiliar data
DateStamp("Processing auxiliar data")
shelf     = read.csv(shelfFile)
coastLine = read.csv(coastFile)
aux       = expand.grid(lon=lon, lat=lat)
xdist     = getSignedDistance(data=aux, ref=shelf, abs=coastLine)
aux$shelf = round(xdist$dist, 3)
aux$dc    = round(xdist$abs, 3)
aux       = data.table(aux, key = names(aux)[2:1])

# seawifs data
DateStamp("Processing seaWIFS data")
nppsD = stl_array(npps, frequency=12, s.window="periodic") # handle more dimensions
seawifs = expand.grid(lon=lon, lat=lat, time=time1)
seawifs$seawifs = as.numeric(npps)
seawifs$sTrend  = as.numeric(nppsD$trend)
seawifs$sSeason = as.numeric(nppsD$seasonal)
seawifs$sRemain = as.numeric(nppsD$remainder)
seawifs = seawifs[complete.cases(seawifs),]
seawifs = data.table(seawifs, key = c("time", "lat", "lon"))

# modis data
DateStamp("Processing MODIS data")
nppmD = stl_array(nppm, frequency=12, s.window="periodic") 
modis = expand.grid(lon=lon, lat=lat, time=time2)
modis$modis = as.numeric(nppm)
modis$mTrend  = as.numeric(nppmD$trend)
modis$mSeason = as.numeric(nppmD$seasonal)
modis$mRemain = as.numeric(nppmD$remainder)
modis = modis[complete.cases(modis),]
modis = data.table(modis, key = c("time", "lat", "lon"))

# model and prediction (full) data sets
DateStamp("Creating full data frame")
fullGrid = data.table(expand.grid(lon=lon, lat=lat, time=time), 
                      key = c("time", "lat", "lon"))

fullData  = aux[modis[seawifs[fullGrid]], on=c("lat", "lon")]
modelData = aux[modis[seawifs, nomatch=0], on=c("lat", "lon")]

# save data.frames
DateStamp("Saving data...")
save(fullData, file = fullDataFile)
save(modelData, file = modelDataFile)
