library(ncdf4)
library(kali)
library(lubridate)
library(data.table)
library(stlplus)
library(mgcv)

source("paperPP_functions.R")

# File names
seawifsFile = "output/vgpm.s.peru.1997-2009.nc4"
modisFile   = "output/vgpm.m.peru.2002-2016.nc4"

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

# Problem: difference between seaWIFS and MODIS, unified time series for analysis

nppsD = stl_array(npps) # handle more dimensions
nppmD = stl_array(nppm) # handle more dimensions

npp = array(dim=c(dim(npps)[1:2], length(jtime), 2))
npp[,,,1] = nppm[,, time2 %in% time1]
npp[,,,2] = npps[,, time1 %in% time2]

# lat, lon
dnpp  = npp[, , , 1] - npp[, , , 2]
mnpp  = 0.5*(npp[, , , 1] + npp[, , , 2])
drnpp = dnpp/mnpp
rnpp  = apply(npp, c(1,2), 
              FUN=function(x) cor(x[,1], x[,2], use="pair"))


# lat, dc
par(mfrow=c(1,4), mar=c(3,3,1,3))
image.map(lon, lat, 1e-3*apply(npps, 1:2, mean, na.rm=TRUE),
          zlim=c(0, 7))
image.map(lon, lat, 1e-3*apply(nppm, 1:2, mean, na.rm=TRUE),
          zlim=c(0, 7))
image.map(lon, lat, apply(drnpp, 1:2, mean, na.rm=TRUE))
image.map(lon, lat, rnpp)

zm = apply(z1, 1:2, FUN=mean)
zm2 = apply(z1+z2+z3, 1:2, FUN=mean, na.rm=TRUE)

par(mfrow=c(3,4), mar=c(4,4,1,1))
for(i in 1:12) 
  image.plot(z2[,,i]/zm2, zlim=c(-0.6, 1.1))
