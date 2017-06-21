library(kali)
library(lubridate)
library(data.table)
library(mgcv)
library(parallel)

source("paperPP_functions.R")

load("modelData.RData")

ncores = 4 
cl = makeCluster(ncores)

modelData = modelData[complete.cases(modelData),]
modelData$shelfF = as.factor(sign(modelData$shelf))
modelData$month = as.factor(month(num2date(modelData$time)))
modelData$monthn = month(num2date(modelData$time))

N = 1e4
ind = sample(nrow(modelData), N)
demo = modelData[ind, ]

val = setdiff(sample(nrow(modelData), 
                     min(10*N, nrow(modelData))), ind)

n = min(3*N, length(val))
test1 = modelData[sample(val, n), ]
# Base model --------------------------------------------------------------

DateStamp()
m0 = bam(modis ~ s(seawifs, k=40, bs="cr"), data=demo,
         method="REML")

DateStamp()
m0a = bam(modis ~ s(lon, by=seawifs, bs="ad"), data=demo,
          method="REML")
DateStamp()
m0a2 = bam(modis ~ te(lon, seawifs, k=25), data=demo,
          method="REML")
DateStamp()
m0b = bam(modis ~ s(lat, by=seawifs, bs="ad"), data=demo,
          method="REML")
DateStamp()
m0b2 = bam(modis ~ te(seawifs, lat, k=25), data=demo,
          method="REML")
DateStamp()
m0c = bam(modis ~ s(lon, lat, by=seawifs, k=169), data=demo,
         method="REML")
DateStamp()
m0c2 = bam(modis ~ te(lon, lat, by=seawifs, k=25), data=demo,
          method="REML")

DateStamp()
m0d = bam(modis ~ s(shelf, by=seawifs, k=40, bs="cr"), data=demo,
          method="REML")
DateStamp()
m0e = bam(modis ~ te(shelf, lat, by=seawifs, k=30), data=demo,
          method="REML")
DateStamp()

explDev(m0, m0a, m0a2, m0b, m0b2, m0c, m0c2, m0d, m0e)
# 
# par(mfrow=c(2,2), mar=c(3,3,1,1))
# x = plot(m1d, se=FALSE)[[1]]
# image.map(x$x, x$y, matrix(x$fit, ncol=length(x$y)),
#           zlim=c(0.2,1.5))
# x = plot(m0b, se=FALSE)[[1]]
# image.map(x$x, x$y, matrix(x$fit, ncol=length(x$y)),
#           zlim=c(0.2,1.5))
# 
# x = plot(m0e, se=FALSE)[[1]]
# image.plot(x$x, x$y, matrix(x$fit, ncol=length(x$y)))
# 
# DateStamp()
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# gam.check(m1)


# Using decomposition -----------------------------------------------------

m1 = bam(modis ~ s(sTrend, k=30, bs="cr") + 
           s(sSeason, k=30, bs="cr") + 
           s(sRemain, k=30, bs="cr"), data=demo)

m1a = bam(modis ~ s(lon, by=sTrend, k=100, bs="cr") + 
           s(lon, by=sSeason, k=70, bs="cr") + 
           s(lon, by=sRemain, k=70, bs="cr"), data=demo)

m1a2 = bam(modis ~ te(lon, sTrend, k=12) + 
             te(lon, sSeason, k=12) + 
             te(lon, sRemain, k=12), data=demo)

m1b = bam(modis ~ s(lat, by=sTrend, k=50, bs="cr") + 
            s(lat, by=sSeason, k=50, bs="cr") + 
            s(lat, by=sRemain, k=50, bs="cr"), data=demo)

m1b2 = bam(modis ~ te(lat, sTrend, k=12) + 
             te(lat, sSeason, k=12) + 
             te(lat, sRemain, k=12), data=demo)

m1c = bam(modis ~ s(shelf, by=sTrend, k=50, bs="cr") + 
            s(shelf, by=sSeason, k=50, bs="cr") + 
            s(shelf, by=sRemain, k=50, bs="cr"), data=demo)

m1c2 = bam(modis ~ te(shelf, sTrend, k=12) + 
             te(shelf, sSeason, k=12) + 
             te(shelf, sRemain, k=12), data=demo)

m1d = bam(modis ~ s(lon, lat, by=sTrend, k=200) + 
            s(lon, lat, by=sSeason, k=200) + 
            s(lon, lat, by=sRemain, k=200), data=demo)

m1d2 = bam(modis ~ te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)

m1d2s = bam(modis ~ te(lon, lat, by=sTrend, k=15) + 
             te(lon, lat, by=sSeason, k=15) + 
             te(lon, lat, by=sRemain, k=15), data=demo,
            select=TRUE)

m1e = bam(modis ~ te(lon, lat, sTrend, k=10) + 
              te(lon, lat, sSeason, k=10) + 
              te(lon, lat, sRemain, k=10), data=demo,
          cluster=cl)


explDev(m0, m0c2, m1, m1a, m1a2, m1b, m1b2, 
        m1c, m1c2, m1d, m1d2, m1d2s, m1e)

# adding covariates

m2a = bam(modis ~ s(shelf, k=50, bs="cr") + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)

m2b = bam(modis ~ te(shelf, lat, k=15) + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)

m2c = bam(modis ~ s(lon, lat, k=144) + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)

DateStamp()
m2d = bam(modis ~ s(lon, lat, by=month, k=144) + month +  
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)

DateStamp()
m2d2 = bam(modis ~ s(lon, lat, by=month, k=144, id="month") + month + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)
DateStamp()

# intercept covariates

DateStamp()
m3a = bam(modis ~ ti(lon, by=month) + ti(lat, by=month) +
            ti(lon, lat, by=month, k=15, id="month") + month + 
             te(lon, lat, by=sTrend, k=15) + 
             te(lon, lat, by=sSeason, k=15) + 
             te(lon, lat, by=sRemain, k=15), data=demo)
DateStamp()

DateStamp()
m3b = bam(modis ~ ti(lon, by=month) + ti(lat, by=month) +
            ti(lon, lat, by=month, k=15, id="month") + month +
            s(shelf, by=month, k=50, bs="cr") +
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)
DateStamp()

DateStamp()

DateStamp()
m3c = bam(modis ~ ti(lon, by=month) + ti(lat, by=month) +
            ti(lon, lat, by=month, k=15, id="month") + month +
            s(shelf, by=month, k=50, bs="cr") +
            s(sTrend, by=shelfF, k=50, bs="cr") +
            s(sSeason, by=shelfF, k=50, bs="cr") + 
            s(sRemain, by=shelfF, k=50, bs="cr") + shelfF + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo)
DateStamp()

DateStamp()
m3c2 = bam(modis ~ te(lon, lat, by=month, k=15, id="month") + 
            month + s(shelf, by=month, k=50, bs="cr") +
            s(sTrend, by=shelfF, k=50, bs="cr") +
            s(sSeason, by=shelfF, k=50, bs="cr") + 
            s(sRemain, by=shelfF, k=50, bs="cr") + shelfF + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo,
          select=TRUE)
DateStamp()
m3d = bam(modis ~ te(lon, lat, by=month, k=15, id="month") + 
             month + s(shelf, k=50, bs="cr") +
             s(sTrend, by=shelfF, k=50, bs="cr") +
             s(sSeason, by=shelfF, k=50, bs="cr") + 
             s(sRemain, by=shelfF, k=50, bs="cr") + shelfF + 
             te(lon, lat, by=sTrend, k=15) + 
             te(lon, lat, by=sSeason, k=15) + 
             te(lon, lat, by=sRemain, k=15), data=demo,
          cluster = cl)
DateStamp()
m3e = bam(modis ~ te(lon, lat, k=15) + 
             month + s(shelf, by=month, k=50, bs="cr") +
             s(sTrend, by=shelfF, k=50, bs="cr") +
             s(sSeason, by=shelfF, k=50, bs="cr") + 
             s(sRemain, by=shelfF, k=50, bs="cr") + shelfF + 
             te(lon, lat, by=sTrend, k=15) + 
             te(lon, lat, by=sSeason, k=15) + 
             te(lon, lat, by=sRemain, k=15), data=demo,
           cluster = cl)
DateStamp()
m3f = bam(modis ~ te(lon, lat, k=15) + 
            month + s(shelf, k=50, bs="cr") +
            s(sTrend, by=shelfF, k=50, bs="cr") +
            s(sSeason, by=shelfF, k=50, bs="cr") + 
            s(sRemain, by=shelfF, k=50, bs="cr") + shelfF + 
            te(lon, lat, by=sTrend, k=15) + 
            te(lon, lat, by=sSeason, k=15) + 
            te(lon, lat, by=sRemain, k=15), data=demo,
          cluster = cl)
DateStamp()


explDev(m0, m0c2, m1, m1a, m1a2, m1b, m1b2, 
        m1c, m1c2, m1d, m1d2, m2a, m2b, m2c, m2d,
        m2d2, m3a, m3b, m3c, m3c2, m3d, m3e, m3f,
        m5i, m5j)


xp = plot(m2d, plot=FALSE)

par(mfrow=c(3,4), mar=c(3,3,1,1), oma=c(1,1,1,1))
for(i in 1:12) {
  x = xp[[i]]
  image.map(x$x, x$y, matrix(x$fit, ncol=length(x$y)),
            zlim=c(-800, 800))
}

par(mfrow=c(1,3), mar=c(3,3,1,1), oma=c(1,1,1,1))
for(i in 13:15) {
  x = xp[[i]]
  image.map(x$x, x$y, matrix(x$fit, ncol=length(x$y)),
            zlim=c(-0.5,1.6))
}

par(mar=c(3,3,1,1), mfrow=c(1,2))
check4(m1a2)
check4(m1a3)


gam.check(m1a2)

DateStamp()
# adding covariates
m2a = bam(modis ~ s(sTrend) + s(sSeason) + s(sRemain) +
               s(lat) + s(lon), data=demo)
DateStamp()
m2b = bam(modis ~ s(sTrend) + s(sSeason) + s(sRemain) +
               s(lat) + s(dc), data=demo)
DateStamp()
m2c = bam(modis ~ s(sTrend) + s(sSeason) + s(sRemain) +
               s(lat) + s(shelf), data=demo)
DateStamp()
m2d = bam(modis ~ s(sTrend) + s(sSeason) + s(sRemain) +
              s(lat) + s(lon) + s(dc) + s(shelf), data=demo)
DateStamp()
# interaction with shelf: trend
m3a = bam(modis ~ te(sTrend, shelf) + s(sSeason) + s(sRemain) +
              s(lat) + s(lon), data=demo)
DateStamp()
m3b = bam(modis ~ s(sTrend, by=shelfF) + s(sSeason) + s(sRemain) +
              shelfF + s(lat) + s(lon), data=demo)
DateStamp()
# interaction with shelf: all variables
m4a = bam(modis ~ s(sTrend, by=shelfF) + 
               s(sSeason, by=shelfF) + 
               s(sRemain, by=shelfF) + shelfF +
               s(lat) + s(lon), data=demo)
DateStamp()
m4b = bam(modis ~ s(sTrend, shelf) + 
               s(sSeason, shelf) + 
               s(sRemain, shelf) +
               s(lat) + s(lon), data=demo)
DateStamp()
m4c = bam(modis ~ s(sTrend, by=shelfF) + 
               s(sSeason, by=shelfF) + 
               s(sRemain, by=shelfF) + shelfF, data=demo)
DateStamp()
# adding seasonality
m5a = bam(modis ~ s(sTrend, shelf, by=month) + 
               s(sSeason) + 
               s(sRemain) + month +
               s(lat) + s(lon), data=demo)
DateStamp()
m5b = bam(modis ~ s(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               s(sRemain) + month + shelfF + 
               s(lat) + s(lon), data=demo)
DateStamp()
m5c = bam(modis ~ s(sTrend, shelf, by=month) + 
               s(sSeason, shelf) + 
               s(sRemain, shelf) + month +  
               s(lat) + s(lon), data=demo)
DateStamp()
m5d = bam(modis ~ s(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               s(sRemain, by=month) + month + shelfF + 
               s(lat) + s(lon), data=demo)
DateStamp()
m5e = bam(modis ~ s(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               s(sRemain, shelf, by=month) + month + shelfF + 
               s(lat) + s(lon), data=demo)
DateStamp()
m5f = bam(modis ~ te(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               te(sRemain, shelf, by=month) + month + shelfF + 
               s(lat) + s(lon), data=demo)
DateStamp()
m5g = bam(modis ~ te(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               te(sRemain, shelf, by=month) + month + shelfF + 
               s(lat, by=month) + s(lon), data=demo)
DateStamp()
m5h = bam(modis ~ te(sTrend, monthn, by=shelfF) + 
               s(sSeason, by=shelfF) + 
               te(sRemain, monthn, by=shelfF) +
               shelfF + te(lat, lon, by=month), data=demo)
DateStamp()
m5i = bam(modis ~ te(sTrend, shelf, by=month) + 
               s(sSeason, by=shelfF) + 
               te(sRemain, shelf, by=month) + month + shelfF + 
               te(lat, lon, by=month), data=demo, cluster=cl)

DateStamp()
m5j = bam(modis ~ te(sTrend, shelf, by=month) + 
               te(sSeason, shelf, by=month) + 
               te(sRemain, shelf, by=month) + month + 
               te(lat, lon, by=month), data=demo, cluster=cl)
DateStamp()

explDev(m0, m0b, m0c)
explDev(m0, m1, m2a, m2b, m2c, m2d)
explDev(m0, m1, m3a, m3b)
explDev(m0, m1, m4a, m4b, m4c)
explDev(m0, m1, m5a, m5b, m5c)
explDev(m0, m1, m5d, m5e, m5f)
explDev(m0, m1, m5g, m5h, m5i, m5j)





z0 = predict(m0, type="response")
z1 = predict(m1, type="response")
# z2 = predict(m2, type="response")
z3a = predict(m3a, type="response")
z3b = predict(m3b, type="response")
z3c = predict(m3c, type="response")


par(mfrow=c(2,3), mar=c(3,3,1,1))
plot(m3c)

# predict modis using seaWIFS, concantenate merged series (ncdf)

# lat: 20ºS - 6ºN
# lon: 93ºW - 70ºW

