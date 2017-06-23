library(kali)
library(lubridate)
library(data.table)
library(mgcv)
library(parallel)

source("auxiliar_functions.R")

load("input/modelData.RData")

ncores = 4 
cl = makeCluster(ncores)

modelData = modelData[complete.cases(modelData),]
modelData$shelfF = as.factor(sign(modelData$shelf))
modelData$month = as.factor(month(num2date(modelData$time)))
modelData$monthn = month(num2date(modelData$time))

DateStamp()

# hyp: small data size, less predictive power 
# and more dispersion. Big data size? Less dispersion
# (one combination for full data). Replacement?

# experiment 1: split data for validation and training.
# sequential split of training data, same validation set
# 10 replicates per sample size: 1k, 5k, 10k, 50k, 100k, 
# 500k, 1M, 3M

set.seed(1234)

N = 3e6

tra = sample(nrow(modelData), N)
val = setdiff(1:nrow(modelData), tra)

training   = modelData[tra, ]
validation = modelData[val, ]

k = 100
sizes = 1e3*round(10^seq(0.8, 3, length=20))

# dim(output) = (val, replicates, model)
output = array(dim=c(length(val), k, length(sizes)))
models = list()

for(iSize in seq_along(sizes)) {
  DateStamp("Fitting model for n =", sizes[iSize])
  smodels = list()
  pb = txtProgressBar(style=3)
  for(i in seq_len(k)) {
    cat("\nReplicate", i, "\n")
    demo = training[sample(nrow(training), sizes[iSize]),]
    mod  = bam(modis ~ s(seawifs, k=40, bs="cr"), data=demo,
               cluster=cl)
    smodels[[i]] = trimGAM(mod)
    output[, i, iSize] = predict(mod, newdata=validation, 
                                 type="response", cluster=cl)
    setTxtProgressBar(pb, (k*(iSize-1) + i)/(length(sizes)*k))
  }
  names(smodels) = paste0("rep", 1:k)
  models[[as.character(sizes[iSize])]] = smodels
}

save(list = c("models", "output"), 
     file="output/exp1.RData")

err = validation$modis - output
mse = sqrt(apply(err^2, 2:3, FUN=mean))
rsq = apply(output, 2:3, FUN=cor, 
            y=validation$modis)

plot(sizes, colMeans(mse), log="x", type="b")
par(new=TRUE)
plot(sizes, colMeans(rsq), log="x", type="b",
     col="red", axes=FALSE)
axis(4)

xx = apply(do.call(rbind,
             lapply(models, FUN = function(x) t(sapply(x, FUN=explDev)))),
             2, FUN=function(x) {class(x)="numeric"; x})             
boxplot(expDev ~ n, data=xx)
boxplot(sqrt(MSE) ~ n, data=xx)

plot(sqrt(MSE) ~ n, data=xx, log="x")

