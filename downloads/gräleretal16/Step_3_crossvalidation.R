## cross-validation
crossStat <- function(var1, var2="PM10", STxDF=DE_RB_2005) {
  diff <- STxDF[,,var1,drop=F]@data[[1]] - STxDF[,,var2,drop=F]@data[[1]]
  RMSE <- sqrt(mean(diff^2))
  MAE <- mean(abs(diff))
  ME <- mean(diff)
  COR <- cor(STxDF[,,var1,drop=F]@data[[1]], STxDF[,,var2,drop=F]@data[[1]])
  res <- c(RMSE, MAE, ME, COR)
  names(res) <- c("RMSE", "MAE", "ME", "COR")
  return(res)
}

# purely spatial:
pureSp <- NULL
for(i in 1:365) { # i <- 1
  pureSp <- c(pureSp, krige.cv(PM10~1,DE_RB_2005[,i,"PM10"],model=spVgmMod,nmax=10)$var1.pred)
}

DE_RB_2005@data$pureSp <- pureSp

pureSp <- NULL
for(i in 1:365) { # i <- 1
  pureSp <- c(pureSp, krige.cv(PM10~1,DE_RB_2005[,i,"PM10"],model=spVgmMod,nmax=50)$var1.pred)
}

DE_RB_2005@data$pureSp50Nghbr <- pureSp

plot(DE_RB_2005@data$pureSp,DE_RB_2005@data$PM10)
abline(0,1,col="red")

crossStat("pureSp")
crossStat("pureSp50Nghbr")

## spatio-temporal LOOCV

target <- as(DE_RB_2005[,,"PM10"],"STFDF")

## seprable model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) {
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                                    newdata=DE_RB_2005[loc,,drop=F], 
                                                    fitSepModel, nmax=10, progress=F,
                                                    stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$sepModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSepModel, nmax=50, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$sepModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## product-sum model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitProdSumModel, nmax=10, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$psModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitProdSumModel, nmax=50, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$psModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitMetricModel, nmax=10, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$metricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitMetricModel, nmax=50, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$metricModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]

## sum-metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSumMetricModel, nmax=10, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$sumMetricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]
crossStat("sumMetricModel10Nghbr")

# 50 neighbours
res <- array(NA, c(length(DE_RB_2005@sp), 365,2))
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 38
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"],] <- as.matrix(krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                                     newdata=DE_RB_2005[loc,,drop=F], 
                                                     fitSumMetricModel, nmax=50, progress=F,
                                                     computeVar=T,
                                                     stAni=linStAni*1000/24/3600)@data[,c("var1.pred","var1.var")])
}

DE_RB_2005@data$sumMetricModel50Nghbr <- as.vector(res[,,1])[!is.na(target@data)]
DE_RB_2005@data$sumMetricModel50NghbrVar <- as.vector(res[,,2])[!is.na(target@data)]
DE_RB_2005@data$sumMetricModel50Nghbr95u <- apply(DE_RB_2005@data, 1, 
                                                  function(x) {
                                                    qnorm(0.975, x["sumMetricModel50Nghbr"],
                                                          sqrt(x["sumMetricModel50NghbrVar"]))
                                                  })
DE_RB_2005@data$sumMetricModel50Nghbr95l <- apply(DE_RB_2005@data, 1, 
                                                  function(x) {
                                                    qnorm(0.025, x["sumMetricModel50Nghbr"],
                                                          sqrt(x["sumMetricModel50NghbrVar"]))
                                                  })

## simple sum-metric model
# 10 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSimpleSumMetricModel, nmax=10, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$simpleSumMetricModel10Nghbr <- as.vector(res)[!is.na(as.vector(res))]

# 50 neighbours
res <- matrix(NA, length(DE_RB_2005@sp), 365)
for(loc in 1:length(DE_RB_2005@sp)) { # loc <- 1
  cat("Location", loc, "\n")
  res[loc,!is.na(target[loc,])[,"PM10"]] <- krigeST(PM10~1, data=DE_RB_2005[(1:length(DE_RB_2005@sp))[-loc],], 
                                           newdata=DE_RB_2005[loc,,drop=F], 
                                           fitSimpleSumMetricModel, nmax=50, progress=F,
                                           stAni=linStAni*1000/24/3600)$var1.pred
}

DE_RB_2005@data$simpleSumMetricModel50Nghbr <- as.vector(res)[!is.na(as.vector(res))]


###
# cross-stats

round(rbind(
crossStat("pureSp"),
crossStat("sepModel10Nghbr"),
crossStat("psModel10Nghbr"),
crossStat("metricModel10Nghbr"),
crossStat("sumMetricModel10Nghbr"),
crossStat("simpleSumMetricModel10Nghbr"),


crossStat("pureSp50Nghbr"),
crossStat("sepModel50Nghbr"),
crossStat("psModel50Nghbr"),
crossStat("metricModel50Nghbr"),
crossStat("sumMetricModel50Nghbr"),
crossStat("simpleSumMetricModel50Nghbr")),2)

texRow <- function(x) {
  paste(paste(x,collapse=" & ")," \\\\ \n")
}

cat(apply(round(rbind(crossStat("pureSp"),
                      
                      crossStat("sepModel10Nghbr"),
                      crossStat("psModel10Nghbr"),
                      crossStat("metricModel10Nghbr"),
                      crossStat("sumMetricModel10Nghbr"),
                      crossStat("simpleSumMetricModel50Nghbr"),
                      
                      crossStat("pureSp50Nghbr"),
                      
                      crossStat("sepModel50Nghbr"),
                      crossStat("psModel50Nghbr"),
                      crossStat("metricModel50Nghbr"),
                      crossStat("sumMetricModel50Nghbr"),
                      crossStat("simpleSumMetricModel50Nghbr")
                      ),
                2),1,texRow))


attr(fitSepModel, "optim.output")$value # 6.82
attr(fitProdSumModel, "optim.output")$value # 6.88
attr(fitMetricModel, "optim.output")$value # 10.05
attr(fitSumMetricModel, "optim.output")$value # 3.31
attr(fitSimpleSumMetricModel, "optim.output")$value # 3.31

## cross-plots
library(lattice)
xyplot(pureSp ~ PM10, DE_RB_2005@data, asp=1)
xyplot(sepModel10Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(psModel10Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(metricModel10Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(sumMetricModel10Nghbr ~ PM10, DE_RB_2005@data, asp=1)

xyplot(pureSp50Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(sepModel50Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(psModel50Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(metricModel50Nghbr ~ PM10, DE_RB_2005@data, asp=1)
xyplot(sumMetricModel50Nghbr ~ PM10, DE_RB_2005@data, asp=1)

## cross stat per station
for (model in colnames(DE_RB_2005@data)[-1]) {
  cvStation <- numeric(length(DE_RB_2005@sp))
  
  for (i in 1:length(DE_RB_2005@sp)) {
    diff <- DE_RB_2005[i,,"PM10",drop=F]@data[[1]] - DE_RB_2005[i,,model,drop=F]@data[[1]]
    cvStation[i] <- sqrt(mean(diff^2))
  }
  
  DE_RB_2005@sp@data[[paste(model,"RMSE",sep="_")]] <- cvStation
}

spplot(DE_RB_2005@sp, "sumMetricModel50Nghbr_RMSE", cuts=2+0:6*2,
       sp.layout = list("sp.polygons", gadm), key.space="right",
       pch=c(16,16,16,16,8,8))

spplot(DE_RB_2005@sp, "annual_mean_PM10",
       sp.layout = list("sp.polygons", gadm), key.space="right")

cvAllVarStation <- as.matrix(DE_RB_2005@sp@data[,10:19])
colnames(cvAllVarStation) <- c("pur_10", "pur_50", "sep_10","sep_50","prs_10","prs_50","met_10","met_50","sum_10","sum_50")
rownames(cvAllVarStation) <- NULL
round(cvAllVarStation,2)
levelplot(cvAllVarStation[,c(1,3,5,7,9,2,4,6,8,10)], 
          main="RMSE per Station", xlab="station", ylab="variogram model", 
          col.regions=rev(heat.colors(20))[-c(1:3)],asp=0.5)


loc <- 38 # sample(68,1) # 15
tw <- "2005-01-15/2005-04-15"

if(paper) 
  png("singleStationTimeSeries.png", 9, 4, "in", bg="white", res = 149)
plot(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr"], main=paste("Location", 
                                                              DE_RB_2005@sp@data$station_european_code[loc]),
     ylim=c(0,70))
points(DE_RB_2005[loc,tw][,"PM10"], type="l", col="darkgreen", lty=1)
points(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr95u"], type="l", col="darkgrey", lty=2)
points(DE_RB_2005[loc,tw][,"sumMetricModel50Nghbr95l"], type="l", col="darkgrey", lty=2)
legend("topright",legend = c("observed","sum-metric","95 % prediction band"), lty=c(1,1,2),
     col=c("darkgreen", "black", "darkgrey") )
if(paper)
  dev.off()
