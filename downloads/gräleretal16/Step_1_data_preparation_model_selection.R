paper <- FALSE # TRUE
animate <- FALSE
set.seed(123)

# libraries
library(sp)
library(spacetime)
library(gstat)
library(rgdal)

# load intial working set (contained in the submitted zip-file):
load("DE_RB_2005.RData")
DE_RB_2005@sp <- spTransform(DE_RB_2005@sp, CRS("+init=epsg:32632"))

# load German boundaries
load(url("http://biogeo.ucdavis.edu/data/gadm2/R/DEU_adm1.RData"))
gadm <- spTransform(gadm, CRS("+init=epsg:32632"))
plot(gadm)

# station wise coverage
barplot(sort(table(DE_RB_2005@index[,1])),  
        main="reported days per station",
        ylab="number of days", xaxt="n")

# annual mean PM10 concentration
spplot(DE_RB_2005@sp,"annual_mean_PM10", col.regions=bpy.colors(),
       sp.layout = list("sp.polygons", gadm), scales=list(draw=T),
       key.space="right", colorkey=T,
       main=expression(paste("annual mean ","PM"[10]," concentration")))

# acf
acf(DE_RB_2005[sample(68,1),,drop=F]@data)
var(DE_RB_2005@data$PM10)

# a few daily snapshots
smplDays <- sort(sample(365,8))
if(paper)
  png("daily_means_PM10.png", width=9, height=6, "in", res=150)
stplot(as(DE_RB_2005[,smplDays],"STFDF"),
       col.regions=bpy.colors(120)[-(1:20)],
       sp.layout = list("sp.polygons", gadm), scales=list(draw=F), 
       key.space="right", colorkey=T, cuts=0:70,
       main=NULL) #expression(paste("daily mean ","PM"[10]," concentration")))
if(paper)
  dev.off()

# number of stations
length(DE_RB_2005@sp)

# calculate the empirical variogram
empVgm <- variogramST(PM10~1, DE_RB_2005, tlags=0:6)
plot(empVgm, wireframe=T, scales=list(arrows=F), col.regions=bpy.colors())
plot(empVgm)

# fit of theoretical purely spatial models #
############################################
spEmpVgm <- empVgm[empVgm$timelag == 0,]
class(spEmpVgm) <- c("gstatVariogram", "data.frame")
spEmpVgm <- spEmpVgm[-1,1:3]
spEmpVgm$dir.hor <- 0
spEmpVgm$dir.ver <- 0
spVgmMod <- fit.variogram(spEmpVgm, vgm(80,"Exp",300000,20))
plot(spEmpVgm, spVgmMod)

# fit of theoretical spatio-temporal models #
#############################################

linStAni <- estiStAni(empVgm, c(50000,200000))

plot(gamma~dist, empVgm[empVgm$timelag == 0,], ylim=c(0,100), xlim=c(0,800000))
points(empVgm[empVgm$spacelag == 0,]$timelag*linStAni, empVgm[empVgm$spacelag == 0,]$gamma, col="red")

##
# rescale empVgm and linStAni to km for estimation
empVgm$dist  <- empVgm$dist/1000
empVgm$avgDist  <- empVgm$avgDist/1000
empVgm$spacelag <- empVgm$spacelag/1000

linStAni <- linStAni/1000

# separable
separableModel <- vgmST("separable", 
                        space=vgm(0.9,"Exp", 200, 0.1),
                        time =vgm(0.9,"Sph", 3.5, 0.1),
                        sill=120)
fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
                               stAni = linStAni, method = "L-BFGS-B", 
                               control = list(parscale=c(100,1,10,1,100)),
                               lower = c(10,0,.1,0,0.1), 
                               upper = c(2000,1,12,1,200))
attr(fitSepModel, "optim.output")$value
# Exp+Exp: 9.87, Exp+Sph: 6.82, Sph+Exp: 10.42, Sph+Sph: 7.50
plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F),
     zlim=c(0,135))

# product-sum
prodSumModel <- vgmST("productSum",
                      space=vgm(10, "Exp", 200, 1),
                      time= vgm(10, "Sph",   2, 1), 
                      k=2)
fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                   stAni = linStAni, method = "L-BFGS-B", 
                                   control = list(parscale=c(1,10,1,1,0.1,1,10)),
                                   lower = rep(0.0001,7))
attr(fitProdSumModel, "optim.output")$value
# Exp+Exp: 10.09, Exp+Sph: 6.88, Sph+Exp: 10.64, Sph+Sph: 7.59
plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),
     zlim=c(0,135))

# metric
metricModel <- vgmST("metric",
                     joint = vgm(60, "Mat", 150, 10, kappa=0.6),
                     stAni = 60)
fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                                  stAni=linStAni, method = "L-BFGS-B",
                                  control = list(parscale = c(10,20,5,10)),
                                  lower=c(80,50,5,50),
                                  upper=c(200,1500,60,300))
attr(fitMetricModel, "optim.output")$value 
# Exp: 10.25, Sph: 10.59,
# Gau: 21.32, Mat 5: 18.20, Mat 2: 14.43, Mat 1.25: 12.04,
# Mat 1: 11.07, Mat 0.75: 10.23, Mat 0.6: 10.05
plot(empVgm, fitMetricModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))

# simplified sumMetric model?
sumMetricFromsimpleSumMetric <- function(vgm) {
  vgmST("sumMetric",
        space=vgm(vgm$space$psill, vgm$space$model, vgm$space$range, vgm$nugget/3),
        time =vgm(vgm$time$psill, vgm$time$model, vgm$time$range, vgm$nugget/3),
        joint=vgm(vgm$joint$psill, vgm$joint$model, vgm$joint$range, vgm$nugget/3),
        stAni=vgm$stAni)
}

pars.simple.l <- c(sill.s = 0, range.s = 10,
                   sill.t = 0, range.t = 0.1,
                   sill.st= 0, range.st= 10,
                   nugget=0, anis = 40)
pars.simple.u <- c(sill.s = 200,  range.s = 500,
                   sill.t = 200,  range.t = 20,
                   sill.st= 200, range.st = 5000,
                   nugget = 1E2, anis = 1E3)

pars.l <- c(sill.s = 0,  range.s = 10,  nugget.s = 0,
            sill.t = 0,  range.t = 0.1,   nugget.t = 0,
            sill.st= 0, range.st = 10, nugget.st = 0, 
            anis = 40)
pars.u <- c(sill.s = 200,  range.s = 1E3,  nugget.s = 20,
            sill.t = 200,  range.t = 75,   nugget.t = 20,
            sill.st= 200, range.st = 5E3, nugget.st = 20,
            anis = 500)

simpleSumMetricModel <- vgmST("simpleSumMetric",
                              space=vgm(120,"Sph", 150),
                              time =vgm(120,"Exp", 10),
                              joint=vgm(120,"Sph", 150),
                              nugget=10, stAni=150)
fitSimpleSumMetricModel <- fit.StVariogram(empVgm, simpleSumMetricModel,  
                                           fit.method = 7, stAni=linStAni,
                                           method = "L-BFGS-B", trace=F,
                                           lower=pars.simple.l, upper=pars.simple.u,
                                           control=list(parscale=c(1,10,
                                                                   1,1,
                                                                   1,100,
                                                                   1,10),
                                                        maxit=1e4))
attr(fitSimpleSumMetricModel, "optim.output")$value
# Exp+Exp+Exp: 4.10 Sph+Sph+Sph: 3.56
# Sph+Exp+Exp: 3.94 Exp+Sph+Exp: 3.60 Exp+Exp+Sph: 3.74
# Exp+Sph+Sph: 3.98 Sph+Exp+Sph: 3.31 Sph+Sph+Exp: 3.32
plot(empVgm,fitSimpleSumMetricModel, wireframe = T, scales = list(arrows = F), all = T , zlim=c(0,130))

# sum-metric
sumMetricModel <- sumMetricFromsimpleSumMetric(fitSimpleSumMetricModel)

sumMetricModel <- vgmST("sumMetric",
                        space = vgm(20, "Sph", 150, 1),
                        time = vgm(10, "Exp", 2, 0.5),
                        joint = vgm(80, "Sph", 1500, 2.5),
                        stAni = 120)

fitSumMetricModel <- fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                                     method = "L-BFGS-B", trace=F,
                                     lower=pars.l, upper=pars.u,
                                     control=list(parscale=c(1, 100,  1,
                                                             1, 0.5, 1,
                                                             1, 100,  1,
                                                             100),
                                                  maxit=1e4))
attr(fitSumMetricModel, "optim.output")$value
# Exp+Exp+Exp: 4.10 Sph+Sph+Sph: 3.36
# Sph+Exp+Exp: 3.89 Exp+Sph+Exp: 3.60 Exp+Exp+Sph: 3.74
# Exp+Sph+Sph: 3.73 Sph+Exp+Sph: 3.31 Sph+Sph+Exp: 3.32

plot(empVgm,fitSumMetricModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,130))


if(paper)
  png("allVgmsWireframe.png", 9, 6, "in", bg="white", res = 150)
plot(empVgm, list(fitSepModel, fitProdSumModel, fitMetricModel,
                  fitSumMetricModel, fitSimpleSumMetricModel), 
     wireframe=T, all=T, zlim=c(0,140), ylim=c(0,6.1), xlim=c(0,300),
     scales=list(arrows = F, cex=.8,
                 x=list(at=0:3*100), 
                 y=list(at=0:6, labels=c("0 ","","2 ","","4 ","","6 ")), 
                 z=list(at=0:5*25, labels=c("0  ","","50   ","","100    ",""))),
     at=0:100*1.4,
     xlab=list("space [km]", rot=27, cex=0.8), 
     ylab=list("time [days]", rot=-40, cex=0.8),
     zlab=list(NULL, rot=94, cex=0.8))
if(paper)
  dev.off()


library(lattice)
spacelag <- rep(0:300, 13)
timelag <- rep(0:12/2,each=301)  

if(paper)
  png("vgmVsMetricDist.png", 9, 6, "in", bg="white", res = 150)
contourplot(model~spacelag+timelag|type, 
            rbind(cbind(variogramSurface(fitSumMetricModel, 
                                         data.frame(spacelag=spacelag,
                                                    timelag=timelag)),
                        data.frame(type = rep("variogram of the sum-metric model", length(spacelag)))),
                  data.frame(spacelag=spacelag,
                             timelag=timelag,
                             model=sqrt(spacelag^2+116^2*timelag^2)/10,
                             type="metric distance [10 km]")),
            at=0:15*10,
            xlab="space [km]",
            ylab="timelag [days]")
if(paper)
  dev.off()
