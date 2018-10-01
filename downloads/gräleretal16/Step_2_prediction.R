# interpolation, build a grid over Germany
spDE <- DE_RB_2005@sp

ceiling(apply(spDE@bbox,1,diff)/10000)
spDE@bbox[,1]%/%10000*10000

gridDE <- SpatialGrid(GridTopology(c(305000, 5255000), c(10000,10000),
                                   cells.dim=c(60,84)))
proj4string(gridDE) <- CRS("+init=epsg:32632")
fullgrid(gridDE) <- F
gridDE <- gridDE[!is.na(over(geometry(gridDE), geometry(gadm)))]
plot(gridDE)

# back scale vgms:
fitSepModel$space$range <- fitSepModel$space$range*1000

fitProdSumModel$space$range <- fitProdSumModel$space$range*1000

fitMetricModel$joint$range <- fitMetricModel$joint$range*1000
fitMetricModel$stAni <- fitMetricModel$stAni*1000

fitSumMetricModel$space$range <- fitSumMetricModel$space$range*1000
fitSumMetricModel$joint$range <- fitSumMetricModel$joint$range*1000
fitSumMetricModel$stAni <- fitSumMetricModel$stAni*1000

fitSimpleSumMetricModel$space$range <- fitSimpleSumMetricModel$space$range*1000
fitSimpleSumMetricModel$joint$range <- fitSimpleSumMetricModel$joint$range*1000
fitSimpleSumMetricModel$stAni <- fitSimpleSumMetricModel$stAni*1000

# animation
if(animate) {
    DE_pred <- STF(gridDE, DE_RB_2005@time)
    
    predMat <- matrix(NA,0,2)
    for (rd in 15:180) {
      predMat <- rbind(predMat,
                       krigeST(PM10~1, data=DE_RB_2005[,rd+(-5:5)], # start: 12:24
                           newdata=DE_pred[,rd,drop=F], computeVar=T,
                           fitSumMetricModel, # nmax=50,
                           stAni=linStAni*1000/24/3600)@data)
      cat(rd,"\n")
    }
    
    DE_pred_winter <- DE_pred[,15:180]
    DE_pred_winter <- addAttrToGeom(DE_pred_winter, predMat)
    
    stplot(DE_pred_winter[,,"var1.pred"])
    stplot(DE_pred_winter[,,"var1.var"])
    
    for(i in 1:length(DE_pred_winter@time)) { # i <- 1
      pnt <- spplot(DE_pred_winter[,i], "var1.pred", col.regions=bpy.colors(160)[-(1:9)], scales=list(draw=F),
                    at=c(-2.8,0:150), 
                    sp.layout = list(list("sp.polygons", gadm, first=FALSE, col=gray(0.5)),
                                     list("sp.points", DE_RB_2005[,i+14], col=gray(0.25), pch=3, cex=.5)),
                    main=as.character(index(DE_pred_winter[,i,drop=F]@time)))
      
      png(file=paste("animate/pred",i,".png", sep=""), width=6, height=6, "in", res=150)
      print(pnt)
      dev.off()
    }
}

# prediction of sample days
DE_pred <- STF(gridDE, DE_RB_2005@time[smplDays])
tIDS <- unique(pmax(1,pmin(as.numeric(outer(-5:5, smplDays, "+")), 365)))

sepPred <- krigeST(PM10~1, data=DE_RB_2005[,tIDS], bufferNmax = 1,
                   newdata=DE_pred, fitSepModel, nmax=50,
                   stAni=linStAni*1000/24/3600)
psPred <- krigeST(PM10~1, data=DE_RB_2005[,tIDS], 
                  newdata=DE_pred, fitProdSumModel, nmax=50,
                  stAni=linStAni*1000/24/3600)
metPred <- krigeST(PM10~1, data=DE_RB_2005[,tIDS], 
                   newdata=DE_pred, fitMetricModel, nmax=50,
                   stAni=linStAni*1000/24/3600)
sumPred <- krigeST(PM10~1, data=DE_RB_2005[,tIDS], 
                   newdata=DE_pred, fitSumMetricModel, nmax=50,
                   stAni=linStAni*1000/24/3600)
smplSumPred <- krigeST(PM10~1, data=DE_RB_2005[,tIDS], # start: 12:24
                       newdata=DE_pred,  fitSimpleSumMetricModel, nmax=50,
                       stAni=linStAni*1000/24/3600)

sepPred@sp <- gridDE
psPred@sp <- gridDE
metPred@sp <- gridDE
sumPred@sp <- gridDE
smplSumPred@sp <- gridDE

pureSpPred <- matrix(NA, nrow=length(gridDE), length(smplDays))
col <- 1
for(i in smplDays) { # i <- 1
  pureSpPred[,col] <- krige(PM10~1, DE_RB_2005[,i],
                            gridDE, model=spVgmMod, nmax=50)$var1.pred
  col <- col+1
}

pureSpPred <- STFDF(gridDE, DE_RB_2005@time[smplDays], data.frame(var1.pred = as.numeric(pureSpPred)))
stplot(pureSpPred[,,"var1.pred"], main="pure spatial interpolation")

if(paper) {
  png(file="pred_daily_means_PM10.png", width=9, height=6, "in", res=150)
  stplot(smplSumPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
         main=NULL, at=0:70, # "spatio-temporal sum-metric model"
         sp.layout = list(list("sp.polygons", gadm, first=FALSE, col=gray(0.5)),
                          list("sp.points", DE_RB_2005@sp, col=gray(0.25), pch=3, cex=.5)))
  dev.off()
} else {
  stplot(pureSpPred, col.regions=bpy.colors, scales=list(draw=F),
         main="pure spatial daily prediction")
  stplot(sepPred, col.regions=bpy.colors(), scales=list(draw=F),
         main="spatio-temporal separable model")
  stplot(psPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal product-sum model")
  stplot(metPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal metric model")
  stplot(sumPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal sum-metric model")
  stplot(smplSumPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal simple sum-metric  model")
}
