## 4/22/2018
## Reproduction code for the paper:
## spBayesSurv: Fitting Bayesian Spatial Survival Models Using R
## R version 3.3.3 (2017-03-06) -- "Another Canoe"
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## spBayesSurv version: 1.1.3

## NOTE: To save computation time, one can use small chains to check 
##       for errors. For example, replace
##       mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
##       with 
##       mcmc <- list(nburn = 500, nsave = 200, nskip = 0, ndisplay = 1000)
##       There 8 of them. 
##       However, such replacement will not reproduce the results presented in the paper

#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the PO model with ICAR frailties 
## in Section 2.4.1 of the paper
#####-----------------------------------------------------######
rm(list = ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data("LeukSurv")
d <- LeukSurv[order(LeukSurv$district), ]
head(d)
## get adjacency matrix
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", 
                                  package="spBayesSurv"))
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)

##-------------PO-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
#prior <- list(maxL = 15, beta0 = rep(0,3), S0 = diag(10,3), theta0 = rep(0,2), V0 = diag(10,2),
#           a0 = 1, b0 = 1, taua0 = 1, taub0 = 1, phia0 = 1, phib0 = 1)
prior <- list(maxL = 15)
ptm <- proc.time()
res1 <- survregbayes(formula = Surv(time, cens) ~ age + sex + wbc + tpi +
                       frailtyprior("car", district), data = d, survmodel = "PO",
                     dist = "loglogistic", mcmc = mcmc, prior = prior, Proximity = E)
(systime1 <- proc.time() - ptm) 
(sfit1 <- summary(res1)) 

save.image("Results_PO_car.RData")

############# BFs for testing loglogistic baseline #############.
res1$BF.baseline

############# Trace plots #############
pdf(file = "Leukemia-PO-CAR-trace.pdf", paper = "special", width = 10, height = 5)
par(mfrow = c(2,3))
par(cex = 1, mar = c(2.5, 4.1, 1, 1))
traceplot(mcmc(res1$beta[1,]), xlab = "", main = "age")
traceplot(mcmc(res1$beta[2,]), xlab = "", main = "sex")
traceplot(mcmc(res1$beta[3,]), xlab = "", main = "wbc")
traceplot(mcmc(res1$beta[4,]), xlab = "", main = "tpi")
traceplot(mcmc(res1$tau2), xlab = "", main = "tau^2")
traceplot(mcmc(res1$alpha), xlab = "", main = "alpha")
dev.off()

############# Cox-Snell residuals #############
set.seed(1)
pdf(file = "Leukemia-PO-CAR-Cox-Snell.pdf", paper="special", width=8, height=8)
cox.snell.survregbayes(res1, ncurves = 10)
dev.off()

############# Survival Curves #############
## plot
pdf(file ="Leukemia-PO-CAR-age.pdf", paper="special", width=8, height=8);
par(mfrow=c(1, 1)) 
tgrid <- seq(0.1, 5000, length.out = 300)
# age effect
xpred <- data.frame(age = c(49, 65, 74),
                    sex = c(0, 0, 0),
                    wbc = c(38.59, 38.59, 38.59),
                    tpi = c(0.3398, 0.3398, 0.3398),
                    row.names = c("age=49", "age=65", "age=74"))
plot(res1, xnewdata = xpred, tgrid = tgrid, cex = 2)
dev.off()

########### frailty maps ##################
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))
frail0 <- (rowMeans(res1$v))
frail <- frail0[as.integer(names(nwengland))]
values <- cbind(as.integer(names(nwengland)), frail) 
op <- par(no.readonly = TRUE)
pdf(file = "Leukemia-PO-CAR-Map.pdf", paper = "special", width = 8, height = 8)
par(mar = c(3, 0, 0, 0))
plotmap(nwengland, x = values, col = (gray.colors(10, 0.3, 1))[10:1], 
        pos = "bottomleft", width = 0.5, height = 0.04)
dev.off()

#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the PO model with GFR frailties
## in Section 2.4.2 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data(LeukSurv);
d <- LeukSurv[order(LeukSurv$district), ]
head(d)

##-------------PO-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000) 
#prior=list(maxL=15,beta0=rep(0,3),S0=diag(10,3),theta0=rep(0,2),V0=diag(10,2),
#           a0=1,b0=1,taua0=1,taub0=1,phia0=1,phib0=1)
prior <- list(maxL = 15, nu = 1, taua0 = 2,taub0 = 10,
              nknots = 100, nblock = 1043) 
d$ID=1:nrow(d);
locations=cbind(d$xcoord,d$ycoord);
ptm <- proc.time()
res2 = survregbayes(formula=Surv(time,cens)~age+sex+wbc+tpi
                    +frailtyprior("grf",ID), data=d, survmodel="PO",
                    dist="loglogistic", mcmc=mcmc, prior=prior,
                    Coordinates=locations);
(systime2 <- proc.time() - ptm) 
(sfit2 <- summary(res2)) 

save.image("Results_PO_grf.RData")

############# BFs for testing loglogistic baseline #############.
res2$BF.baseline

############# Trace plots #############
pdf(file ="Leukemia-PO-GRF-trace.pdf", paper="special", width=10, height=5)
par(mfrow=c(2,3));
par(cex=1,mar=c(2.5,4.1,1,1))
traceplot(mcmc(res2$beta[1,]), xlab="", main="age")
traceplot(mcmc(res2$beta[2,]), xlab="", main="sex")
traceplot(mcmc(res2$beta[3,]), xlab="", main="wbc")
traceplot(mcmc(res2$beta[4,]), xlab="", main="tpi")
traceplot(mcmc(res2$tau2), xlab="", main="tau^2")
traceplot(mcmc(res2$phi), xlab="", main="phi")
dev.off()

############# Cox-Snell residuals #############
set.seed(1)
pdf(file = "Leukemia-PO-GRF-Cox-Snell.pdf", paper="special", width=8, height=8)
cox.snell.survregbayes(res2, ncurves = 10)
dev.off()

############# Survival Curves #############
## plot
pdf(file ="Leukemia-PO-GRF-age.pdf", paper="special", width=8, height=8);
par(mfrow=c(1, 1)) 
tgrid <- seq(0.1, 5000, length.out = 300)
# age effect
xpred <- data.frame(age = c(49, 65, 74),
                    sex = c(0, 0, 0),
                    wbc = c(38.59, 38.59, 38.59),
                    tpi = c(0.3398, 0.3398, 0.3398),
                    row.names = c("age=49", "age=65", "age=74"))
plot(res2, xnewdata = xpred, tgrid = tgrid, cex = 2)
dev.off()

########### frailty maps ##################
nwengland=read.bnd(system.file("otherdata/nwengland.bnd", package="spBayesSurv"));
frail= round((rowMeans(res2$v)),3); nclust=5;
frail.cluster = cut(frail, breaks = nclust);
frail.names = names(table(frail.cluster))
rbPal <- colorRampPalette(c('blue','red'))
frail.colors=rbPal(nclust)[as.numeric(frail.cluster)]
pdf(file ="Leukemia-PO-GRF-Map.pdf", paper="special", width=8, height=8)
par(mar=c(3,0,0,0))
plot(nwengland)
points(cbind(d$xcoord,d$ycoord), col=frail.colors)
legend("topright",title="frailty values",legend=frail.names,
       col=rbPal(nclust),pch=20,cex=1.7)
dev.off()

#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for variable selection under the PO model with ICAR frailties
## in Section 2.5 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data("LeukSurv")
d <- LeukSurv[order(LeukSurv$district), ]
head(d)
## get adjacency matrix
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", 
                                  package="spBayesSurv"))
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)

##-------------PO-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
#prior=list(maxL=15,beta0=rep(0,3),S0=diag(10,3),theta0=rep(0,2),V0=diag(10,2),
#           a0=1,b0=1,taua0=1,taub0=1,phia0=1,phib0=1)
prior=list(maxL=15);
ptm<-proc.time()
res3 = survregbayes(formula=Surv(time,cens)~age+sex+wbc+tpi
                    +frailtyprior("car",district),data=d,survmodel="PO",
                    dist="loglogistic",mcmc=mcmc,prior=prior,Proximity=E,
                    selection=TRUE);
(sfit3=summary(res3)) 
(systime3=proc.time()-ptm) 

save.image("Results_PO_car_selection.RData")


#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the log-logistic PO model with ICAR frailties 
## in Section 2.6 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data("LeukSurv")
d <- LeukSurv[order(LeukSurv$district), ]
head(d)
## get adjacency matrix
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", 
                                  package="spBayesSurv"))
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)

##-------------PO-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
#prior=list(maxL=15,beta0=rep(0,3),S0=diag(10,3),theta0=rep(0,2),V0=diag(10,2),
#           a0=1,b0=1,taua0=1,taub0=1,phia0=1,phib0=1)
prior=list(maxL=15, a0=-1,thete0=rep(0,2),V0=diag(1e10,2));
state=list(alpha=Inf);
ptm<-proc.time()
res11 = survregbayes(formula=Surv(time,cens)~age+sex+wbc+tpi
                     +frailtyprior("car",district),data=d,survmodel="PO",
                     dist="loglogistic",mcmc=mcmc,prior=prior,state=state,
                     Proximity=E,InitParamMCMC=FALSE);
(sfit11=summary(res11)) 
(systime11=proc.time()-ptm) 

save.image("Results_PO_car_parametric.RData")


#####-----------------------------------------------------######
## Analysis of the PBC Data with time-dependent covariates
## Below is R code for fitting the PH model 
## in Section 2.7 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")

############## Read data ##############
temp <- subset(pbc, id <= 312, select=c(id:sex, stage)) # baseline data
pbc2 <- tmerge(temp, temp, id=id, endpt = event(time, status))
pbc2 <- tmerge(pbc2, pbcseq, id=id, ascites = tdc(day, ascites),
               bili = tdc(day, bili), albumin = tdc(day, albumin),
               protime = tdc(day, protime), alk.phos = tdc(day, alk.phos))
pbc2 = pbc2[,c("id","tstart","tstop","endpt","bili","protime")];
head(pbc2);
coxph(Surv(tstart,tstop,endpt==2)~log(bili)+log(protime), data=pbc2)

##-------------PH-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
ptm<-proc.time()
fit1 = survregbayes(Surv(tstart,tstop,endpt==2)~log(bili)+log(protime), 
                    data=pbc2, survmodel="PH", dist="loglogistic", mcmc=mcmc,
                    subject.num=id);
fit1
systime1=proc.time()-ptm; systime1;
save.image("Results_PBC_PH.RData")

## another way
pbc2$tleft=pbc2$tstop; pbc2$tright=pbc2$tstop;
pbc2$tright[which(pbc2$endpt!=2)]=NA;
pbc2$ltruncation=pbc2$tstart;
head(pbc2)
ptm<-proc.time()
fit11 = survregbayes(Surv(tleft,tright,type="interval2")~log(bili)
                     +log(protime), data=pbc2, survmodel="PH", 
                     dist="loglogistic", mcmc=mcmc,
                     truncation_time=ltruncation, subject.num=id);
sfit11=summary(fit11); sfit11
systime11=proc.time()-ptm; systime11;


#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the GAFT model with ICAR frailties 
## in Section 3.4 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data("LeukSurv")
d <- LeukSurv[order(LeukSurv$district), ]
head(d)
## get adjacency matrix
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", 
                                  package="spBayesSurv"))
adj.mat <- bnd2gra(nwengland)
E <- diag(diag(adj.mat)) - as.matrix(adj.mat)

##-------------GAFT-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior=list(maxL=4, a0=5, b0=1);
ptm<-proc.time()
res1 = frailtyGAFT(formula=Surv(time,cens)~age + sex + wbc + tpi + 
                     baseline(age, sex, wbc, tpi) +
                     frailtyprior("car",district),
                   data=d,mcmc=mcmc,prior=prior,Proximity=E);
sfit1=summary(res1); sfit1
systime1=proc.time()-ptm; systime1;

save.image("Result_GAFT_car.RData")

############# Trace plots #############
pdf(file ="Leukemia-GAFT-CAR-trace.pdf", paper="special", width=8, height=8)
par(mfrow=c(2,3));
par(cex=1,mar=c(2.5,4.1,1,1))
traceplot(mcmc(res1$beta[1,]), xlab="", main="intercept")
traceplot(mcmc(res1$beta[2,]), xlab="", main="age")
traceplot(mcmc(res1$beta[3,]), xlab="", main="sex")
traceplot(mcmc(res1$beta[4,]), xlab="", main="wbc")
traceplot(mcmc(res1$beta[5,]), xlab="", main="tpi")
traceplot(mcmc(res1$tau2), xlab="", main="tau^2")
dev.off()

############# Survival Curves #############
pdf(file ="Leukemia-GAFT-CAR-age.pdf", paper="special", width=8, height=8);
par(mfrow=c(1, 1)) 
tgrid <- seq(0.1, 5000, length.out = 300)
# age effect
xpred <- data.frame(age = c(49, 65, 74),
                    sex = c(0, 0, 0),
                    wbc = c(38.59, 38.59, 38.59),
                    tpi = c(0.3398, 0.3398, 0.3398),
                    row.names = c("age=49", "age=65", "age=74"))
plot(res1, xnewdata = xpred, xtfnewdata = xpred, tgrid = tgrid, cex = 2);
dev.off()

########### frailty maps ##################
nwengland=read.bnd(system.file("otherdata/nwengland.bnd", package="spBayesSurv"));
frail0=-(rowMeans(res1$v));
frail = frail0[as.integer(names(nwengland))];
values = cbind(as.integer(names(nwengland)), frail) 
op <- par(no.readonly = TRUE)
pdf(file ="Leukemia-GAFT-CAR-Map.pdf", paper="special", width=8, height=8)
par(mar=c(3,0,0,0))
plotmap(nwengland, x=values, col=(gray.colors(10,0.3,1))[10:1], 
        pos = "bottomleft",width = 0.5, height = 0.04)
dev.off()


#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the PH copula model 
## in Section 4.3.1 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data("LeukSurv")
d <- LeukSurv[order(LeukSurv$district), ]
head(d) 

##-------------PH-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior=list(M=20, nknots=100, nblock=1043);
ptm<-proc.time()
res1 = spCopulaCoxph(formula=Surv(time,cens)~age+sex+wbc+tpi,data=d,
                     mcmc=mcmc,prior=prior,
                     Coordinates=cbind(d$xcoord,d$ycoord));
sfit1=summary(res1); sfit1
systime1=proc.time()-ptm; systime1;

save.image("Result_PH_copula.RData")

############# Trace plots #############
pdf(file ="Leukemia-PH-copula-trace.pdf", paper="special", width=10, height=5)
par(mfrow=c(3,2));
par(cex=1,mar=c(2.5,4.1,1,1))
traceplot(mcmc(res1$beta[1,]), xlab="", main="age")
traceplot(mcmc(res1$beta[2,]), xlab="", main="sex")
traceplot(mcmc(res1$beta[3,]), xlab="", main="wbc")
traceplot(mcmc(res1$beta[4,]), xlab="", main="tpi")
traceplot(mcmc(res1$theta[1,]), xlab="", main="partial sill")
traceplot(mcmc(res1$theta[2,]), xlab="", main="range")
dev.off()

############# Survival Curves #############
pdf(file ="Leukemia-PH-copula-age.pdf", paper="special", width=8, height=8);
par(mfrow=c(1, 1)) 
tgrid = seq(0.1, 5000, length.out=300);
# age effect
xpred <- data.frame(age = c(49, 65, 74),
                    sex = c(0, 0, 0),
                    wbc = c(38.59, 38.59, 38.59),
                    tpi = c(0.3398, 0.3398, 0.3398),
                    row.names = c("age=49", "age=65", "age=74"))
plot(res1, xnewdata = xpred, tgrid = tgrid, cex = 2);
dev.off()

########### frailty maps ##################
nwengland=read.bnd(system.file("otherdata/nwengland.bnd", package="spBayesSurv"));
frail= round((rowMeans(res1$Zpred)),3); nclust=5;
frail.cluster = cut(frail, breaks = nclust);
frail.names = names(table(frail.cluster))
rbPal <- colorRampPalette(c('red','blue'))
frail.colors=rbPal(nclust)[as.numeric(frail.cluster)]
pdf(file ="Leukemia-PH-copula-Map.pdf", paper="special", width=8, height=8)
par(mar=c(3,0,0,0))
plot(nwengland)
points(cbind(d$xcoord,d$ycoord), col=frail.colors)
legend("topright",title="z values",legend=frail.names,
       col=rbPal(nclust),pch=20, cex=1.7)
dev.off()


#####-----------------------------------------------------######
## Analysis of the Leukemia Data
## Below is R code for fitting the LDDPM copula model 
## in Section 4.3.2 of the paper
#####-----------------------------------------------------######
rm(list=ls())
library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")

############## Read data ##############
data(LeukSurv);
attach(LeukSurv);
d = LeukSurv[order(district),]; n = nrow(d); detach(LeukSurv);
head(d); 

##-------------LDDPM-------------------##
set.seed(1)
# MCMC parameters
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000)
prior=list(N=10, nknots=100, nblock=1043);
ptm<-proc.time()
res1 = spCopulaDDP(formula=Surv(time,cens)~age+sex+wbc+tpi,data=d,
                   mcmc=mcmc,prior=prior,
                   Coordinates=cbind(d$xcoord,d$ycoord));
sum(log(res1$cpo)); ## LPML
systime1=proc.time()-ptm; systime1;

save.image("Result_LDDPM_copula.RData")

############# Trace plots #############
pdf(file ="Leukemia-LDDPM-copula-trace.pdf", paper="special", width=8, height=8)
par(mfrow=c(2,1));
par(cex=1,mar=c(2.5,4.1,1,1))
traceplot(mcmc(res1$theta[1,]), xlab="", main="partial sill")
traceplot(mcmc(res1$theta[2,]), xlab="", main="range")
dev.off()

############# Survival Curves #############
pdf(file ="Leukemia-LDDPM-copula-age.pdf", paper="special", width=8, height=8);
par(mfrow=c(1, 1)) 
tgrid = seq(0.1,5000,length.out=300);
# age effect
xpred <- data.frame(age = c(49, 65, 74),
                    sex = c(0, 0, 0),
                    wbc = c(38.59, 38.59, 38.59),
                    tpi = c(0.3398, 0.3398, 0.3398),
                    row.names = c("age=49", "age=65", "age=74"))
plot(res1, xnewdata = xpred, tgrid = tgrid, cex = 2);
dev.off()

########### frailty maps ##################
nwengland=read.bnd(system.file("otherdata/nwengland.bnd", package="spBayesSurv"));
frail= round((rowMeans(res1$Zpred)),3); nclust=5;
frail.cluster = cut(frail, breaks = nclust);
frail.names = names(table(frail.cluster))
rbPal <- colorRampPalette(c('red','blue'))
frail.colors=rbPal(nclust)[as.numeric(frail.cluster)]
pdf(file ="Leukemia-LDDPM-copula-Map.pdf", paper="special", width=8, height=8)
par(mar=c(3,0,0,0))
plot(nwengland)
points(cbind(d$xcoord,d$ycoord), col=frail.colors)
legend("topright",title="z values",legend=frail.names,
       col=rbPal(nclust),pch=20, cex=1.7)
dev.off()


