
###############
# simulate data

library("spatsurv")
set.seed(11)
n <- 300
DIST <- weibullHaz()
OMEGA <- c(0.5, 2)
COVMODEL <- ExponentialCovFct()
COVPARS <- c(0.7, 0.1)


dat <- simsurv(X = cbind(age = runif(n, 5, 50), sex = rbinom(n, 1, 0.5),
                         cancer = rbinom(n, 1, 0.2)),
               beta = c(0.0296, 0.0261, 0.035), dist = DIST, omega = OMEGA,
               cov.parameters = COVPARS, cov.model = COVMODEL,
               mcmc.control = mcmcpars(nits = 110000, burn = 10000, thin = 100))


plot(dat$T[, 178], type = "s")
plot(apply(dat$T, 2, function(x){acf(x, plot = FALSE)$acf[2]}),
     xlab = "Subject Index", ylab = "Lag 1 autocorrelation")


survtimes <- dat$survtimes
censtimes <- rexp(n, 1 / mean(survtimes))
survdat <- gencens(survtimes, censtimes)


library("sp")
coords <- dat$coords
X <- as.data.frame(dat$X)
spatdat <- SpatialPointsDataFrame(coords, data = as.data.frame(X))
spatdat$ss <- survdat



##########
# plotting

plotsurv(spp = spatdat, ss = spatdat$ss, maxcex = 3)




###################
# specifying priors

betaprior <- betapriorGauss(mean = 0, sd = 10)
omegaprior <- omegapriorGauss(mean = 0, sd = 10)
etaprior <- etapriorGauss(mean = c(0, -2), sd = 0.5)


priors <- mcmcPriors(betaprior = betaprior, omegaprior = omegaprior,
                     etaprior = etaprior, call = indepGaussianprior,
                     derivative = derivindepGaussianprior)



############################
# Running the MCMC Algorithm

mod <- survspat(formula = ss ~ age + sex + cancer, data = spatdat,
                dist = DIST, cov.model = COVMODEL,
                mcmc.control = mcmcpars(nits = 500000, burn = 10000, thin = 490),
                priors = priors)



##################                        
# MCMC Diagnostics

library("mcmcplots")
mcmcplot(cbind(mod$betasamp, mod$omegasamp, mod$etasamp))


frailtylag1(mod)
plot(mod$tarrec, ylab = "log-Posterior", type = "s")


mod

#################
# Post-Processing   

# Prior and posterior
priorposterior(mod)

# Baseline Hazard
baselinehazard(mod)

# BSpatial covariance
posteriorcov(mod)

# Posterior Hazard
predict(mod, type = "hazard", indx = 260)

# Posterior Survival
predict(mod, type = "survival", indx = 260)

# Posterior Density
predict(mod, type = "density", indx = 260)

# DensityQuantiles
dq <- predict(mod, type = "densityquantile")


hazardexceedance <- function(threshold, direction = "upper"){
  fun <- function(beta, omega, eta, Y){
    EY <- exp(Y)
    d <- length(Y)
    len <- length(threshold)
    A <- matrix(NA, len, d)
    
    for(i in 1:len){
      if(direction == "upper"){
        A[i, ] <- as.numeric(EY > threshold[i])
      }
      else{
        A[i, ] <- as.numeric(EY < threshold[i])
      }
    }
    return(A)
  }
  attr(fun, "threshold") <- threshold
  attr(fun, "direction") <- direction
  return(fun)
}

# Monte Carlo Expectations
func <- hazardexceedance(c(1.5, 2, 3))
mchaz <- MCE(mod, func)
mchaz[, 1:10]





##################################################################
# London Fires Analysis
##################################################################

library("spatsurv")
library("spatstat")
library("sp")
library("survival")
library("rgeos")
library("rgdal")
library("leaflet")
library("fields")
library("raster")
library("lubridate")

data("fstimes")


for(i in 1:4){
  fstimes@data[paste("s", i, sep="")] <-
    sin(i * 2 * pi * fstimes$timenumeric / 24)
  fstimes@data[paste("c", i, sep="")] <-
    cos(i * 2 * pi * fstimes$timenumeric / 24)
}


data("fs")
fs <- fs[fs$Description=="Fire Station", ]
fscoords <- cbind(fs$Easting,fs$Northing)
chull <- convexhull.xy(rbind(coordinates(fstimes), fscoords))
win <- expandwinPerfect(chull, TRUE, 1.2)
fsppp <- ppp(x=fscoords[, 1], fscoords[, 2], window=win)
fsintens <- density.ppp(fsppp)
fsintens <- raster(fsintens)
proj4string(fsintens) <- CRS("+init=epsg:27700")
fstimes$fsintens <- raster::extract(fsintens, fstimes) * 100000000


FORM <- S ~ fsintens + s1 + c1 + s2 + c2 + s3 + c3 + s4 + c4
betaprior <- betapriorGauss(mean = 0, sd = 100)
omegaprior <- omegapriorGauss(mean = 0, sd = 10)
etaprior <- etapriorGauss(mean = log(c(1, 1000)), sd = 0.5)
priors <- mcmcPriors(betaprior = betaprior, omegaprior = omegaprior,
                     etaprior = etaprior, call = indepGaussianprior,
                     derivative = derivindepGaussianprior)


mod <- survspat(formula = FORM, data = fstimes,
                dist = BsplineHaz(fstimes$S[, 1]), cov.model = ExponentialCovFct(),
                mcmc.control = mcmcpars(nits = 500000, burn = 10000, thin = 490),
                priors = priors,
                control = inference.control(gridded = TRUE, cellwidth = 1000))


tseq <- seq(0, 24, length.out = 100)
getpred <- function(par){
  mat <- rbind( sin(2 * pi * tseq / 24),
                cos(2 * pi * tseq / 24),
                sin(2 * pi * tseq * 2 / 24),
                cos(2 * pi * tseq * 2 / 24),
                sin(2 * pi * tseq * 3 / 24),
                cos(2 * pi * tseq * 3 / 24),
                sin(2 * pi * tseq * 4 / 24),
                cos(2 * pi * tseq * 4 / 24))
  return(exp(colSums(par*mat)))
}
timetrend <- apply(mod$betasamp[, -1], 1, getpred)
qts <- t(apply(timetrend, 1, quantile, probs=c(0.025, 0.5, 0.975)))
matplot(tseq, qts, type = "l", col = c("purple", "black", "blue"),
        lty = c("dashed", "solid", "dashed"), xlab = "Time of Day",
        ylab = "Relative Risk")
legend("topleft", lty = c("dashed", "solid", "dashed"),
       col = rev(c("purple", "black", "blue")),
       legend = rev(c(0.025, 0.5, 0.975)))


func <- hazardexceedance(c(0.9, 0.8, 0.7), direction = "lower")
mchaz <- MCE(mod, func)


gr <- getGrid(mod)
gr$exceed <- mchaz[1, ]
m <- spplot1(gr, "exceed", useLeaflet = TRUE,
             palette = rev(brewer.pal(5, "RdBu")))
m


mod <- survspat(formula = FORM, data = cancerdata, dist = DIST,
                cov.model = COVMODEL, mcmc.control = mcmcpars(nits = 500000,
                                                              burn = 10000, thin = 490), priors = priors, shape = stateshp,
                ids = list(shpid = "FIPS", dataid = "STCOUNTY"))


ExponentialCovFct <- function(){
  ans <- list()
  ans$npar <- 2
  ans$parnames <- c("sigma", "phi")
  ans$itrans <- exp
  ans$trans <- log
  ans$eval <- function(u, pars){
    ans<- pars[1] ^ 2 * exp(-u / pars[2])
    return(ans)
  }
  class(ans) <- c("covmodel", "fromUserFunction")
  return(ans)
}


weibullHaz <- function(){
  
  flist <- list()
  
  flist$distinfo <- function(){
    retlist <- list()
    retlist$npars <- 2
    retlist$parnames <- c("alpha", "lambda")
    retlist$trans <- log
    retlist$itrans <- exp
    retlist$jacobian <- exp
    retlist$hessian <- list(exp, exp)
    return(retlist)
  }
  
  flist$basehazard <- function(pars){
    fun <- function(t){
      return(pars[2] * pars[1] * t ^ (pars[1] - 1))
    }
    return(fun)
  }
  
  flist$gradbasehazard <- function(pars){
    fun <- function(t){
      return(t ^ (pars[1] - 1) *
               cbind(pars[2] * (1 + pars[1] * log(t)), pars[1]))
    }
    return(fun)
  }
  
  flist$hessbasehazard <- function(pars){
    funfun <- function(t, pars){
      m <- matrix(0, 2, 2)
      m[1, 2] <- m[2, 1] <- t ^ (pars[1] - 1) *
        (1 + pars[1] * log(t))
      m[1, 1] <- pars[2] * t ^ (pars[1] - 1) *
        log(t) * (2 + pars[1] * log(t))
      return(m)
    }
    
    fun <- function(t){
      return(lapply(t, funfun, pars = pars))
    }
    return(fun)
  }
  
  flist$cumbasehazard <- function(pars){
    fun <- function(t){
      return(pars[2] * t ^ (pars[1]))
    }
    return(fun)
  }
  
  flist$gradcumbasehazard <- function(pars){
    fun <- function(t){
      return(t ^ (pars[1]) * cbind(pars[2] * log(t), 1))
    }
    return(fun)
  }
  
  flist$hesscumbasehazard <- function(pars){
    funfun <- function(t, pars){
      m <- matrix(0, 2, 2)
      other <- log(t) * t ^ pars[1]
      m[1, 2] <- m[2, 1] <- other
      m[1, 1] <- pars[2] * other * log(t)
      return(m)
    }
    
    fun <- function(t){
      return(lapply(t, funfun, pars = pars))
    }
    return(fun)
  }
  
  flist$densityquantile <- function(pars, other){
    fun <- function(probs){
      return((-log(1 - probs) /
                (pars[2] * other$expXbetaplusY)) ^ (1 / pars[1]))
    }
    return(fun)
  }
  
  class(flist) <- c("basehazardspec", "list")
  return(flist)
}






##################################################################
# Simulation Study
##################################################################

set.seed(1)   

n <- 500

omega_list <- list()
omega_list[[1]] <- 0.25
omega_list[[2]] <- 0.5
omega_list[[3]] <- 0.75
omega_list[[4]] <- 1
omega_list[[5]] <- 1.25
omega_list[[6]] <- 1.5
omega_list[[7]] <- 1.75
omega_list[[8]] <- 2
omega_list[[9]] <- c(0.25, 1)
omega_list[[10]] <- c(0.5, 1)
omega_list[[11]] <- c(0.75, 1)
omega_list[[12]] <- c(1, 1)
omega_list[[13]] <- c(0.25, 2)
omega_list[[14]] <- c(0.5, 2)
omega_list[[15]] <- c(0.75, 2)
omega_list[[16]] <- c(1, 2)

cov_pars_list <- list()
cov_pars_list[[1]] <- c(0.5, 0.02)
cov_pars_list[[2]] <- c(1, 0.02)
cov_pars_list[[3]] <- c(1.5, 0.02)
cov_pars_list[[4]] <- c(2, 0.02)
cov_pars_list[[5]] <- c(0.5, 0.05)
cov_pars_list[[6]] <- c(1, 0.05)
cov_pars_list[[7]] <- c(1.5, 0.05)
cov_pars_list[[8]] <- c(2, 0.05)
cov_pars_list[[9]] <- c(0.5, 0.1)
cov_pars_list[[10]] <- c(1, 0.1)
cov_pars_list[[11]] <- c(1.5, 0.1)
cov_pars_list[[12]] <- c(2, 0.1)
cov_pars_list[[13]] <- c(0.5, 0.15)
cov_pars_list[[14]] <- c(1, 0.15)
cov_pars_list[[15]] <- c(1.5, 0.15)
cov_pars_list[[16]] <- c(2, 0.15)

beta_list <- list()
beta_list[[1]] <- c(0.03, 0.03)
beta_list[[2]] <- c(0.06, 0.03)
beta_list[[3]] <- c(0.09, 0.03)
beta_list[[4]] <- c(0.12, 0.03)
beta_list[[5]] <- c(0.03, 0.09)
beta_list[[6]] <- c(0.06, 0.09)
beta_list[[7]] <- c(0.09, 0.09)
beta_list[[8]] <- c(0.12, 0.09)
beta_list[[9]] <- c(0.03, 0.15)
beta_list[[10]] <- c(0.06, 0.15)
beta_list[[11]] <- c(0.09, 0.15)
beta_list[[12]] <- c(0.12, 0.15)
beta_list[[13]] <- c(0.03, 0.21)
beta_list[[14]] <- c(0.06, 0.21)
beta_list[[15]] <- c(0.09, 0.21)
beta_list[[16]] <- c(0.12, 0.21)



for(i in 1:16){
  print(i)
  
  haz <- exponentialHaz()
  if(i>8){
    haz <- weibullHaz()
  }
  
  dat <- simsurv(X = cbind( age = runif(n, 5, 50), sex = rbinom(n, 1, 0.5)), 
                 beta = beta_list[[i]], 
                 dist = haz, 
                 omega = omega_list[[i]], 
                 cov.parameters = cov_pars_list[[i]], 
                 cov.model = ExponentialCovFct(), 
                 mcmc.control = mcmcpars(nits = 100000, burn = 10000, thin = 1000))  
  
  survtimes <- dat$survtimes
  censtimes <- rexp(n, 1/mean(survtimes))                                    
  survdat <- gencens(survtimes, censtimes) 
  
  coords <- dat$coords
  X <- as.data.frame(dat$X) # covariates 
  spatdat <- SpatialPointsDataFrame(coords, data = as.data.frame(X))
  spatdat$ss <- survdat
  
  assign(paste("spatdat", i, sep = ""), spatdat)
  assign(paste("dat", i, sep = ""), dat)
  assign(paste("censtimes", i, sep = ""), censtimes)
}


betaprior <- betapriorGauss(mean = 0, sd = 10)
omegaprior <- omegapriorGauss(mean = 0, sd = 10)
etaprior <- etapriorGauss(mean = c(0, -3), sd = c(0.5, 0.4))

priors <- mcmcPriors(betaprior = betaprior, 
                     omegaprior = omegaprior, 
                     etaprior = etaprior, 
                     call = indepGaussianprior, 
                     derivative = derivindepGaussianprior)

for(SEED in c(1, 2, 3)){
  set.seed(SEED)
  
  for(i in 1:16){
    
    print(paste("Analysing dataset", i, "..."))
    
    haz <- exponentialHaz()
    if(i>8){
      haz <- weibullHaz()
    }
    
    mod <- survspat(formula = ss ~ age + sex, 
                    data = get(paste("spatdat", i, sep="")), 
                    dist = haz, 
                    cov.model = ExponentialCovFct(), 
                    mcmc.control = mcmcpars(nits = 1000000, 
                                            burn = 10000, 
                                            thin = 990), 
                    control = inference.control(gridded = TRUE, cellwidth = 0.034), 
                    priors = priors)
  }
}
