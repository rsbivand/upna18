#######################################################################################################################  
# BEGIN REML/ML 
########################################################################################################################
# Code to fit SLM-M and SEM-M according to Suesse (2017) in Computational Statistics

library(spdep)
source("SAR.M.sparse.r")


  data(house)
  hform <- formula(log(price) ~ age + I(age^2) + I(age^3) +
                     + log(lotsize) + rooms + log(TLA) + beds + syear)
  hmm0 <- model.matrix(hform, data = house)
  Y<-log(house$price)
  hlw <- nb2listw(LO_nb)
  hsn <- listw2sn(hlw)
  W<- as(hlw, "CsparseMatrix")
  X<-hmm0
  n<-dim(X)[1]
  rownames(W)<-colnames(W)<-1:n
  rownames(X)<-1:n
  W<- as(W, "CsparseMatrix")
  
  

  



# arguments
# Y response vector
# W contiguity matrix
# X design matrix
# model: either "SAR" (SLM-M) or "SARerr" (SEM-M)
# REML: TRUE or FALSE, if FALSE then ML
# se: TRUE or FALSE, if FALSE no standard errors
# rho0 starting value
# omega.eps0: starting value of sigma2.eps
# alternatively omega.Y0: starting value of sigma2.y

  model<-  "SAR" 
  
  # Fit model on single core
  
  system.time(
  resREML<-try(AR.error(Y,W,X,model=model,se=TRUE,REML=FALSE,rho0=0.5,omega.eps0=NULL,omega.Y0=NULL))
  )

  beta.hat<-resREML$beta  # beta estimates
  rho.hat<-  resREML$rho  # rho estimates
  sigma2.eps.hat<- resREML$omega.eps # measurement error variance
  rho.hat.se<-resREML$rho.se  # s.e. of rho
  beta.hat.se<- resREML$beta.se # s.e. of beta
  omega.eps.hat.se<-resREML$omega.eps.se # s.e. of measurement error variance estimate
  sigma2.y.hat<-resREML$omega.Y # SAR model variance
  sigma2.y.hat.se<-resREML$omega.Y.se # s.e. of SAR model variance estimate 
  L<-resREML$L  # log-likelihood
  
#######################################################################################################################  
  # END REML/ML 
########################################################################################################################  

 # compare with INLA for SAR (SLM-M)
 # to install INLA use:  install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
  
  library(INLA)
  
  e_extr <- c(lextrW(hlw))
  rho.max <- 1/e_extr[2]
  rho.min <- 1/e_extr[1]
  rho <- mean(c(rho.min, rho.max))
  n <- length(Y)
  house$idx <- 1:n
  zero.variance <- list(prec=list(initial = 25, fixed=TRUE)) 
  
  betaprec <- .0001
  Q.beta = Diagonal(n=ncol(X), betaprec)
  
  ## some Priors on the hyperparameters
  ## to be modified
  hyper.slm = list(
    prec = list(
      prior = "loggamma",
      param = c(0.01, 0.01)),
    rho = list(
      initial=0,
      prior = "logitbeta",
      param = c(1,1)))
  
  args.slm<- list(
    rho.min = rho.min,
    rho.max = rho.max,
    W=W,
    X=X,
    Q.beta=Q.beta)
  
  
  ## Fit model (multiple cores)
  system.time(
    slmm1<-inla(update(hform, . ~ . + f(idx, model="slm", args.slm=args.slm, hyper=hyper.slm)), data=as(house, "data.frame"), family="gaussian", control.family = list(hyper=zero.variance), control.compute=list(dic=TRUE, cpo=TRUE))
    )
  
  
  

  summary(slmm1)
  ## Summary of the coefficients (at the end of the vector of random effects)
  slmm1$summary.random$idx[n+1:ncol(X),]
  
  ## Re-scale rho to real scale
  rhomarg <- inla.tmarginal(function(x){rho.min+x*(rho.max-rho.min)},
                            slmm1$marginals.hyperpar[[2]])
  inla.zmarginal(rhomarg)
  
  
