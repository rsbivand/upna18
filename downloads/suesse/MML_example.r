# make sure all files are in same directory, i.e. ARmodels.r, ARmodels.sparse.r and ARmodels.sparse.add.fcts.r

# House data set

#required R libraries  spdep, Matrix, numDeriv

# select systematic sample of size ns=1015 or ns=5072 ??
#choice=1: ns=1015
#choice=2: ns=5072

choice<-1

model<-"SAR"
# SAR means SAM
# SARerr means SEM

if(model=="SARerr"){ARmod<-2}
if(model=="SAR"){ARmod<-1}


library(spdep)
data(house)
hform <- formula(log(price) ~ age + I(age^2) + I(age^3) +
                   + log(lotsize) + rooms + log(TLA) + beds + syear)
hmm0 <- model.matrix(hform, data = house)
hlw <- nb2listw(LO_nb)
hsn <- listw2sn(hlw)
W<- as(hlw, "CsparseMatrix")
X<-hmm0
n<-dim(W)[1]
Y<-log(house$price)
if(choice==1){
s<- seq(1,n,by=5)  # sample
}else{
  s<- seq(1,n,by=25) #sample
}
ns<-length(s)
u<-setdiff(1:n,s)  # non-sample units


# full data set, all n=25357 observations
# using spdep to fit full data set without missing observations


# SAM
system.time(
  house.SAM <- lagsarlm(hform, data = house, listw = hlw,method = "Matrix")
)
house.SAM$rho # 0.5228141

# alternative the SEM 
system.time(
  house.SEM <- errorsarlm(hform, data = house, listw = hlw,method = "Matrix")
)
house.SEM$lambda # 0.6194028


WtW<-t(W)%*%W
WplusWt<-W+t(W)

source("ARmodels.sparse.r")




# extract sample information and re-order X and W, such that first ns units refer to sample (observed)
# and remaing nu=n-ns units refer to non-sample (missing)
  
  Ys<-Y[s]
  X1<-X[c(s,u),]
  W1<-W[c(s,u),c(s,u)] 
  
  
# fit SAM model (SAR)
# function:  AR.pop.exact.sparse  
# arguments 
# Ys vector of ns observed units of response variable
# W1 n*n contiguity matrix, where first ns units refer to Ys
# X1 n*p design matrix, where first ns units refer to Ys
# model="SAR" (SAM) or model="SARerr" (SEM)
# se=TRUE or FALSE: calculate standard errors yes or no, default TRUE
# approx.se= TRUE or FALSE (approximate standard errors), default=TRUE
# approx=TRUE or FALSE (use approximate method or not), default =FALSE
# rho.limits limits of rho, default=c(-1,1)
# rho0 starting value of rho, default=NULL  
  
# rho.limits can also be found with
#  eigen.values<-eigen(W)$values # takes too much time and memory
#  rho.limits<-c(1/min(eigen.values),1/max(eigen.values))
#  -1 and 1 OK, upper bound +1 known, as W is row-normalised
  

# by defult in does optimisation for positive rho and then seperately for negative values and then
# chooses the rho that yields a larger log-likelihood
# this is mostly unnecessary (and not done in Suesse 2017) as it increases the time by a factor of 2, but it ensures that
# global maxima is not missed
  
  
# SAM
    
 lag.exact<-try(AR.pop.exact.sparse(Ys,W1,X1,model="SAR",se=T,approx.se=F,approx=F,rho.limits=c(-1,1),tol=.Machine$double.eps));
 
 lag.exact$L  # log-lik
 
 # estimates
 
 lag.exact$rho;  # rho
 lag.exact$beta; #beta
 lag.exact$sigma2; #sigma^2

 
 # standard errors
 lag.exact$rho.se  # rho
 lag.exact$beta.se #beta
 lag.exact$sigma2.se #sigma^2
 
 # SEM
 
 lag.exact<-try(AR.pop.exact.sparse(Ys,W1,X1,model="SARerr",se=T,approx.se=F,approx=F,rho.limits=c(-1,1),tol=.Machine$double.eps));
 
 lag.exact$L  # log-lik
 
 # estimates
 lag.exact$rho;  # rho
 lag.exact$beta; #beta
 lag.exact$sigma2; #sigma^2
 
 # standard errors
 lag.exact$rho.se  # rho
 lag.exact$beta.se #beta
 lag.exact$sigma2.se #sigma^2
  
  

# approximate method based on Taylor series approximation of order K
 
# for first use W-matrices are calculated, for second use these can be used without recalculating them
# by obtainng the Wlists entry, i.e. Wlists=lag.approx$Wlists

lag.approx<-try(AR.pop.approx.sparse(Ys,W1,X1,model="SARerr",se=TRUE,tol=1e-6,K=10,Wlists=NULL,rho.limits=c(-1,1)))

lag.approx$L  # approximated log-lik

# estimates

lag.approx$L
lag.approx$rho
lag.approx$beta

# standard errors

lag.approx$sigma2
lag.approx$rho.se
lag.approx$beta.se
lag.approx$sigma2.se

# obtain exact log-lik for approximate method
AR.exact.log.lik.fct(par=lag.approx$rho,Ys,W1,X1,WplusWt[c(s,u),c(s,u)],WtW[c(s,u),c(s,u)],ARmod=ARmod,pr=TRUE)


