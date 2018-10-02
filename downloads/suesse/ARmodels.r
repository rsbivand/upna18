#rm(list=ls())
#setwd("E:/Thomas/UoW/R/Ray")
#source("ARmodels.R")
library(Matrix)


log.lik.AR<-function(rho,lambda,tyMy,ytMyL,yLtMyL,n,grad=TRUE){
s2<-tyMy-2*rho*ytMyL+(rho^2)*yLtMyL
res<- -2/n*Re(sum(log(1-rho*lambda)))+log(s2)
if(grad){
attr(res, "gradient") <-  2/n*Re(sum(lambda/(1-rho*lambda)))+2*(rho*yLtMyL - ytMyL)/s2
}
res
}#end function

log.lik.AR.L<-function(rho,lambda,tyMy,ytMyL,yLtMyL,n){log.lik.AR(rho,lambda,tyMy,ytMyL,yLtMyL,n,grad=FALSE)}
log.lik.AR.gr<-function(rho,lambda,tyMy,ytMyL,yLtMyL,n){attr(log.lik.AR(rho,lambda,tyMy,ytMyL,yLtMyL,n,grad=TRUE),"gradient")}


################################################################################
### AR model  ##################################################################
################################################################################

AR<-function(W,Y,X,rho0,lambda){
n<-length(Y)
I<-diag(n)
M1<-solve(t(X)%*%X)%*%t(X)
M<-I-X%*%M1
YL<-W%*%Y
tyMy<-t(Y)%*%M%*%Y
ytMyL<-t(Y)%*%M%*%YL
yLtMyL<-t(YL)%*%M%*%YL

abs.lambda<-abs(lambda)
#lambda.min<-min(abs.lambda)
lambda.max<-max(abs.lambda)
#find ML estimate for rho
#for(i in 1:1e3){
#r<-optimize(f = log.lik.AR, interval = c(-1/lambda.max,1/lambda.max),lambda,tyMy,ytMyL,yLtMyL,n,maximum = FALSE)
#}

#rho.hat<-r$minimum
#for(i in 1:1e3){
r<-optim(par=rho0,fn=log.lik.AR.L , gr =log.lik.AR.gr ,lambda,tyMy,ytMyL,yLtMyL,n,method = c("L-BFGS-B"),lower = -1/lambda.max+0.001, upper = 1/lambda.max-0.001)
     # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
#}
#rho.hat<-r$par
#for(i in 1:1e3){
#r<-nlm(f=log.lik.AR, p=rho0,lambda,tyMy,ytMyL,yLtMyL,n, hessian = FALSE)
#}
if(r$convergence==0){
rho.hat<-r$par
counts<-r$counts
}else{
#return error
return(NULL)
}

p<-dim(X)[2]
#now compute beta and sigma
A<-I-rho.hat*W
A.inv<-solve(A)
z<-A%*%Y
beta.hat<-M1%*%z
N<-  n
sigma.hat<-c(t(z)%*%M%*%z)/N




#now compute variance matrix of parameter estimates
B<-A.inv%*%W
muL<-B%*%X%*%beta.hat
alpha<-Re(sum(lambda^2/(1-rho.hat*lambda)^2))
#first beta, rho, sigma
V11<- sigma.hat*(t(X)%*%X)   #beta, beta
V12<- sigma.hat*t(X)%*%muL  #beta,rho
V13<- rep(0,p)                    #beta, sigma
V22<-sigma.hat^2*tr.prod(t(B),B) + sigma.hat*t(muL)%*%muL - alpha*sigma.hat^2   #rho, rho
V23<-sigma.hat*tr(B)     #rho, sigma
V33<-n/2



V<-rbind(cbind(V11,V12,V13),cbind(t(V12),V22,V23),cbind(t(V13),t(V23),V33))
acvm<-sigma.hat^2*solve(V)
lnlik.model<- - n/2*log(2*pi) - N/2 - n/2*log(sigma.hat)

par.names<-c(paste("beta",1:p),"rho","sigma2")
colnames(acvm)<-rownames(acvm)<-par.names

sigma<-sqrt(sigma.hat)
beta<-beta.hat
rho1<-rho.hat

#gradient
grad<-log.lik.AR.gr(rho1,lambda,tyMy,ytMyL,yLtMyL,n)


return(list(beta=beta,sigma=sigma,rho1=rho1,lnlik.model=lnlik.model,acvm=acvm,grad=grad,counts=counts))

}#end function  AR

################################################################################
### AR-error model  ############################################################
################################################################################

log.lik.ARerr<-function(rho,lambda,Y,X,W,n,grad=TRUE,REML=REML){
p<-dim(X)[2]
I<-diag(n)
A<-I-rho*W
Vinv<-t(A)%*%A  #V<-Ainv%*%t(Ainv), i.e. Vinv<-t(A)%*%A
XtVinv<-t(X)%*%Vinv
M<-X%*%solve(XtVinv%*%X)%*%XtVinv
M.hat<-Vinv%*%(I - M)
YtM.hatY<-t(Y)%*%M.hat%*%Y
#that's plus L
#res<-  n/2*log(2*pi)    - Re(sum(log(1-rho*lambda))) + 0.5*n*log(YtM.hatY) + 0.5*n - 0.5*n*log(n) 
if(REML){
AX<-A%*%X
H<-t(AX)%*%AX
N<-n-p
res<- n/2*log(2*pi)    - Re(sum(log(1-rho*lambda))) + 0.5*N*log(YtM.hatY) + 0.5*N - 0.5*n*log(n) +  0.5*log(det(H))
}else{
res<-  n/2*log(2*pi)    - Re(sum(log(1-rho*lambda))) + 0.5*n*log(YtM.hatY) + 0.5*n - 0.5*n*log(n) 
} 

#Re(sum(log(1-rho*lambda)))

#print(Re(sum(log(1-rho*lambda))))
#print(log(YtM.hatY))

#print(res)


if(grad){
dVinvdrho<-2*rho*t(W)%*%W - W - t(W)
dM.hatdrho<-  dVinvdrho%*%(I-M) + t(M)%*%dVinvdrho%*%(M-I)
YtdM.hatdrhoY<-t(Y)%*%dM.hatdrho%*%Y
hilf<-  Re(sum(lambda/(1-rho*lambda)))+ n/2*YtdM.hatdrhoY/YtM.hatY
if(REML){
Hinv<-solve(H)
WX<-W%*%X
WXtAX<-t(WX)%*%AX
dHdrho<- -  WXtAX - t(WXtAX)
                # derivative of  0.5*N*log(YtM.hatY)   derivative of 0.5*log(det(H))
hilf<-hilf - (p/2)*YtdM.hatdrhoY/YtM.hatY             + 0.5*tr.prod(Hinv,dHdrho) 
}
attr(res, "gradient") <- hilf 
}#end grad
res 
}#end function

log.lik.ARerr.L<-function(rho,lambda,Y,X,W,n,REML=REML){log.lik.ARerr(rho,lambda,Y,X,W,n,grad=FALSE,REML=REML)}
log.lik.ARerr.gr<-function(rho,lambda,Y,X,W,n,REML=REML){attr(log.lik.ARerr(rho,lambda,Y,X,W,n,grad=TRUE,REML=REML),"gradient")}



ARerr<-function(W,Y,X,rho0,lambda,REML=TRUE){
n<-length(Y)
I<-diag(n)
#find ML estimate for rho
#test log.lik.ARerr(rho0,lambda,Y,X,W,n)
abs.lambda<-abs(lambda)
lambda.max<-max(abs.lambda)
#find ML estimate for rho

# test log.lik.ARerr(rho0,lambda,Y,X,W,n)
#for(i in 1:1e1){
#r<-optimize(f = log.lik.ARerr, interval = c(-1/lambda.max,1/lambda.max),lambda,Y,X,W,n,maximum = FALSE)
#}

#rho.hat<-r$minimum
#for(i in 1:1e2){
r<-optim(par=rho0,log.lik.ARerr.L, gr = log.lik.ARerr.gr,lambda,Y,X,W,n,REML=REML,method = c("L-BFGS-B"),
lower = -1/lambda.max+0.001, upper = 1/lambda.max-0.001)
#}     # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
#for(i in 1:1e1){
#r<-nlm(f=log.lik.ARerr, p=rho0,lambda,Y,X,W,n, hessian = FALSE)
#}

if(r$convergence==0){
rho.hat<-r$par
counts<-r$counts
}else{
#return error
return(NULL)
}              

p<-dim(X)[2]
#now compute beta and sigma
A<-I-rho.hat*W
A.inv<-solve(A)
Vinv<-t(A)%*%A  #V<-Ainv%*%t(Ainv), i.e. Vinv<-t(A)%*%A
XtVinv<-t(X)%*%Vinv
XtVinvX<-XtVinv%*%X
XtVinvX.inv<-solve(XtVinvX)
XtVinvX.inv.XtVinv<-XtVinvX.inv%*%XtVinv
M.hat<-Vinv%*%(I - X%*%XtVinvX.inv.XtVinv)
beta.hat <- XtVinvX.inv.XtVinv%*%Y
if(REML){N<-n-p}else{N<-n}
sigma.hat<-c(t(Y)%*%M.hat%*%Y)/N


#now compute variance matrix of parameter estimates
B<- A.inv%*%W
AX<-A%*%X
muL<-B%*%X%*%beta.hat
alpha<-Re(sum(lambda^2/(1-rho.hat*lambda)^2))
#first beta, rho, sigma
V11<- sigma.hat*(t(AX)%*%AX)   #beta, beta
V12<- rep(0,p)                       #beta,rho
V13<- rep(0,p)                    #beta, sigma
V22<-sigma.hat^2*tr.prod(t(B),B)  - alpha*sigma.hat^2   #rho, rho
V23<-sigma.hat*tr(B)     #rho, sigma
V33<-n/2                 #sigma,sigma

V<-rbind(cbind(V11,V12,V13),cbind(t(V12),V22,V23),cbind(t(V13),t(V23),V33))/sigma.hat^2


if(REML){
#add terms for REML, actually only for rho and for omega
WX<-W%*%X
WXtAX<-t(WX)%*%AX
dHdrho<- -WXtAX- t(WXtAX)
d2Hdrho2<- 2*t(WX)%*%WX
H<-t(AX)%*%AX
Hinv<-solve(H)
HinvdHdrho<- Hinv%*%dHdrho


V12<-V12 +  0   #beta,rho
V13<-V13 +  0   #beta, sigma
V22<-V22 +  tr.prod(HinvdHdrho,HinvdHdrho) - tr.prod(Hinv,d2Hdrho2)   #rho, rho
V23<-V23 +  0   #rho, sigma
V33<-V33 +  p/2*(1/sigma.hat)    #sigma,sigma
}#end


acvm<-solve(V)
par.names<-c(paste("beta",1:p),"rho","sigma2")
colnames(acvm)<-rownames(acvm)<-par.names
lnlik.model<- - n/2*log(2*pi) - n/2 - n/2*log(sigma.hat)

if(REML){lnlik.model <- lnlik.model + p/2*log(sigma.hat) - log(det(H))}

sigma<-sqrt(sigma.hat)
beta<-beta.hat
rho2<-rho.hat

#gradient
grad<-log.lik.ARerr.gr(rho2,lambda,Y,X,W,n,REML=REML)

return(list(beta=beta,sigma=sigma,rho2=rho2,lnlik.model=lnlik.model,acvm=acvm,grad=grad,counts=counts))

}#end function  AR


################################################################################
### AR pop  ##################################################
################################################################################

AR.pop.exact<-function(x,W,Ys,X,ARmod=1,grad=FALSE,hess=FALSE,pr=TRUE,WtW=t(W)%*%W,WplusWt=W+t(W)){

#K is order
#exact says whether exact inverse or approximation of inverse
ns<-length(Ys)
hilf<-dim(X)
n<-hilf[1]
p<-hilf[2]
p1<-p+1
I<-diag(n)
indNs<-1:ns
indNr<-(ns+1):n


if(pr){
beta<-x[1:p]
rho<-x[p1]
}else{
omega<-abs(x[1])
beta<-x[2:p1]
rho<-x[p+2]
}
Is<-diag(ns)
Ir<-diag(n-ns)

A<-I-rho*W
M<-I-rho*WplusWt+rho^2*WtW
Mss<-M[indNs,indNs]
Msr<-M[indNs,indNr]
Mrs<-M[indNr,indNs]
Mrr<-M[indNr,indNr]
Ass <- A[indNs,indNs]
Asr <- A[indNs,indNr]
Ars <- A[indNr,indNs]
Arr <- A[indNr,indNr]


Mrr.inv<-solve(Mrr)

if(ARmod==1){
# (Ainv)ss= (Ass - Asr*Arrinv*Ars)^{-1}
#therefore   (Ainv)ss^{-1}=Ass - Asr*Arrinv*Ars

Arr.inv<-solve(Arr)
Ass.inv<-solve(Ass)
Ainv.ss<-solve(Ass-Asr%*%Arr.inv%*%Ars)

#similar for (Ainv)rr
Ainv.rr<-solve(Arr-Ars%*%Ass.inv%*%Asr)
# (Ainv)sr=  - Ass.inv * Asr * (Ainv)rr 
Ainv.sr<- - Ass.inv%*%Asr%*%Ainv.rr

Bs<- cbind(Ainv.ss,Ainv.sr)
BsX<-Bs%*%X
}else{
BsX<-X[indNs,]
}
#same formula applies for V= (t(A)*A)^{-1}, M=t(A)*A
Vss<- Mss - Msr %*% Mrr.inv %*% Mrs



mu<-BsX%*%beta
ry<-Ys-mu
Vssry<-Vss%*%ry
rytVssry<-t(ry)%*%Vssry

#   terms                    T1               T2           T3               T4
if(pr){
omega<-c(rytVssry/ns)
}
#it is minus log(detV) to avoid the inverse
detVss<-det(Vss)
L<- 0.5*ns*log(2*pi) + 0.5*ns*log(omega)   - 0.5*log(detVss) + 0.5*rytVssry/omega            
attr(L, "omega") <-   omega


if(grad){

dMdrho<- -WplusWt+2*rho*WtW
dMssdrho<-  dMdrho[indNs,indNs]
dMsrdrho<-  dMdrho[indNs,indNr]
dMrsdrho<-  dMdrho[indNr,indNs]
dMrrdrho<-  dMdrho[indNr,indNr]

Vssinv<-solve(Vss)
BsXtVssry<-t(BsX)%*%Vssry


Mrr.invMrs <- Mrr.inv%*%Mrs 
MsrMrr.inv <- Msr%*%Mrr.inv

if(ARmod==1){
Xbeta<-X%*%beta
Ainv<-solve(A)
AinvW<-Ainv%*%W
AinvWAinv<-AinvW%*%Ainv               
dAinvdrho<-  AinvWAinv  #correct, minus and minus gives plus
 
dBsdrho<-    dAinvdrho[indNs,]
drytdrho <-  - t(dBsdrho%*%Xbeta) #=dBsdrho*drytdBs and  drytdBs=- X*beta


dLdrho <-  drytdrho%*%Vssry/omega
}else{
drytdrho<-matrix(0,1,ns)
dLdrho<-0
}#end if(ARmod==1){
# Vss<- Mss - Msr %*% Mrr.inv %*% Mrs
dVssdrho <- dMssdrho - dMsrdrho%*%Mrr.invMrs + MsrMrr.inv%*%(dMrrdrho%*%Mrr.invMrs - dMrsdrho)
rytdVssdrhory<-t(ry)%*%dVssdrho%*%ry

dLdrho<- dLdrho -0.5*tr.prod(Vssinv,dVssdrho) + 0.5*rytdVssdrhory/omega 

dLdbeta<- - BsXtVssry/omega


par.names<-c("omega",paste("beta",1:p),"rho")

if(pr){
gr.pr<-c(dLdbeta,dLdrho)
names(gr.pr)<-par.names[-1]
attr(L, "gradient") <-   gr.pr
}else{ 
dLdomega<- 0.5*n/omega - 0.5/(omega^2)*rytVssry
gr<-c(dLdomega,dLdbeta,dLdrho)
names(gr)<-par.names
attr(L, "gradient") <-   gr
}

}#end grad



if(grad & hess){
VssinvdVssdrho<- Vssinv%*%dVssdrho
Mrr.invdMrrdrho<- Mrr.inv%*%dMrrdrho
Mrr.invdMrrdrhoMrr.inv<- Mrr.invdMrrdrho%*%Mrr.inv
MsrMrr.inv<-Msr%*%Mrr.inv
Mrr.invMrs<-Mrr.inv%*%Mrs
dMrrdrhoMrr.invdMrrdrho<- dMrrdrho%*%Mrr.invdMrrdrho

d2Mdrho2<- 2*WtW
d2Mssdrho2<-  d2Mdrho2[indNs,indNs]
d2Msrdrho2<-  d2Mdrho2[indNs,indNr]
d2Mrsdrho2<-  d2Mdrho2[indNr,indNs]
d2Mrrdrho2<-  d2Mdrho2[indNr,indNr]

d2Vssdrho2<- d2Mssdrho2 - d2Msrdrho2%*%Mrr.invMrs + 2*dMsrdrho%*%Mrr.invdMrrdrhoMrr.inv%*%Mrs
d2Vssdrho2<-d2Vssdrho2 + 2*Msr%*%Mrr.invdMrrdrhoMrr.inv%*%dMrsdrho - dMsrdrho%*%Mrr.inv%*%dMrsdrho - MsrMrr.inv%*%d2Mrsdrho2 
d2Vssdrho2<-d2Vssdrho2 - 2*MsrMrr.inv%*%Mrr.invdMrrdrhoMrr.inv%*%Mrr.invMrs + MsrMrr.inv%*%d2Mrrdrho2%*%Mrr.invMrs
 
#   "omega"  "beta "   "rho"  
V11 <- 0.5*ns/omega^2
V22 <- t(BsX)%*%Vss%*%BsX/omega
               #first term -1/2*tr(Vssinv*dVssdrho*Vssinv*dVssdrho) - 1/2*tr(Vssin*d2Vssdrho2) +1/2*tr(d2Vssdrho2*Vssinv)/omega + drytdrho*Vss*drydrho/omega
V33 <-  0.5*tr.prod(VssinvdVssdrho,VssinvdVssdrho)         #           - 0.5*(1-1/omega)*tr.prod(Vssinv,d2Vssdrho2) +             drytdrho%*%Vss%*%t(drytdrho)/omega

#omega and beta
V12 <- matrix(0,1,p)
#omega and rho
V13<- - 0.5*tr(VssinvdVssdrho)/omega
#beta and rho
V23<-matrix(0,p,1)

if(ARmod==1){
#additional term because for AR model mean depends on rho
V33 <- V33  + 0.5*drytdrho%*%Vss%*%t(drytdrho)/omega 
V23 <- t(drytdrho%*%Vss%*%BsX)/omega
}



Hess<-rbind(cbind(V11,V12,V13),cbind(t(V12),V22,V23),cbind(t(V13),t(V23),V33))


                  
rownames(Hess)<-colnames(Hess)<- par.names


attr(L, "hessian")  <-    Hess
}#end if(hess){

#return L with possible attributes "grad" and "hess"
L
}#end function  AR.pop.exact




#################################################################################

get.list<-function(L,coeff){
Lnew<-list()
for(k in 1:length(coeff)){
Lnew[[k]]<-coeff[k]*L[[k]]
}
Lnew
}#end get.list


################################################################################
### AR pop  ##################################################
################################################################################

AR.pop.approx<-function(x,WX.list.sample=NULL,Ys,Xs,ARmod=1,grad=FALSE,Info=FALSE,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=length(WplusWt.list.sample)){
  
  if(exists("counter")){counter<<-counter+1;cat("iteration",counter,"\n")}
  
#K is order  of Taylor
# W, WtW and  WplusWt are lists of length K, where kth element refers to powers of W

#exact says whether exact inverse or approximation of inverse
ns<-length(Ys)
hilf<-dim(Xs)

#print(dim(Xs))

#ns<-hilf[1]
p<-hilf[2]
p1<-p+1
Is<-Diagonal(ns,x=rep(1,ns))
#indNs<-1:ns



if(pr){
#beta<-x[1:p]
rho<-x
}else{
omega<-abs(x[1])
beta<-x[2:p1]
rho<-x[p+2]
}


#A<-I-rho*W.list[[1]]


seqK <- 1:K
seqK2<- 2:(2*K)


rho.seq<- rho^seqK 
rho.seq2<- rho^seqK2  


Vssinv <- Is + Reduce("+",get.list(WplusWt.list.sample,rho.seq)) + Reduce("+",get.list(WWt.list.sample,rho.seq2))
Vss<-solve(Vssinv)


#M<-(t(A)%*%A)
#Vss.inv1<-solve(M)[indNs,indNs]

if(ARmod==1){
seqK3<- 1:length(WX.list.sample)
rho.seq3<- rho^seqK3 
#Bs<- cbind(Is,matrix(0,ns,nr))+ Reduce("+",get.list(W.list.sample,rho.seq))[indNs,]
BsX<-Xs+Reduce("+",get.list(WX.list.sample,rho.seq3))
}else{
#cat("ARmod==2")
#  print(Xs)  
BsX<-Xs
#print(BsX)
}


#   terms                    T1               T2           T3               T4
if(pr){
  #print(BsX)
  #print(dim(Vss))
  BsXtVss<-t(BsX)%*%Vss
  beta<-  solve(BsXtVss%*%BsX)%*%BsXtVss%*%Ys

}


mu<-BsX%*%beta
ry<-Ys-mu
Vssry<-Vss%*%ry
rytVssry<-as.numeric(t(ry)%*%Vssry)

if(pr){
  omega<-c(rytVssry/ns)
}

#print(omega)
#print(log(omega))

#it is minus log(detV) to avoid the inverse
detVss<-det(Vss)
L<- c(0.5*ns*log(2*pi) + 0.5*ns*log(omega)   - 0.5*log(detVss) + 0.5*rytVssry/omega)            
attr(L, "omega") <-   omega
attr(L, "beta") <- as.matrix(beta)

if(grad){



BsXtVssry<-t(BsX)%*%Vssry

if(ARmod==1){
d1rho.seq3  <- seqK3*rho^(seqK3-1) #derivative for A^{-1}
#Xbeta<-X%*%beta
#Bs<- cbind(Is,matrix(0,ns,nr))+ Reduce("+",get.list(W.list,rho.seq))[indNs,] 


#dBsdrho<-  Reduce("+",get.list(W.list,d1rho.seq))[indNs,] 
#drytdrho <-  -t(dBsdrho%*%Xbeta)  #=dBsdrho*drytdBs and  drytdBs=- X*beta
drytdrho <-  - t(Reduce("+",get.list(WX.list.sample,d1rho.seq3))%*%beta)


dLdrho <-  drytdrho%*%Vssry/omega

}else{
drytdrho<-matrix(0,1,ns)
dLdrho<-0
}#end if(ARmod==1){
d1rho.seq <- seqK*rho^(seqK-1)  #needed for derivative of Vssinv for WplusWt
d1rho.seq2 <- seqK2*rho^(seqK2-1)  #needed  for derivative of Vssinv for WtW

#Vssinv <- Is + Reduce("+",get.list(WplusWt.list,rho.seq))[indNs,indNs] + Reduce("+",get.list(WWt.list,rho.seq2))[indNs,indNs]
dVssinvdrho<-Reduce("+",get.list(WplusWt.list.sample,d1rho.seq)) + Reduce("+",get.list(WWt.list.sample,d1rho.seq2))
dVssdrho <- - Vss%*%dVssinvdrho%*%Vss
rytdVssdrhory<-t(ry)%*%dVssdrho%*%ry

dLdrho<- as.numeric(dLdrho -0.5*tr.prod(Vssinv,dVssdrho) + 0.5*rytdVssdrhory/omega )

dLdbeta<- - c(as.matrix(BsXtVssry)/omega)


par.names<-c("omega",paste("beta",1:p),"rho")

#if(pr){
#gr.pr<-c(dLdrho)
#names(gr.pr)<-"rho"
#attr(L, "gradient") <-   gr.pr
#}else{ 
dLdomega<- 0.5*ns/omega - 0.5*rytVssry/(omega^2)
#print(dLdomega)
#print(dLdrho)
gr<-c(dLdomega,c(dLdbeta),dLdrho)
names(gr)<-par.names
attr(L, "gradient") <-   gr
#}

}#end grad



if(grad & Info){
#d2rho.seq  <- seqK*(seqK-1)*rho^(seqK-2)
#d2rho.seq2 <- seqK2*(seqK2-1)*rho^(seqK2-2) 

VssinvdVssdrho<- Vssinv%*%dVssdrho


#d2Vssinvdrho2<-  Reduce("+",get.list(WplusWt.list.sample,d1rho.seq)) + Reduce("+",get.list(WWt.list.sample,d1rho.seq2))


#d2Vssdrho2<-  Vss%*%(2*dVssinvdrho%*%Vss%*%dVssinvdrho - d2Vssinvdrho2)%*%Vss

 
#   "omega"  "beta "   "rho"  
V11 <- 0.5*ns/omega^2  #omega,omega
V22 <- t(BsX)%*%Vss%*%BsX/omega  #beta,beta

              #first term -1/2*tr(Vssinv*dVssdrho*Vssinv*dVssdrho) - 1/2*tr(Vssin*d2Vssdrho2) +1/2*tr(d2Vssdrho2*Vssinv)/omega + drytdrho*Vss*drydrho/omega
V33 <-  0.5*tr.prod(VssinvdVssdrho,VssinvdVssdrho)     #               - 0.5*(1-1/omega)*tr.prod(Vssinv,d2Vssdrho2) +             drytdrho%*%Vss%*%t(drytdrho)/omega
#rho,rho


#omega and beta
V12<-Matrix(0,1,p)
#omega and rho
V13<- -dLdrho/omega
#beta and rho
V23<-Matrix(0,p,1)

if(ARmod==1){
#additional term because for AR model mean depends on rho
V33 <- V33  + drytdrho%*%Vss%*%t(drytdrho)/omega  #rho,rho
V23 <- - t(drytdrho%*%Vss%*%BsX)/omega    #beta,rho
}


#print(V12)
#print(V22)
#print(V23)

#print(dim(V12))
#print(dim(V22))
#print(dim(V23))

#print(cBind(V11,V12,V13))
#print(cBind(t(V12),V22,V23))
#print(cBind(t(V13),t(V23),V33))

Hess<-as.matrix(rBind(cBind(V11,V12,V13),cBind(t(V12),V22,V23),cBind(t(V13),t(V23),V33)))


                  
rownames(Hess)<-colnames(Hess)<- par.names


attr(L, "Info")  <-    Hess
}#end if(hess){

#return L with possible attributes "grad" and "hess"
L
}#end function  AR.pop




################################################################################
### AR model with random effects ##################################################
################################################################################


AR.random<-function(x,Y,X,W,Z1,Z2=NULL,lambda,grad=TRUE,hess=FALSE,simple=TRUE,pr=TRUE){
#x: first omega, then beta, then rho,  then random effects theta1 and theta2
#Z1 is list of block matrices for random effects (HH)
#Z2 is list of block matrices for random effetcs (area)
#Model: Y=theta*W*Y+X*beta+Z*u+epsilon
# epsilon~N(0,omega), u~N(0,omega*D)
if(exists("counter")){counter<<-counter+1;cat("Iter:",counter,"\n")}
  
hilf<-dim(X)
p<-hilf[2]
p1<-p+1
n<-hilf[1]
#print(p)
#print(n)
isnullZ2<-is.null(Z2)
if(isnullZ2){q<-1}else{q<-2}
if(pr){
rho<-x[1] #rho is first value
theta<-x[2:(1+q)]  #2nd and 3rd are theta
#print(beta);print(theta);print(rho)
}else{
omega<-abs(x[1])
beta<-x[2:p1]
rho<-x[p+2]
theta<-abs(x[(p+3):length(x)]) #theta is of length 2, associated with Z1 and Z2
}
#theta<-pmax(theta,0)
Y<-as.matrix(Y)

#- loglik=n/2*log(2*pi)+2/n*log(omega)- log |A| + 1/2*log|V| + rzt*Vinv*rz/omega
diagn<-diag(n) 
A<-diagn - rho*W 
AY<-A%*%Y


hilf<- get.Vadd(Z1,Z2,theta)
V<-hilf$V
Vinv<-hilf$Vinv

logdetV<-l.det.V(V)

if(pr){
   #obtain beta
  XtVinv<-prod.vector.list.vector(X=X,Z1=Vinv)
  XtVinvAY<-XtVinv%*%AY
  XtVinvX<-XtVinv%*%X
  beta<-as.matrix(solve(XtVinvX,XtVinvAY))
}
Xbeta<-X%*%beta

rz<- as.matrix(AY - Xbeta)
log.detA<-Re(sum(log(1-rho*lambda)))
rztVinv<-prod.vector.list.vector(X=rz,Z1=Vinv)
rztVinvrz<- sum(t(rztVinv)*rz)

#   terms                    T1               T2           T3               T4
if(pr){
omega<-rztVinvrz/n
}
L<- 0.5*n*log(2*pi) + 0.5*n*log(omega) - log.detA  + 0.5*logdetV + 0.5*rztVinvrz/omega            
attr(L, "omega") <-   omega
attr(L,"beta")<-beta

#now compute gradient
##############################################
if(grad){
WY<-W%*%Y
if(!pr){
XtVinv<-prod.vector.list.vector(X=X,Z1=Vinv)
}
XtVinvrz<-XtVinv%*%rz
dlog.detAdrho<- - Re(sum(lambda/(1-rho*lambda)))
WYtVinv<-prod.vector.list.vector(X=WY,Z1=Vinv)
WYtVinvrz<-sum(t(WYtVinv)*rz)
C1<- prod.list.list.list(Z1=Vinv,Z2=Z1,Z3=Vinv)
rztC1rz<-prod.vector.list.vector(X=rz,Z1=C1,Z2=NULL,Y=rz)
tr.VinvZ1<-tr.prod.list(V1=Vinv,V2=Z1)

#standard gradient

dLdomega <- (n/2)*(1/omega) - 0.5*rztVinvrz/(omega^2)
dLdbeta  <- - XtVinvrz/omega
dLdrho   <-  -dlog.detAdrho - WYtVinvrz/omega
dLdtheta <- as.numeric(0.5*tr.VinvZ1 -1/(2*omega)*rztC1rz)
if(!isnullZ2){
C2<- prod.list.list.list(Z1=Vinv,Z2=Z2,Z3=Vinv)
rztC2rz<-prod.vector.list.vector(X=rz,Z1=C2,Z2=NULL,Y=rz)
tr.VinvZ2<-tr.prod.list(V1=Vinv,V2=Z2)
dLdtheta <- c(dLdtheta,as.numeric(0.5*tr.VinvZ2 -1/(2*omega)*rztC2rz))
}#end if(!isnullZ2){
par.names<-c("omega",paste("beta",1:p),"rho",paste("theta",1:q))

gr<-c(as.numeric(dLdomega),as.numeric(dLdbeta),as.numeric(dLdrho),as.numeric(dLdtheta))

names(gr)<-par.names
attr(L, "gradient") <-   gr

if(pr){
gr.pr<-c(as.numeric(dLdrho),as.numeric(dLdtheta))
names(gr.pr)<-c("rho",paste("theta",1:q))
attr(L, "gradient") <-   gr.pr
}



}#end if grad


#now compute hessian
##############################################
if(hess & grad){
#Ainv<-solve(A)
XtVinvX<- XtVinv%*%X
AinvXbeta<- solve(A,Xbeta)
XtVinvWAinvXbeta<-  (XtVinv%*%W)%*%AinvXbeta
d2log.detAdrho2<- - Re(sum(lambda^2/(1-rho*lambda)^2))



d2Ldomega2<-(n/2)*(1/omega^2)
d2Ldomegadbeta<- rep(0,p)
d2Ldomegadrho<- (1/omega) * tr.prod(W,Ainv)
d2Ldomegadtheta<-(1/(2*omega))*tr.VinvZ1

d2Ldbeta2      <- XtVinvX/omega
d2Ldbetadrho   <- (1/omega)*XtVinvWAinvXbeta
d2Ldbetadtheta <- matrix(0,p,q)

#here just take observed information, easier to compute
if(simple){
WYtVinvWY<-sum(t(WYtVinv)*WY)
WYtC1rz<- prod.vector.list.vector(X=WY,Z1=C1,Z2=NULL,Y=rz)
d2Ldrho2 <-  - d2log.detAdrho2 + (1/omega)*WYtVinvWY
d2Ldrhodtheta <- as.numeric((1/omega)*WYtC1rz)
}else{
Ainv<-solve(A)
WtVinvW<-prod.vector.list.vector(X=W,Z1=Vinv,Z2=NULL,Y=W)
AinvVAinvt<- prod.vector.list.vector(X=t(Ainv),Z1=V,Z2=NULL,Y=t(Ainv))    

AinvZ1Vinv<-prod.vector.list.vector(X=t(Ainv),Z1=Z1,Z2=Vinv,Y=NULL)

d2Ldrho2 <-  - d2log.detAdrho2 + (1/omega)*tr.prod(WtVinvW,omega*AinvVAinvt + AinvXbeta%*%t(AinvXbeta))
d2Ldrhodtheta <- (1/omega)*tr.prod(AinvZ1Vinv,W)
}

#here again expected observation matrix
d2Ldtheta2 <- 0.5*tr.prod.list(V1=C1,V2=Z1)


if(!isnullZ2){

d2Ldomegadtheta<-c(d2Ldomegadtheta,1/(2*omega)*tr.VinvZ2)

if(simple){
WYtC2rz<- prod.vector.list.vector(X=WY,Z1=C2,Z2=NULL,Y=rz)
d2Ldrhodtheta <- c(d2Ldrhodtheta,(1/omega)*as.numeric(WYtC2rz))
}else{
AinvZ2Vinv<-prod.vector.list.vector(X=t(Ainv),Z1=Z2,Z2=Vinv,Y=NULL)
d2Ldrhodtheta <- c(d2Ldrhodtheta,(1/omega)*tr.prod(AinvZ2Vinv,W))
}

h12<-0.5*tr.prod.list(V1=C1,V2=Z2)
h22<-0.5*tr.prod.list(V1=C2,V2=Z2)
d2Ldtheta2 <- matrix(c(d2Ldtheta2,h12,h12,h22),2,2,byrow=TRUE)
                                     
}#end if(!isnullZ2){

Hess1<-  c(d2Ldomega2     , d2Ldomegadbeta,d2Ldomegadrho,d2Ldomegadtheta)
Hess2<-  as.matrix(cBind(d2Ldomegadbeta , d2Ldbeta2     ,d2Ldbetadrho ,d2Ldbetadtheta))
Hess3<-  c(d2Ldomegadrho  , as.numeric(d2Ldbetadrho)  , d2Ldrho2     ,as.numeric(d2Ldrhodtheta))
#print(dim())
if(!isnullZ2){
Hess4<-  cBind(d2Ldomegadtheta, t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}else{
Hess4<-  cBind(d2Ldomegadtheta, t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}

Hess<-  rBind(Hess1,Hess2,Hess3,Hess4)
Hess<-as.matrix(Hess)                
                  



rownames(Hess)<-colnames(Hess)<- par.names


attr(L, "hessian")  <-    Hess
}#end if(hess){

L
}#end  AR.random
################################################################################
# end AR random
################################################################################


################################################################################
### ARerr1 model with random effects ##################################################
################################################################################


ARerr1.random<-function(x,Y,X,W,Z1,Z2=NULL,lambda,grad=TRUE,hess=FALSE,simple=TRUE,pr=TRUE,REML=TRUE){
#x: first omega, then beta, then rho,  then random effects theta1 and theta2
#Z1 is list of block matrices for random effects (HH)
#Z2 is list of block matrices for random effetcs (area)
#Model: Y=X*beta+epsilon
# epsilon=theta*W*epsilon + Z*u + v, v~N(0,omega), u~N(0,omega*D)
if(exists("counter")){counter<<-counter+1;cat("Iter:",counter,"\n")}
hilf<-dim(X)
p<-hilf[2]
p1<-p+1
n<-hilf[1]
#print(p)
#print(n)
isnullZ2<-is.null(Z2)
if(isnullZ2){q<-1}else{q<-2}
if(pr){
rho<-x[1]
theta<-abs(x[2:length(x)])
}else{
omega<-abs(x[1])
beta<-x[2:p1]
rho<-x[p+2]
theta<-abs(x[(p+3):length(x)]) #theta is of length 2, associated with Z1 and Z2
}
#theta<-pmax(theta,0)
Y<-as.matrix(Y)

#- loglik=n/2*log(2*pi)+2/n*log(omega)- log |A| + 1/2*log|V| + rzt*Vinv*rz/omega
diagn<-Diagonal(n)
A<-diagn - rho*W 

AX<-A%*%X
AY<-A%*%Y

hilf<- get.Vadd(Z1,Z2,theta)
V<-hilf$V
Vinv<-hilf$Vinv

logdetV<-l.det.V(V)
if(pr){
  #obtain beta
  XtAtVinv<-prod.vector.list.vector(X=AX,Z1=Vinv)
  XtAtVinvAX<-XtAtVinv%*%AX
  XtAtVinvAY<-XtAtVinv%*%AY
  beta<-as.matrix(solve(XtAtVinvAX,XtAtVinvAY))
}
Xbeta<-X%*%beta

ry<- as.matrix(Y - Xbeta)
rz<-A%*%ry
log.detA<-Re(sum(log(1-rho*lambda)))
rztVinv<-prod.vector.list.vector(X=rz,Z1=Vinv)
rztVinvrz<- sum(t(rztVinv)*rz)

#   terms                    T1               T2           T3               T4
if(pr){
omega<-rztVinvrz/n
if(REML){
omega<-rztVinvrz/(n-p)
}
}
L<- 0.5*n*log(2*pi) + 0.5*n*log(omega) - log.detA  + 0.5*logdetV + 0.5*rztVinvrz/omega            

if(REML){
#Sigma=Ainv%*%V%*%t(Ainv)
#Sigmainv=A%*%Vinv%*%t(A)
#H=t(X)%*%Sigmainv%*%X
#AX<-A%*%X
H<- prod.vector.list.vector(X=AX,Z1=Vinv,Z2=NULL,Y=AX)
L<-L - 0.5*p*log(omega)+0.5*log(det(H))
}


attr(L, "omega") <-   omega
attr(L, "beta") <-   beta
#now compute gradient
##############################################
if(grad){
  
AX<-A%*%X   
Wry<-W%*%ry
AXtVinv<-prod.vector.list.vector(X=AX,Z1=Vinv)
AXtVinvrz<-AXtVinv%*%rz   
dlog.detAdrho<- - Re(sum(lambda/(1-rho*lambda)))
WrytVinv<-prod.vector.list.vector(X=Wry,Z1=Vinv)
WrytVinvrz<-sum(t(WrytVinv)*rz)


C1<- prod.list.list.list(Z1=Vinv,Z2=Z1,Z3=Vinv)
rztC1rz<-prod.vector.list.vector(X=rz,Z1=C1,Z2=NULL,Y=rz)
tr.VinvZ1<-tr.prod.list(V1=Vinv,V2=Z1)

#standard gradient

dLdomega <- (n/2)*(1/omega) - 0.5*rztVinvrz/(omega^2)
dLdbeta  <- as.matrix(- AXtVinvrz/omega)
dLdrho   <-  -dlog.detAdrho - WrytVinvrz/omega
dLdtheta <-  0.5*tr.VinvZ1 -1/(2*omega)*as.numeric(rztC1rz)
if(!isnullZ2){
C2<- prod.list.list.list(Z1=Vinv,Z2=Z2,Z3=Vinv)
rztC2rz<-prod.vector.list.vector(X=rz,Z1=C2,Z2=NULL,Y=rz)
tr.VinvZ2<-tr.prod.list(V1=Vinv,V2=Z2)
dLdtheta <- c(dLdtheta,0.5*tr.VinvZ2 -1/(2*omega)*as.numeric(rztC2rz))
}#end if(!isnullZ2){


if(REML){
Hinv<-solve(H)
WX<-W%*%X
dHdrho<-  prod.vector.list.vector(X=-AX,Z1=Vinv,Z2=NULL,Y=WX)        
dHdrho<-  dHdrho + t(dHdrho)



dLdomega <-  dLdomega  -  (p/2)*(1/omega)   
dLdrho   <-  dLdrho    +  0.5*tr.prod(Hinv,dHdrho)

if(!isnullZ2){
dHdtheta1<- prod.vector.list.vector(X=-AX,Z1=C1,Z2=NULL,Y=AX)       
dHdtheta2<- prod.vector.list.vector(X=-AX,Z1=C2,Z2=NULL,Y=AX)       
dLdtheta <- dLdtheta +  0.5*c(tr.prod(Hinv,dHdtheta1),tr.prod(Hinv,dHdtheta2))
}else{
dHdtheta1<-  prod.vector.list.vector(X=-AX,Z1=C1,Z2=NULL,Y=AX)       
dLdtheta <-  dLdtheta  +  0.5*tr.prod(Hinv,dHdtheta1)
}
}#end REML


par.names<-c("omega",paste("beta",1:p),"rho",paste("theta",1:q))

gr<-c(dLdomega,dLdbeta,dLdrho,dLdtheta)

names(gr)<-par.names
attr(L, "gradient") <-   gr

if(pr){
gr.pr<-c(dLdrho,dLdtheta)
names(gr.pr)<-c("rho",paste("theta",1:q))
attr(L, "gradient") <-   gr.pr
}
}#end if grad


#now compute hessian
##############################################
if(hess & grad){
AXtVinvAX<- AXtVinv%*%AX
AinvXbeta<- solve(A,Xbeta)
#XtVinvWAinvXbeta<-  (XtVinv%*%W)%*%AinvXbeta
d2log.detAdrho2<- - Re(sum(lambda^2/(1-rho*lambda)^2))



d2Ldomega2<-(n/2)*(1/omega^2)
d2Ldomegadbeta<- rep(0,p)
d2Ldomegadrho<- (1/omega) * tr.prod(W,Ainv)
d2Ldomegadtheta<-(1/(2*omega))*tr.VinvZ1

d2Ldbeta2      <- AXtVinvAX/omega
d2Ldbetadrho   <- rep(0,p)
d2Ldbetadtheta <- matrix(0,p,q)

#here just take observed information, easier to compute
if(simple){
WrytVinvWry<-sum(t(WrytVinv)*Wry)
WrytC1rz<- prod.vector.list.vector(X=Wry,Z1=C1,Z2=NULL,Y=rz)
d2Ldrho2 <-  - d2log.detAdrho2 + (1/omega)*WrytVinvWry
d2Ldrhodtheta <- (1/omega)*as.numeric(WrytC1rz)
}else{
Ainv<-solve(A)
WtVinvW<-prod.vector.list.vector(X=W,Z1=Vinv,Z2=NULL,Y=W)
AinvVAinvt<- prod.vector.list.vector(X=t(Ainv),Z1=V,Z2=NULL,Y=t(Ainv))    

AinvZ1Vinv<-prod.vector.list.vector(X=t(Ainv),Z1=Z1,Z2=Vinv,Y=NULL)

d2Ldrho2 <-  - d2log.detAdrho2 + tr.prod(WtVinvW,AinvVAinvt)
d2Ldrhodtheta <- (1/omega)*tr.prod(AinvZ1Vinv,W)
}

#here again expected observation matrix
d2Ldtheta2 <- 0.5*tr.prod.list(V1=C1,V2=Z1)


if(!isnullZ2){

d2Ldomegadtheta<-c(d2Ldomegadtheta,1/(2*omega)*tr.VinvZ2)

if(simple){
WrytC2rz<- prod.vector.list.vector(X=Wry,Z1=C2,Z2=NULL,Y=rz)
d2Ldrhodtheta <- c(d2Ldrhodtheta,(1/omega)*as.matrix(WrytC2rz))
}else{
AinvZ2Vinv<-prod.vector.list.vector(X=t(Ainv),Z1=Z2,Z2=Vinv,Y=NULL)
d2Ldrhodtheta <- c(d2Ldrhodtheta,(1/omega)*tr.prod(AinvZ2Vinv,W))
}

h12<-0.5*tr.prod.list(V1=C1,V2=Z2)
h22<-0.5*tr.prod.list(V1=C2,V2=Z2)
d2Ldtheta2 <- matrix(c(d2Ldtheta2,h12,h12,h22),2,2,byrow=TRUE)
                                     
}#end if(!isnullZ2){


if(REML){
d2Ldomega2<-d2Ldomega2 + (p/2)/omega^2


if(!isnullZ2){


#theta1,theta2
HinvdHdtheta1<-Hinv%*%dHdtheta1
HinvdHdtheta2<-Hinv%*%dHdtheta2


C1VC1<- prod.list.list.list(Z1=C1,Z2=V,Z3=C1)
C1VC2<- prod.list.list.list(Z1=C1,Z2=V,Z3=C2)
C2VC2<- prod.list.list.list(Z1=C2,Z2=V,Z3=C2)


d2Hdtheta12<- 2*prod.vector.list.vector(X=AX,Z1=C1VC1,Z2=NULL,Y=AX) 
d2Hdtheta1dtheta2<- 2*prod.vector.list.vector(X=AX,Z1=C1VC2,Z2=NULL,Y=AX) 
d2Hdtheta22<- 2*prod.vector.list.vector(X=AX,Z1=C2VC2,Z2=NULL,Y=AX) 


#print(dim(HinvdHdtheta1))
#print(dim(Hinv))
#print(dim(HinvdHdtheta2))
#print(dim(d2Hdtheta1dtheta2))
#print(dim(d2Hdtheta12))
#print(dim(d2Hdtheta22))

h11<-tr.prod(-HinvdHdtheta1,HinvdHdtheta1) + tr.prod(Hinv,d2Hdtheta12)
h12<-tr.prod(-HinvdHdtheta1,HinvdHdtheta2) + tr.prod(Hinv,d2Hdtheta1dtheta2)
h22<-tr.prod(-HinvdHdtheta2,HinvdHdtheta2) + tr.prod(Hinv,d2Hdtheta22)
d2Ldtheta2<-d2Ldtheta2 + 0.5*matrix(c(h11,h12,h12,h22),2,2)

#print(matrix(c(h11,h12,h12,h22),2,2))


#rho,theta1
HinvdHdrho<-Hinv%*%dHdrho
WXtC1AX<-prod.vector.list.vector(X=WX,Z1=C1,Z2=NULL,Y=AX)
WXtC2AX<-prod.vector.list.vector(X=WX,Z1=C2,Z2=NULL,Y=AX)  
d2Hdrhodtheta12<- (WXtC1AX+t(WXtC1AX))
d2Hdrhodtheta22<- (WXtC2AX+t(WXtC2AX)) 


 
h1<-0.5*tr.prod(-HinvdHdrho,HinvdHdtheta1) + 0.5*tr.prod(Hinv,d2Hdrhodtheta12) 

h2<-0.5*tr.prod(-HinvdHdrho,HinvdHdtheta2) + 0.5*tr.prod(Hinv,d2Hdrhodtheta22)

d2Ldrhodtheta<-d2Ldrhodtheta +  c(h1,h2)

#rho2
Sigmainv<-prod.vector.list.vector(X=A,Z1=Vinv,Z2=NULL,Y=A)
AinvW<-Ainv%*%W
AinvWAinv<-AinvW%*%Ainv
AinvV<-prod.vector.list.vector(X=t(Ainv),Z1=V,Z2=NULL,Y=NULL) # Ainv%*%V
AinvWAinvWAinv<-AinvW%*%AinvWAinv
AinvWAinvWAinvVAinvt<-AinvWAinvWAinv%*%t(AinvV)
AinvWAinvVAinvtWtAinvt<-prod.vector.list.vector(X=t(AinvWAinv),Z1=V,Z2=NULL,Y=t(AinvWAinv))
AinvWAinvVAinvt<-AinvWAinv%*%t(AinvV)

dSigmadrho<-  AinvWAinvVAinvt + t(AinvWAinvVAinvt)
d2Sigmadrho2<-2*(AinvWAinvWAinvVAinvt+AinvWAinvVAinvtWtAinvt+t(AinvWAinvWAinvVAinvt))
d2Hdrho2<-t(X)%*%(Sigmainv%*%(2*dSigmadrho%*%Sigmainv%*%dSigmadrho- d2Sigmadrho2)%*%Sigmainv)%*%X  
d2Ldrho2<-d2Ldrho2 + 0.5*tr.prod(-HinvdHdrho,HinvdHdrho) + 0.5*tr.prod(Hinv,d2Hdrho2) 



}else{ 
#theta1,theta1
C1VC1<- prod.list.list.list(Z1=C1,Z2=V,Z3=C1)
d2Hdtheta12<- 2*prod.vector.list.vector(X=AX,Z1=C1VC1,Z2=NULL,Y=AX)


HinvdHdtheta1<-Hinv%*%dHdtheta1
d2Ldtheta2<-d2Ldtheta2 + 0.5*tr.prod(-HinvdHdtheta1,HinvdHdtheta1)+0.5*tr.prod(Hinv,d2Hdtheta12)



#rho,theta1
HinvdHdrho<-Hinv%*%dHdrho
WXtC1AX<-prod.vector.list.vector(X=WX,Z1=C1,Z2=NULL,Y=AX) 
d2Hdrhodtheta12<- (WXtC1AX+t(WXtC1AX)) 
d2Ldrhodtheta<-d2Ldrhodtheta + 0.5*tr.prod(-HinvdHdrho,HinvdHdtheta1) + 0.5*tr.prod(Hinv,d2Hdrhodtheta12) 

#rho2
Sigmainv<-prod.vector.list.vector(X=A,Z1=Vinv,Z2=NULL,Y=A)
AinvW<-Ainv%*%W
AinvWAinv<-AinvW%*%Ainv
AinvV<-prod.vector.list.vector(X=t(Ainv),Z1=V,Z2=NULL,Y=NULL) # Ainv%*%V
AinvWAinvWAinv<-AinvW%*%AinvWAinv
AinvWAinvWAinvVAinvt<-AinvWAinvWAinv%*%t(AinvV)
AinvWAinvVAinvtWtAinvt<-prod.vector.list.vector(X=t(AinvWAinv),Z1=V,Z2=NULL,Y=t(AinvWAinv))
AinvWAinvVAinvt<-AinvWAinv%*%t(AinvV)

dSigmadrho<-  AinvWAinvVAinvt + t(AinvWAinvVAinvt)
d2Sigmadrho2<-2*(AinvWAinvWAinvVAinvt+AinvWAinvVAinvtWtAinvt+t(AinvWAinvWAinvVAinvt))
d2Hdrho2<-t(X)%*%(Sigmainv%*%(2*dSigmadrho%*%Sigmainv%*%dSigmadrho- d2Sigmadrho2)%*%Sigmainv)%*%X  
d2Ldrho2<-d2Ldrho2 + 0.5*tr.prod(-HinvdHdrho,HinvdHdrho) + 0.5*tr.prod(Hinv,d2Hdrho2) 

}#end  if(!isnullZ2){

}#end REML




Hess1<-  c(d2Ldomega2     , d2Ldomegadbeta,d2Ldomegadrho,d2Ldomegadtheta)
Hess2<-  cbind(d2Ldomegadbeta , as.matrix(d2Ldbeta2)     ,d2Ldbetadrho ,d2Ldbetadtheta)
Hess3<-  c(d2Ldomegadrho  , d2Ldbetadrho  ,d2Ldrho2     ,d2Ldrhodtheta)
if(!isnullZ2){
Hess4<-  cbind(d2Ldomegadtheta, t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}else{
Hess4<-  c(d2Ldomegadtheta, t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}

Hess<-  rbind(Hess1,Hess2,Hess3,Hess4)
                  
                  



rownames(Hess)<-colnames(Hess)<- par.names


attr(L, "hessian")  <-    Hess
}#end if(hess){

L
#########################################################################
}#end  ARerr1.random
################################################################################

################################################################################
### ARerr2 model with random effects ##################################################
################################################################################


ARerr2.random<-function(x,Y,X,W,Z1,Z2=NULL,lambda,grad=TRUE,hess=FALSE,simple=TRUE,pr=TRUE,REML=TRUE){
#x: first omega, then beta, then rho,  then random effects theta1 and theta2
#Z1 is list of block matrices for random effects (HH)
#Z2 is list of block matrices for random effetcs (area)
  if(exists("counter")){counter<<-counter+1;cat("Iter:",counter,"\n")}
hilf<-dim(X)
p<-hilf[2]
p1<-p+1
n<-hilf[1]
#print(p)
#print(n)
isnullZ2<-is.null(Z2)
if(isnullZ2){q<-1}else{q<-2}
if(pr){
rho<-x[1]
theta<-abs(x[2:(q+1)])
cat("rho:",rho,"\n")
cat("theta:",theta,"\n")
}else{
omega<-abs(x[1])
beta<-x[2:p1]
rho<-x[p+2]
theta<-abs(x[(p+3):length(x)]) #theta is of length 2, associated with Z1 and Z2
}
#theta<-pmax(theta,0)
Y<-as.matrix(Y)

#cat("rho:",rho,"\n")
#cat("theta:",theta,"\n")

#- loglik=n/2*log(2*pi)+2/n*log(omega)- 1/2*log|V| + ryt*Vinv*ry/omega
diagn<-Diagonal(n) 
A<-diagn - rho*W 
Xbeta<-X%*%beta
#Ainv<-solve(A)
#AinvAinvt<-Ainv%*%t(Ainv)

hilf<- get.Vadd(Z1,Z2,theta)
V<-hilf$V
V<-get.Vblock(V)
diag(V)<-diag(V)-1

Vnew<-forceSymmetric(A%*%V%*%A + Diagonal(n))
#Vinv<-solve(V)
#detV<-det(V)
#cat("rho:",rho,"\n")
#cat("theta:",theta,"\n")
chol.Vnew<-chol(Vnew)
AtA<-forceSymmetric(t(A)%*%A)
chol.AtA<-chol(AtA)

if(pr){
AX<-A%*%X
AY<-A%*%Y
#beta<- solve(t(AX)%*%solve(Vnew)%*%AX)%*%t(AX)%*%solve(Vnew)%*%AY
b<-solve(t(chol.Vnew),AX)
a<-solve(t(chol.Vnew),AY)
beta<- as.matrix(solve(t(b)%*%b)%*%t(b)%*%a) 
}

ry<- as.matrix(Y - Xbeta)
Ary<-A%*%ry


hilf<- solve(t(chol.Vnew),Ary)
rytVinvry<- as.numeric(t(hilf)%*%hilf)

#   terms                    T1               T2              T3              
if(pr){
omega<-rytVinvry/n
if(REML){
omega<-rytVinvry/(n-p)
}
}

#print(log(omega))
#print(detV)
#print(rytVinvry)

L<- 0.5*n*log(2*pi) + 0.5*n*log(omega)  + 0.5*2*sum(log(diag(chol.Vnew)))-0.5*2*sum(log(diag(chol.AtA))) + 0.5*rytVinvry/omega            
if(REML){
H<- forceSymmetric(t(b)%*%b)
L<-L - 0.5*p*log(omega)+0.5*2*sum(log(diag(chol(H))))
}
attr(L, "omega") <-   omega
attr(L,"beta")<-beta
#print(L)

#now compute gradient
##############################################
if(grad){
Ainv<-solve(A)
AinvAinvt<-Ainv%*%t(Ainv)
V0<-V+AinvAinvt
V0inv<-solve(V0)
AinvW<-Ainv%*%W
AinvWAinvAinvt<-AinvW%*%AinvAinvt

XtVinv<-t(X)%*%V0inv
XtVinvry<-XtVinv%*%ry

Z1<-get.Vblock(Z1)
Z2<-get.Vblock(Z2)

VinvZ1<-V0inv%*%Z1
C1<-VinvZ1%*%V0inv  
rytC1ry<-t(ry)%*%C1%*%ry
tr.VinvZ1<-tr.prod(V0inv,Z1)

dVdrho<-AinvWAinvAinvt+t(AinvWAinvAinvt)
rytVinv<-t(ry)%*%V0inv
#standard gradient

dLdomega <- (n/2)*(1/omega) - 0.5*rytVinvry/(omega^2)
dLdbeta  <-  as.matrix(- XtVinvry/omega)
dLdrho   <-  as.numeric( 0.5*tr.prod(V0inv,dVdrho) -1/(2*omega)*(rytVinv)%*%dVdrho%*%t(rytVinv) )
dLdtheta <-  as.numeric(0.5*tr(VinvZ1) -1/(2*omega)*rytC1ry)



if(!isnullZ2){
VinvZ2<-V0inv%*%Z2
C2<-VinvZ2%*%V0inv 
rytC2ry<-t(ry)%*%C2%*%ry
tr.VinvZ2<-tr.prod(V0inv,Z2)
dLdtheta <- c(dLdtheta,0.5*tr.VinvZ2 -1/(2*omega)*as.numeric(rytC2ry))
}#end if(!isnullZ2){

if(REML){

dVdtheta<- Z1


Hinv<-solve(H)


#dVdrho already defined
dHdrho<- -XtVinv%*%dVdrho%*%t(XtVinv)


dLdomega<- dLdomega -  (p/2)*(1/omega)   
dLdrho<-   dLdrho   +  0.5*tr.prod(Hinv,dHdrho)

#Z1 is dVdtheta1
dHdtheta1 <-  -t(X)%*%C1%*%X  #-XtVinv%*%Z1%*%t(XtVinv)


if(!isnullZ2){
dHdtheta2 <- -t(X)%*%C2%*%X  #-XtVinv%*%Z2%*%t(XtVinv)
dLdtheta<- dLdtheta +  0.5*c(tr.prod(Hinv,dHdtheta1),tr.prod(Hinv,dHdtheta2))
}else{
dLdtheta<- dLdtheta +  0.5*tr.prod(Hinv,dHdtheta1)
}
}#end REML


par.names<-c("omega",paste("beta",1:p),"rho",paste("theta",1:q))

gr<-c(dLdomega,dLdbeta,dLdrho,dLdtheta)

names(gr)<-par.names
attr(L, "gradient") <-   gr

if(pr){
gr.pr<-c(dLdrho,dLdtheta)
names(gr.pr)<-par.names[-(1:(p+1))]
attr(L, "gradient") <-   gr.pr
}
}#end if grad


#now compute hessian
##############################################
if(hess & grad){
AinvWAinv<-AinvW%*%Ainv
AinvWAinvWAinvAinvt<-AinvW%*%AinvWAinv%*%t(Ainv)
XtVinvX<- XtVinv%*%X
VinvdVdrho<-V0inv%*%dVdrho
tr.VinvdVdrho<-tr(VinvdVdrho)
d2Vdrho2<-2*(AinvWAinv%*%t(AinvWAinv)+AinvWAinvWAinvAinvt+t(AinvWAinvWAinvAinvt))

d2Ldomega2<-(n/2)*(1/omega^2)
d2Ldomegadbeta<- rep(0,p)
d2Ldomegadrho<- (1/(2*omega)) * tr.VinvdVdrho
d2Ldomegadtheta<-(1/(2*omega))*tr.VinvZ1

d2Ldbeta2      <- as.matrix(XtVinvX/omega)
d2Ldbetadrho   <- rep(0,p)
d2Ldbetadtheta <- matrix(0,p,q)

#here just take observed information, easier to compute

d2Ldrho2 <-  0.5*tr.prod(VinvdVdrho,VinvdVdrho)-(0.5-1/(2*omega))*tr.prod(V0inv,d2Vdrho2)
d2Ldrhodtheta <-0.5*tr.prod(VinvdVdrho,VinvZ1) 


#here again expected observation matrix
d2Ldtheta2 <-0.5*tr.prod(VinvZ1,VinvZ1) 


if(!isnullZ2){

d2Ldomegadtheta<-c(d2Ldomegadtheta,1/(2*omega)*tr.VinvZ2)
d2Ldrhodtheta <- c(d2Ldrhodtheta,(-0.5+1/omega)*tr.prod(VinvdVdrho,VinvZ2) )


h12<-0.5*tr.prod(VinvZ1,VinvZ2) 
h22<-0.5*tr.prod(VinvZ2,VinvZ2) 
d2Ldtheta2 <- matrix(c(d2Ldtheta2,h12,h12,h22),2,2,byrow=TRUE)
                                     
}#end if(!isnullZ2){


if(REML){
#d2Vdrho2  already defined

#rho,rho
HinvdHdrho<-Hinv%*%dHdrho
d2Hdrho2<- XtVinv%*%(2*dVdrho%*%V0inv%*%dVdrho - d2Vdrho2)%*%t(XtVinv)
#print(dim(HinvdHdrho))
#print(dim(Hinv))
#print(dim(d2Hdrho2))
d2Ldrho2 <- d2Ldrho2  + 0.5*tr.prod(-HinvdHdrho,HinvdHdrho) + 0.5*tr.prod(Hinv,d2Hdrho2) 

d2Ldomega2<-d2Ldomega2 + (p/2)/omega^2


if(!isnullZ2){
HinvdHdtheta1<-Hinv%*%dHdtheta1
HinvdHdtheta2<-Hinv%*%dHdtheta2
dVdtheta1<-Z1
dVdtheta2<-Z2


d2Hdtheta12 <- XtVinv%*%(2*dVdtheta1%*%V0inv%*%dVdtheta1 )%*%t(XtVinv)
d2Hdtheta22 <- XtVinv%*%(2*dVdtheta2%*%V0inv%*%dVdtheta2 )%*%t(XtVinv)
d2Hdtheta1dtheta2 <- XtVinv%*%(2*dVdtheta1%*%V0inv%*%dVdtheta2 )%*%t(XtVinv)
#rho, theta2
d2Hdtheta1drho<- XtVinv%*%(dVdtheta1%*%V0inv%*%dVdrho + dVdrho%*%V0inv%*%dVdtheta1)%*%t(XtVinv)
d2Hdtheta2drho<- XtVinv%*%(dVdtheta2%*%V0inv%*%dVdrho + dVdrho%*%V0inv%*%dVdtheta2)%*%t(XtVinv)
d2Ldrhodtheta <-  d2Ldrhodtheta + 0.5*c(tr.prod(-HinvdHdtheta1,HinvdHdrho) + tr.prod(Hinv,d2Hdtheta1drho),tr.prod(-HinvdHdtheta2,HinvdHdrho) + tr.prod(Hinv,d2Hdtheta2drho))

#theta1,theta2
#print(dim(HinvdHdtheta1))
#print(dim(HinvdHdtheta2))
#print(dim(Hinv))
#print(dim(d2Hdtheta12))
#print(dim(d2Hdtheta22))
#print(dim(d2Hdtheta1dtheta2))

h11<-tr.prod(-HinvdHdtheta1,HinvdHdtheta1) + tr.prod(Hinv,d2Hdtheta12)
h12<-tr.prod(-HinvdHdtheta1,HinvdHdtheta2) + tr.prod(Hinv,d2Hdtheta1dtheta2)
h22<-tr.prod(-HinvdHdtheta2,HinvdHdtheta2) + tr.prod(Hinv,d2Hdtheta22)


d2Ldtheta2<-d2Ldtheta2  + 0.5*matrix(c(h11,h12,h12,h22),2,2)

}else{
#rho,theta1
dVdtheta1<-Z1
d2Hdtheta12<- XtVinv%*%(2*dVdtheta1%*%Vinv%*%dVdtheta1 )%*%t(XtVinv)  #- d2Vdtheta2=0
HinvdHdtheta1<-Hinv%*%dHdtheta1
d2Ldtheta2<-d2Ldtheta2  + 0.5*tr.prod(-HinvdHdtheta1,HinvdHdtheta1) + tr.prod(Hinv,d2Hdtheta12)

d2Hdtheta1drho<- XtVinv%*%(dVdtheta1%*%Vinv%*%dVdrho + dVdrho%*%Vinv%*%dVdtheta1)%*%t(XtVinv)
d2Ldrhodtheta <- d2Ldrhodtheta + 0.5*(tr.prod(-HinvdHdtheta1,HinvdHdrho) + tr.prod(Hinv,d2Hdtheta1drho))
}


}#end REML




Hess1<-  c(d2Ldomega2     , d2Ldomegadbeta,d2Ldomegadrho,d2Ldomegadtheta)
Hess2<-  cbind(d2Ldomegadbeta , d2Ldbeta2     ,d2Ldbetadrho ,d2Ldbetadtheta)
Hess3<-  c(d2Ldomegadrho  , d2Ldbetadrho  ,d2Ldrho2     ,d2Ldrhodtheta)
if(!isnullZ2){
Hess4<-  cbind(d2Ldomegadtheta, t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}else{
Hess4<-  c(d2Ldomegadtheta,t(d2Ldbetadtheta),d2Ldrhodtheta,d2Ldtheta2)
}

Hess<-  rbind(Hess1,Hess2,Hess3,Hess4)
                  
rownames(Hess)<-colnames(Hess)<- par.names


attr(L, "hessian")  <-    Hess
}#end if(hess){

L
#########################################################################
}#end  ARerr2.random
################################################################################


################################################################################
### random intercept model ##################################################
################################################################################


random.intercept<-function(x,Y,X,Z1,Z2=NULL,grad=TRUE,hess=FALSE,pr=TRUE){
#x: first omega, then beta, then rho,  then random effects theta1 and theta2
#Z1 is list of block matrices for random effects (HH)
#Z2 is list of block matrices for random effetcs (area)
hilf<-dim(X)
p<-hilf[2]
p1<-p+1
n<-hilf[1]
#print(p)
#print(n)
isnullZ2<-is.null(Z2)
if(isnullZ2){q<-1}else{q<-2}
if(pr){
beta<-x[1:p]
theta<-abs(x[p1:length(x)])
}else{
omega<-abs(x[1])
beta<-x[2:p1]
theta<-abs(x[(p+2):length(x)]) #theta is of length 2, associated with Z1 and Z2
}
#theta<-pmax(theta,0)
Y<-as.matrix(Y)


diagn<-diag(n) 
Xbeta<-X%*%beta



hilf<- get.Vadd(Z1,Z2,theta)
V<-hilf$V
Vinv<-hilf$Vinv

logdetV<-l.det.V(V)


ry<- as.matrix(Y - Xbeta)
rytVinv<-prod.vector.list.vector(X=ry,Z1=Vinv)
rytVinvry<- sum(t(rytVinv)*ry)


if(pr){
omega<-rytVinvry/n
}
#   terms                    T1                T3               T4
L<- 0.5*n*log(2*pi) + 0.5*n*log(omega)  + 0.5*logdetV + 0.5*rytVinvry/omega            
attr(L, "omega") <-   omega


#now compute gradient
##############################################
if(grad){
XtVinv<-prod.vector.list.vector(X=X,Z1=Vinv)
XtVinvry<-XtVinv%*%ry
C1<- prod.list.list.list(Z1=Vinv,Z2=Z1,Z3=Vinv)
rytC1ry<-prod.vector.list.vector(X=ry,Z1=C1,Z2=NULL,Y=ry)
tr.VinvZ1<-tr.prod.list(V1=Vinv,V2=Z1)

par.names<-c("omega",paste("beta",1:p),paste("theta",1:q))

dLdomega <- (n/2)*(1/omega) - 0.5*rytVinvry/(omega^2)
dLdbeta  <- - XtVinvry
dLdtheta <-  0.5*tr.VinvZ1 -1/(2*omega)*rytC1ry
if(!isnullZ2){
C2<- prod.list.list.list(Z1=Vinv,Z2=Z2,Z3=Vinv)
rytC2ry<-prod.vector.list.vector(X=ry,Z1=C2,Z2=NULL,Y=ry)
tr.VinvZ2<-tr.prod.list(V1=Vinv,V2=Z2)
dLdtheta <- c(dLdtheta,0.5*tr.VinvZ2 -1/(2*omega)*rytC2ry)
}#end if(!isnullZ2){

gr<-c(dLdomega,dLdbeta,dLdtheta)
names(gr)<-par.names
attr(L, "gradient") <-   gr
if(pr){
gr.pr<-c(dLdbeta,dLdtheta)
names(gr.pr)<-par.names[-1]
attr(L, "gradient") <-   gr.pr
}

}#end if grad

#now compute hessian
##############################################
if(hess & grad){

XtVinvX<- XtVinv%*%X

d2Ldomega2<-(n/2)*(1/omega^2)
d2Ldomegadbeta<- rep(0,p)
d2Ldomegadtheta<-(1/(2*omega))*tr.VinvZ1

d2Ldbeta2      <- XtVinvX/omega
d2Ldbetadtheta <- matrix(0,p,q)

#here again expected observation matrix
d2Ldtheta2 <- (-0.5+ (1/omega))*tr.prod.list(V1=C1,V2=Z1)

if(!isnullZ2){
d2Ldomegadtheta<-c(d2Ldomegadtheta,1/(2*omega)*tr.VinvZ2)

h12<-(-0.5+ (1/omega))*tr.prod.list(V1=C1,V2=Z2)
h22<-(-0.5+ (1/omega))*tr.prod.list(V1=C2,V2=Z2)
d2Ldtheta2 <- matrix(c(d2Ldtheta2,h12,h12,h22),2,2,byrow=TRUE)
                                     
}#end if(!isnullZ2){

Hess1<-  c(d2Ldomega2     , d2Ldomegadbeta,d2Ldomegadtheta)
Hess2<-  cbind(d2Ldomegadbeta , d2Ldbeta2 ,d2Ldbetadtheta)

if(!isnullZ2){
Hess3<-  cbind(d2Ldomegadtheta, d2Ldbetadtheta,d2Ldtheta2)
}else{
Hess3<-  c(d2Ldomegadtheta, d2Ldbetadtheta,d2Ldtheta2)
}
Hess<-  rbind(Hess1,Hess2,Hess3)
rownames(Hess)<-colnames(Hess)<- par.names
attr(L, "hessian")  <-    Hess
}#end if(hess){
L
}#end  AR.random
################################################################################
# end random  intercept
################################################################################






#Functions
minus<-function(A){return(-A)}


tr.list<-function(V){
K<-length(V)
traceV<-0
for(k in 1:K){
traceV<-traceV+tr(V[[k]])
}
traceV
}


l.det.V<-function(V){
K<-length(V)
l.detV<-0
for(k in 1:K){
l.detV<-l.detV+2*sum(log(diag(chol(V[[k]]))))
}
l.detV
}

get.Vblock<-function(V){
K<-length(V)
U<-Matrix(0,0,0)
for(k in 1:K){
U<-bdiag(U,V[[k]])
}
U
}

get.Vadd<-function(Z1,Z2=NULL,theta){
K<-length(Z1)
Vinv<-V<-vector("list",K)
if(!is.null(Z2)){
for(k in 1:K){
nk<-dim(Z1[[k]])[1]
V[[k]]<-diag(nk)+theta[1]*Z1[[k]]+theta[2]*Z2[[k]]
Vinv[[k]]<-solve(V[[k]])
}
}else{
for(k in 1:K){
nk<-dim(Z1[[k]])[1]
V[[k]]<-diag(nk)+theta[1]*Z1[[k]]
Vinv[[k]]<-solve(V[[k]])
}
}
return(list(V=V,Vinv=Vinv))
}

prod.vector.list.vector<-function(X,Z1,Z2=NULL,Y=NULL){
K<-length(Z1)
h<-dim(X)
n<-h[1]
p<-h[2]
method<-1  
v<-n
if(!is.null(Z2) & is.null(Y)){method<-2}
if(is.null(Z2) & !is.null(Y)){method<-3;v<-dim(Y)[2]}
if(!is.null(Z2) & !is.null(Y)){method<-4;v<-dim(Y)[2]}
#print(v)
#print(p)
product<-Matrix(0,p,v)
#print(dim(product))
counter<-0;
for(k in 1:K){
nk<-dim(Z1[[k]])[1]
ind<-(counter+1):(counter+nk)
#print(dim(t(X[ind,])))
#print(dim(Z1[[k]]))
#print(dim(Y[ind,]))
#print(ind)
#print(k)
switch(method,
{ product[,ind]<- t(X[ind,])%*%Z1[[k]] },
{ product[,ind]<- t(X[ind,])%*%Z1[[k]]%*%Z2[[k]] },
{ product<-product + t(X[ind,])%*%Z1[[k]]%*%Y[ind,] },
{ product<-product + t(X[ind,])%*%Z1[[k]]%*%Z2[[k]]%*%Y[ind,] }
)
counter<-counter+nk
}#end for
product
}#end  function


prod.list.list.list<-function(Z1,Z2,Z3=NULL){
K<-length(Z1)
product<-vector("list",K)
method<-2
if(is.null(Z3)){method<-1}
for(k in 1:K){
switch(method,
{product[[k]]<-Z1[[k]]%*%Z2[[k]]},
{product[[k]]<-Z1[[k]]%*%Z2[[k]]%*%Z3[[k]]}
)
}#end for
product
}#end  function


tr.prod.list<-function(V1,V2){
#both lists of length K
K<-length(V1)
hilf<-0
for(k in 1:K){
hilf<-hilf+tr.prod(V1[[k]],V2[[k]])
}#end for k
hilf
}#end 


createCounter <- function(value) { function(i) { value <<- value+i} }

