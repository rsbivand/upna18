library(Matrix)
#library(Rsolnp) 
library(numDeriv)
#library(dplyr)
source("ARmodels.sparse.add.fcts.r")
source("ARmodels.r")
#library(sparseinv)



###################################################################################
###########                AR.pop.Exact
###################################################################################
AR.pop.exact.sparse<-function(Ys,W,X,model="SAR",se=TRUE,approx=FALSE,approx.se=TRUE,REML=FALSE,tol=1e-6,rho.limits=c(-1,1),rho0=NULL){
  
#ARmod="SAR"
#ARmod="SARerr"
#ARmod="CAR"  
if(model=="SAR"){ARmod<-1}
if(model=="SARerr"){ARmod<-2}
if(model=="CAR"){ARmod<-3}
if(model=="SARAR"){ARmod<-4}
  
ns<-length(Ys)
n<-dim(W)[1]
nr<-n-ns
p<-dim(X)[2]
p1<-p+1

rownames(W)<-colnames(W)<-1:n
rownames(X)<-1:n
W<- as(W, "CsparseMatrix")

chol.list<-list()
chol.list$n<-n
if(isSymmetric(W)){
  
  W<-forceSymmetric(W)
  
  chol.list$Wsimn<--W
  chol.list$Wsimp<-W
  chol.list$Cn <- Cholesky(chol.list$Wsimn, super = TRUE, Imult = 2) #needed for positive rho's
  chol.list$Cp <- Cholesky(chol.list$Wsimp, super = TRUE, Imult = 2) #needed for negative rho's 
  
  chol.list$sim<-TRUE
  
}else{
  
  Whilf<-1*(W>0)
  rSW<-rowSums(W)
  
  if( isTRUE(all.equal(min(rSW),max(rSW))) & isTRUE(all.equal(min(rSW),1)) && isSymmetric(Whilf) ){
    
    cat("similar\n")
    d<-sqrt(rowSums(Whilf))
    d[d==0]<-1
    Whilf<-Whilf/d
    Whilf<-t(Whilf)/d
    Whilf<-forceSymmetric(Whilf)
    
    chol.list$Wsimn<--Whilf
    chol.list$Wsimp<-Whilf
    chol.list$Cn <- Cholesky(chol.list$Wsimn, super = TRUE, Imult = 2) #needed for positive rho's
    chol.list$Cp <- Cholesky(chol.list$Wsimp, super = TRUE, Imult = 2) #needed for negative rho's 
    
    chol.list$sim<-TRUE
    
  }else{

    cat("not similar\n")
    chol.list$Cp<-chol.list$Cn<-NULL
    chol.list$Wsimn<-NULL
    chol.list$Wsimp<-NULL
    chol.list$sim<-FALSE
    
  }#end if/else rowsum==1
  
}#end else if symm


WplusWt<-t(W)+W  
WtW<-t(W)%*%W
#W<-bdiag(W,W,W,W,W,W,W,W,W,W)
#W<-bdiag(W,W,W,W)
  
  
  I<-Diagonal(n,rep(1,n))
  
  indNs<-1:ns
  indNr<-(ns+1):n
  
  #Is<-Matrix(0,ns,ns)
  #diag(Is)<-1
  #Ir<-Matrix(0,n-ns,n-ns)
  #diag(Ir)<-1

eps<-tol


if(approx){
  
  
  #K=6
  
  
  A <- WplusWt[indNr,indNr]
  B <- WtW[indNr,indNr]
  #6th order Taylor approximation of Muu.inv
  A2<-A%*%A
  A3<-A%*%A2
  A4<-A%*%A3
  A5<-A%*%A4
  A6<-A%*%A5
  
  B2<-B%*%B
  B3<-B%*%B2
  B4<-B%*%B3
  B5<-B%*%B4
  B6<-B%*%B5
  AB<-A%*%B
  A2B<-A%*%AB
  A3B<-A%*%A2B
  A4B<-A%*%A3B
  A5B<-A%*%A4B
  
  AB2<-A%*%B2
  A2B2<-A%*%AB2
  A3B2<-A%*%A2B2
  A4B2<-A%*%A3B2
  AB3 <- A%*%B3
  A2B3<- A%*%AB3 
  A3B3<-A%*%A2B3
  AB4 <-A%*%B4
  A2B4<-A%*%AB4
  AB5<-A%*%B5
  
  
  
  Muu.inv.list<- list(.symDiagonal(nr),A,A2-B,A3-2*AB,A4-3*A2B+B2,A5-4*A3B+3*AB2,
                      A6-5*A4B+6*A2B2-B3,-6*A5B+10*A3B2-4*AB3,15*A4B2-10*A2B3+B4,-20*A3B3+5*AB4,15*A2B4-B5,-6*AB5,B6)
  
  
  Mnew.list<-list()
  Aus<-WplusWt[indNr,indNs]
  Bus<-WtW[indNr,indNs]
                               
  Mnew.list[[1]]<-  t(Aus)%*%Muu.inv.list[[1]]%*%Aus
  
  hilf<- t(Aus)%*%Muu.inv.list[[1]]%*%Bus
  Mnew.list[[2]]<-  t(Aus)%*%Muu.inv.list[[2]]%*%Aus - (hilf+t(hilf))
  
  for(i in 3:13){
    hilf< - t(Aus)%*%Muu.inv.list[[i-1]]%*%Bus
    Mnew.list[[i]]<- t(Aus)%*%Muu.inv.list[[i]]%*%Aus - (hilf+t(hilf)) + t(Bus)%*%Muu.inv.list[[i-2]]%*%Bus
  }#end
  hilf<- t(Aus)%*%Muu.inv.list[[13]]%*%Bus
  Mnew.list[[14]]<- - (hilf+t(hilf)) + t(Bus)%*%Muu.inv.list[[12]]%*%Bus
  Mnew.list[[15]]<- t(Bus)%*%Muu.inv.list[[13]]%*%Bus  
   
  # test
  # M<- I - rho*WplusWt + rho^2*WtW
  # hilf2 <- M[indNs,indNr]%*%solve(M[indNr,indNr])%*%M[indNr,indNs]
  # hilf3 <- Reduce("+",get.list(Wlists[[1]],rho^seq(2,14))) 
  # summary(as.numeric(hilf2-hilf3))

  
  
  Wlists<-list()
  Wlists[[1]]<-Mnew.list
  if(ARmod==1){
    hilf<-get.W.lists(W,X,K=50,Xonly=TRUE)$WX.list #contains WX,W2X,W3X,etc.
    hilf<-make.sample.list(hilf,indNs,1:p) #contains WX[s,],W2X[s,],W3X[s,],etc.
    WX.list<-list()
    WX.list[[1]]<-X[indNs,]
    WX.list<-append(WX.list,hilf)  #new WX.list is Xs and WX[[1]],WX[[2]], etc
    Wlists[[2]]<-WX.list
    rm(hilf,WX.list)
  }#end
  
}#end approx


# SAR, SARerr, CAR (1 parameter family)
if(ARmod<=3){
counter<<-0
#res<-solnp(pars=0.2, fun=log.lik.fct, LB = -1+eps, UB = 1-eps,ARmod=1,pr=T)


# It does minimise -L (and hence maximise L)

# for whole interval
if(rho.limits[1]*rho.limits[2]>0){
cat("For specific interval\n")
res<-try(optimize(AR.exact.log.lik.fct,interval=c(rho.limits[1],rho.limits[2]),Ys,W,X,WplusWt,WtW,
                   ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I,tol = tol^2))


}else{

cat("For positive values only\n")
#print(rho.limits)
res1<-try(optimize(AR.exact.log.lik.fct,interval=c(0,rho.limits[2]-eps),Ys,W,X,WplusWt,WtW,
ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I,tol = tol^2))


cat("Number Iterations:",counter,"\n")
CleanEnvir(counter)
counter<<-0
cat("For negative values only\n")
res2<-try(optimize(AR.exact.log.lik.fct,interval=c(rho.limits[1]+eps,0),Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I,
                  tol = tol^2))
cat("Number Iterations:",counter,"\n")
CleanEnvir(counter)
if(res1$objective<=res2$objective){res<-res1}else{res<-res2}
}

#print(res1)
#print(res2)



rho<-res$minimum
L<-AR.exact.log.lik.fct(par=rho,Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I)


beta<-as.matrix(attr(L,"beta"))
omega<- as.numeric(attr(L,"omega")) 
Lpos<- -c(L)

par1<-c(omega,beta,rho)

}else{
  
  
  # SRAR model
  
  
  #calculate matrices needed
  
  # (A2*A)'*(A2*A) = I - lambda*(W2plusW2t)- rho*W1plusW1t + lambda^2*W2tW2 + rho^2*W1tW1
  # +rho*lambda*(W2W1+W2tW1+t(W2tW1)+t(W2W1)) - lambda^2*rho*(W2tW2W1 + t(W2tW2W1))
  # -rho^2*lambda(W1tW2W1+t(W1tW2W1))+rho^2*lambda^2*( t(W2W1)W2W1
  
  # when W2=W1=W
  #(A2*A)'*(A2*A) = I - (lambda+rho)*WplusWt + (lambda^2+rho^2)*WtW +rho*lambda*(WW+2*WtW+t(WW))
  # - (lambda^2*rho +rho^2*lambda )*( WtWW +t(WtWW)) + rho^2*lambda^2*(WtWtWW)
  
  
  WW<-W%*%W
  WWplus2WtWplusWtWt<-WW+2*WtW+t(WW)
  WtWW<-  t(W)%*%WW
  WtWWplusWtWtW<- WtWW+t(WtWW)
  WtWtWW<-t(WW)%*%WW
  
  
  Wlists<-list(WW=WW,WWplus2WtWplusWtWt=WWplus2WtWplusWtWt,WtWWplusWtWtW=WtWWplusWtWtW,WtWtWW=WtWtWW)
  
   # source("ARmodels.sparse.r")
  # SARAR 2-parameter family
  if(is.null(rho0)){rho0<-c(0,0)}
  
  counter<<-0
  
  res <- try(optim(par=rho0, fn=AR.exact.log.lik.fct, gr = NULL, 
  # arguments of the AR.exact.log.lik.fct    
  Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I,
         method = "L-BFGS-B",
        lower = c(rho.limits[1]+eps,rho.limits[1]+eps), upper = c(1-eps,1-eps),
        control = list(), hessian = FALSE))
  
  if(res$convergence==0){
    rho<-res$par # Y
    cat("Number Iterations:",counter,"\n")
    CleanEnvir(counter)
    
  }else{return(NULL)}
   
  L<-AR.exact.log.lik.fct(par=rho,Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=TRUE,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I)
  
  
  beta<-as.matrix(attr(L,"beta"))
  omega<- as.numeric(attr(L,"omega")) 
  Lpos<- -c(L)
  
  par1<-c(omega,beta,rho[1],rho[2])
  
  
}



# L<-AR.exact.log.lik.fct(par=0.2,Ys,W,X,WplusWt,WtW,ARmod=ARmod,pr=TRUE,approx=approx,Wlists,ns,n,p,indNs,indNr,I)
#L<-res$objective
#res1<-AR.exact.log.lik.fct(0.1,Ys,W,X,WplusWt,WtW,ARmod=ARmod,pr=TRUE,ns,n,p,indNs,indNr,I)
#print(res1)


#par0<-c(as.numeric(attr(L,"omega"))+1,as.numeric(attr(L,"beta"))-0.01,res$minimum+0.05)
#optim(par=par0,fn=log.lik.fct,pr=F)






if(se){
#calculate first derivative
#print(log.lik.fct(par1,W,Ys,X,WplusWt,WtW,ARmod=ARmod,pr=F))

  
  
cat("Calculation Info Matrix\n")
  
if(ARmod<=3){
par.names<-c("sigma2",paste("beta",1:p),"rho")
}else{   
par.names<-c("sigma2",paste("beta",1:p),"rho","lambda")  
  }
#calculate 2nd derivative




if(approx.se){
 
  # approximate method

if(0){  
grad<-as.matrix(c(jacobian(AR.exact.log.lik.fct,par1,method="Richardson", method.args=list(),Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=F,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I)))

#print(grad)  
  
Info<-hessian(AR.exact.log.lik.fct,par1,method="Richardson", method.args=list(),Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=F,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I)
#print(Info)
#print(par.names)
rownames(Info)<-colnames(Info)<- par.names  
rownames(grad)<-par.names

print(grad)
print(Info)

Cov<-try(solve(Info))
if(inherits(Cov,"try-error")){Cov<-ginv(Info)}

# cat("hello\n")
Cov.omega<-Cov[1,1]
Cov.beta<-Cov[2:p1,2:p1]
Cov.rho<-Cov[p+2,p+2]

}else{
  
 
  
#print(grad)

par2<-c(omega,rho) 
#print(par2)
f1<-function(par2,beta,Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=F,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I){
  AR.exact.log.lik.fct(c(par2[1],beta,par2[2:length(par2)]),Ys,W,X,WplusWt,WtW,ARmod=ARmod,chol.list,pr=F,approx=approx,REML=REML,Wlists,ns,n,p,indNs,indNr,I)
  }


#print(par2)

Info<-hessian(f1,x=par2,method="Richardson", method.args=list(),beta=beta,Ys=Ys,W=W,X=X,WplusWt=WplusWt,WtW=WtW,ARmod=ARmod,chol.list=chol.list,pr=F,approx=approx,REML=REML,Wlists=Wlists,ns=ns,n=n,p=p,indNs=indNs,indNr=indNr,I=I)
#print(Info)
hilf<-jacobian(f1,x=par2,method="Richardson", method.args=list(),beta=beta,Ys=Ys,W=W,X=X,WplusWt=WplusWt,WtW=WtW,ARmod=ARmod,chol.list=chol.list,pr=F,approx=approx,REML=REML,Wlists=Wlists,ns=ns,n=n,p=p,indNs=indNs,indNr=indNr,I=I)

#print(hilf)
grad<-as.matrix(c(hilf))


#print(Info)
#obtain these numerically  
if(ARmod==1){
  A  <- I - rho[1]*W
  M=I-rho*WplusWt+rho^2*WtW 
  BsX<-solve(A,X)[indNs,]
}
if(ARmod==2){
  M=I-rho*WplusWt+rho^2*WtW 
  BsX<-X[indNs,]
}
if(ARmod==3){
  M<-I - rho*W 
  BsX<-X[indNs,]
}
if(ARmod==4){
  
  A  <- I - rho[1]*W  # Y process
  A2 <- I - rho[2]*W  # residual process
  BsX<-solve(A,X)[indNs,]
  
  M<-Diagonal(n)- (rho[1]+rho[2])*WplusWt + (rho[1]^2+rho[2]^2)*WtW + rho[1]*rho[2]*WWplus2WtWplusWtWt
  M<-M - (rho[1]^2*rho[2]+rho[1]*rho[2]^2)*WtWWplusWtWtW + rho[1]^2*rho[2]^2*WtWtWW
  
  
}
   #Bs%*%X




#print(isSymmetric(M))
Mss<-M[indNs,indNs]
Msr<-M[indNs,indNr]
Mrr<-M[indNr,indNr]
Mrs<-M[indNr,indNs]
Mrr.inv<-solve(Mrr)

Bss <- Mss - forceSymmetric(Msr%*%Mrr.inv%*%Mrs) #solve(Ainv.ss )

#print(dim(Bss))
Info.beta<- (t(BsX)%*%Bss%*%BsX)/omega

Info.beta.omega<-rep(0,p)





if(ARmod==2 | ARmod==3){
  Info.beta.rho<-rep(0,p)
  
  
}else{
  
  Xbeta<-X%*%beta
  f2<-function(rho,I,W,Xbeta){A<-I -rho[1]*W;return(solve(A,Xbeta))}
  
  
  #drytdrho<-  - t(dBsXdrho%*%Xbeta) #=dBsdrho*drytdBs and  drytdBs=- X*beta
  
  drytdrho <-  - t(f2(rho[1]+tol,I,W,Xbeta)-f2(rho[1]-tol,I,W,Xbeta))/(2*tol)
  
  #print(dim(drytdrho))
  #print(dim(Bss))
  #print(dim(BsX))
  
  Info.beta.rho<-   as.numeric(- t(drytdrho[,indNs]%*%Bss%*%BsX)/omega)
  
  Info.beta.lambda <- rep(0,p)
  
}



if(ARmod<=3){ 

#c(omega,beta,rho) 
Info1<-c(Info[1,1],Info.beta.omega,Info[1,2])
Info2<-cbind(Info.beta.omega,Info.beta,Info.beta.rho)
Info3<-c(Info[1,2],Info.beta.rho,Info[2,2])
Info<-rbind(Info1,Info2,Info3)
  
rownames(Info)<-colnames(Info)<- par.names  


Cov<-try(solve(Info))
if(inherits(Cov,"try-error")){Cov<-ginv(Info)}

Cov.omega<-Cov[1,1]
Cov.beta<-Cov[2:p1,2:p1]
Cov.rho<-Cov[p+2,p+2]
Cov.lambda<-lambda.se<-NULL
rho.se<-sqrt(Cov.rho)
omega.se<-sqrt(Cov.omega)


rs<- Ys - BsX%*%beta
grad.beta <- (t(BsX)%*%Bss%*%rs)/omega 
grad<- c(grad[1],as.numeric(grad.beta),grad[2])
names(grad)<-par.names

}else{
  
  #c(omega,beta,rho) 
  Info1<-c(Info[1,1],Info.beta.omega,Info[1,2],Info[1,3])
  Info2<-cbind(Info.beta.omega,Info.beta,Info.beta.rho,Info.beta.lambda)
  Info3<-c(Info[1,2],Info.beta.rho,Info[2,2],Info[2,3])
  Info4<-c(Info[1,3],Info.beta.lambda,Info[2,3],Info[3,3])
  Info<-rbind(Info1,Info2,Info3,Info4)
  
  rownames(Info)<-colnames(Info)<- par.names  
  
  
  Cov<-try(solve(Info))
  if(inherits(Cov,"try-error")){Cov<-ginv(Info)}
  
  Cov.omega<-Cov[1,1]
  Cov.beta<-Cov[2:p1,2:p1]
  Cov.rho<-Cov[(p+2):(p+3),(p+2):(p+3)]
  rho.se<-sqrt(diag(Cov.rho))
  omega.se<-sqrt(Cov.omega)
  names(rho)<-names(rho.se)<-c("rho","lambda")
  
  
  rs<- Ys - BsX%*%beta
  grad.beta <- as.numeric((t(BsX)%*%Bss%*%rs)/omega) 
  grad<- c(grad[1],grad.beta,grad[2:3])
  
  #print(grad)
  #print(par.names)
  names(grad)<-par.names
}





}# end if(0/1)

# exact
}else{ # if approx.se
  

hilf<-try(AR.exact.Info(Ys,W,X,beta,rho,sigma2=omega,model=model,approx))
#cat("hello") 
if(!inherits(hilf,"try-error")){
  
Info<-hilf$Info
grad<-hilf$grad

# return(list(grad=grad,Cov=Cov,Cov.sigma2=Cov.omega,Cov.beta=Cov.beta,Cov.rho=Cov.rho,Info=Info))
  
  
rownames(Info)<-colnames(Info)<- par.names  
rownames(grad)<-par.names

Cov<-hilf$Cov
Cov.omega<-Cov[1,1]
Cov.beta<-Cov[2:p1,2:p1]
Cov.rho<-Cov[p+2,p+2]

rho.se<-sqrt(Cov.rho)
omega.se<-sqrt(Cov.omega)

}else{
  
  grad<-Info<-Cov<-Cov.rho<-Cov.omega<-NA
  Cov.beta<-NA*diag(p)
  
  
}# if(!inherits(hilf,"try-error")){


}# if approx.se

if(p>1){beta.se<-sqrt(diag(Cov.beta))}else{beta.se<-sqrt(Cov.beta)}

return(list(L=Lpos,grad=grad,sigma2=omega,beta=beta,rho=rho,Cov=Cov,Cov.sigma2=Cov.omega,Cov.beta=Cov.beta,Cov.rho=Cov.rho,Info=Info,
            beta.se=beta.se,rho.se=rho.se,sigma2.se=omega.se,res=res))


}#if se

return(list(L=Lpos,sigma2=omega,beta=beta,rho=rho,res=res))

}#end function AR.pop.exact




#############################################################
#             beginning of log-lik function
#############################################################
AR.exact.log.lik.fct<-function(par,Ys,W,X,WplusWt,WtW,ARmod=1,chol.list=NULL,pr=TRUE,approx=F,REML=F,Wlists=NULL,ns=length(Ys),n=dim(X)[1],p=dim(X)[2],indNs=1:ns,indNr=(ns+1):n,I=NULL){

  #cat("hello")
  if(exists("counter")){counter<<-counter+1;cat("Iteration:",counter,"\n")}
  
  
  if(is.null(I)){
    
    #I<-.symDiagonal(n)
    I<-Diagonal(n)
  }
  
  if(pr){
    rho<-par # rho[1] is SAM parameter, rho[2] is SEM parameter
    
    #only rho
  }else{
    p1<-p+1
    omega<-abs(par[1]) #first omega
    beta<-par[2:p1]  #2nd beta
    rho<-par[(p+2):length(par)]   #last rho,lambda
  }
  #print(dim(A))
  #print(dim(W))
  #print(rho)
  
  
  
  
  
  switch(ARmod,
{
  # SAR model
  # (Ainv)ss= (Ass - Asr*Arrinv*Ars)^{-1}
  #therefore   (Ainv)ss^{-1}=Ass - Asr*Arrinv*Ars
 
   
  if(approx){
    
    
    rho.seq1<-rho^seq(2,16)  #first to MsuMuuinvMus
    
    Msr.Mrr.inv.Mrs<-forceSymmetric(Reduce("+",get.list(Wlists[[1]],rho.seq1)))  
    
    rho.seq2<-rho^seq(0,50) #second sequence refers to WX.list
    BsX<-Reduce("+",get.list(Wlists[[2]],rho.seq2))  
    
    #A<-I-rho*W
    #BsX<-solve(A,X)[indNs,]
    
   
    Mss <- forceSymmetric(I[indNs,indNs] - rho*WplusWt[indNs,indNs]+rho^2*WtW[indNs,indNs])   #M=AtA
   
    Vss<- Mss - Msr.Mrr.inv.Mrs
     
        
  }else{
    
    A<-I-rho*W
    
    M<-forceSymmetric(I-rho*WplusWt+rho^2*WtW)
    
    
    #Ass <- A[indNs,indNs]
    #Asr <- A[indNs,indNr]
    #Ars <- A[indNr,indNs]
    #Arr <- A[indNr,indNr]
  
  #Arr.inv<-solve(Arr)
  #Ass.inv<-solve(Ass)
  #Ainv.ss<-solve(Ass-Asr%*%Arr.inv%*%Ars)
  
  #similar for (Ainv)rr
  #Ainv.rr<-solve(Arr-Ars%*%Ass.inv%*%Asr)
  # (Ainv)sr=  - Ass.inv * Asr * (Ainv)rr 
  #Ainv.sr<- - Ass.inv%*%Asr%*%Ainv.rr
  
  #Bs<- cBind(Ainv.ss,Ainv.sr)
  BsX<-solve(A,X)[indNs,]   #Bs%*%X

  constant<-2
  
  }#end if approx
  
},
{ # SARerr model
  #print(indNs)
  #print(X)
  BsX<-X[indNs,]
  
  
  
  if(approx){
    
    rho.seq1<-rho^seq(2,16)
    
    Msr.Mrr.inv.Mrs<-forceSymmetric(Reduce("+",get.list(Wlists[[1]],rho.seq1)))  
    
    
    
    Mss<-forceSymmetric(I[indNs,indNs]-rho*WplusWt[indNs,indNs]+rho^2*WtW[indNs,indNs])   #M=AtA
    
    Vss<- Mss - Msr.Mrr.inv.Mrs
  }else{
  
  #WtplusW<-t(W)+W
  #WtW<-t(W)%*%W
  M<- forceSymmetric(I-rho*WplusWt+rho^2*WtW)
  
  
  
  A<-I-rho*W
  }
  constant<-2
}
, 
{ #CAR model
  approx<-F
  BsX<-X[indNs,]
  M<-A<- I - rho*W
  constant<-1
  
}
,
{ #SARAR model
  approx<-F
  
  A<- I - rho[1]*W # Y process
  #Bs<- cBind(Ainv.ss,Ainv.sr)
  BsX<-solve(A,X)[indNs,]   #Bs%*%X
  
  A2<- I - rho[2]*W  # residual process
  
  # V = Ainv*A2inv*t(A2inv)*t(Ainv) = ( t(A)*t(A2)*A2*A = t(A2*A)*(A2*A) )
  
  if(is.null(Wlists)){
   M <- A2%*%A
   M <- t(M)%*%M
  }else{
  M<-Diagonal(n)- (rho[1]+rho[2])*WplusWt + (rho[1]^2+rho[2]^2)*WtW + rho[1]*rho[2]*Wlists$WWplus2WtWplusWtWt
  M<-M - (rho[1]^2*rho[2]+rho[1]*rho[2]^2)*Wlists$WtWWplusWtWtW + rho[1]^2*rho[2]^2*Wlists$WtWtWW
  }
  M<-forceSymmetric(M)
  #print(rho)
  
  constant<-2
  
}
  )

#print(isSymmetric(M))

#Vss<-M[indNs,indNs]-M[indNs,indNr]%*%solve(M[indNr,indNr])%*%M[indNr,indNs]
#rs<- Ys - BsX%*%beta

if(approx){
if(pr){
  
  #print(dim(Vss))
  #print(dim(Ys))
  #print(length(Ys))
  #print(dim(BsX))
  beta<-solve(t(BsX)%*%Vss%*%BsX)%*%t(BsX)%*%Vss%*%Ys
  mu<-BsX%*%beta
  ry<-Ys-mu
  Vssry<-Vss%*%ry
  rstVssrs<-as.numeric(t(ry)%*%Vssry)
  omega<-c(rstVssrs/ns)
  
}else{

  mu<-BsX%*%beta
  ry<-Ys-mu
  Vssry<-Vss%*%ry
  rstVssrs<-as.numeric(t(ry)%*%Vssry)   

}

logdetVss<- log(det(Vss))         #sum(log(diag(chol(Vss))))

# not approx, but exact
}else{
  if(pr){
  
  #print(Muu)  
    
  Muu<-M[indNr,indNr]
  
  #cat("hello1\n")
  
  chol.Muu<-chol(Muu)
  
  #cat("hello2\n")
  
  chol.Mss<-chol(M[indNs,indNs])
  #a1<- A[,indNs]%*%BsX  
  a1<- chol.Mss%*%BsX
  a2<- solve(t(chol.Muu),M[indNr,indNs]%*%BsX)
  # a3<- t(M[indNr,indNs]%*%BsX)%*%solve(Muu,M[indNr,indNs]%*%BsX)
  # a3-t(a2)%*%a2
  
  #b1<- A[,indNs]%*%Ys
  b1<- chol.Mss%*%Ys
  b2<- solve(t(chol.Muu),M[indNr,indNs]%*%Ys)
  
  XstVssXs<-t(a1)%*%a1-t(a2)%*%a2
  XstVssys <-t(a1)%*%b1-t(a2)%*%b2
  
  beta<- solve(XstVssXs)%*%XstVssys
   
  #print(beta)
  ystVssys<-  t(b1)%*%b1 - t(b2)%*%b2
  
  #print(ystVssys)
                  #          ystVssBXsbeta
  rstVssrs<- as.numeric(ystVssys  -2*(t(beta)%*%XstVssys) + t(beta)%*%XstVssXs%*%beta)
   
  #print(rstVssrs)
  omega<-rstVssrs/ns
   
  #print(omega)
  }else{
    
   
    
    Muu<-M[indNr,indNr]
    
    
    #print(Muu)
    
    
    chol.Muu<-chol(Muu)
    
    
    
    chol.Mss<-chol(M[indNs,indNs])
    
    
    
    #print(rho)
    #a1<- A[,indNs]%*%BsX  
    a1<- chol.Mss%*%BsX
    a2<- solve(t(chol.Muu),M[indNr,indNs]%*%BsX)
    # a3<- t(M[indNr,indNs]%*%BsX)%*%solve(Muu,M[indNr,indNs]%*%BsX)
    # a3-t(a2)%*%a2
    
    #b1<- A[,indNs]%*%Ys
    b1<- chol.Mss%*%Ys
    b2<- solve(t(chol.Muu),M[indNr,indNs]%*%Ys)
    
    XstVssXs<-t(a1)%*%a1-t(a2)%*%a2
    XstVssys <-t(a1)%*%b1-t(a2)%*%b2
    
    #beta<- solve(XstVssXs)%*%XstVssys
    
    #print(beta)
    ystVssys<-  t(b1)%*%b1 - t(b2)%*%b2
    
    #print(ystVssys)
    #          ystVssBXsbeta
    rstVssrs<- as.numeric(ystVssys  -2*(t(beta)%*%XstVssys) + t(beta)%*%XstVssXs%*%beta)
    
    
    
    #print(rstVssrs)
    #omega<-rstVssrs/ns
  }
  
#print(chol.list)  
#print(ldetA(rho,chol.list))
#print(prod(diag(chol(Muu))))
#factor 2 is there because of cholesky

  #print(!is.null(chol.list) && !is.null(chol.list$sim) && chol.list$sim)
  
if(!is.null(chol.list) && !is.null(chol.list$sim) && chol.list$sim){
  
  logdetVss<- constant*ldetA(rho[1],chol.list) - 2*sum(log(diag(chol.Muu))) 
  if(ARmod==4){ logdetVss<- logdetVss + constant*ldetA(rho[2],chol.list)  }
  
}else{
  
  
  logdetVss<- 2*sum(log(diag(chol(M))))- 2*sum(log(diag(chol.Muu)))    

  
  }  

}#end if appprox



#print(rstVssrs)
#print(logdetVss)

#   terms                    T1               T2           T3               T4

#it is minus log(detV) to avoid the inverse
#print(rho)

L<- c(0.5*ns*log(2*pi) + 0.5*ns*log(omega)   - 0.5*logdetVss + 0.5*rstVssrs/omega)            

if(REML){
  L<-L+ 0.5*2*sum(log(diag(chol(XstVssXs))))
}
#print(L)


attr(L, "omega") <-   omega
attr(L, "beta") <- as.matrix(beta)

#print(L)
return(L)

#############################################################
} #             end of log-lik function of AR exact
#############################################################


#################################################################
# start computation of info matrix
AR.exact.Info<-function(Ys,W,X,beta,rho,sigma2,model="SAR",approx=F){  
omega<-sigma2  
  

#approx not implemented yet

    #ARmod="SAR"
    #ARmod="SARerr"
    #ARmod="CAR"  
    if(model=="SAR"){ARmod<-1}
    if(model=="SARerr"){ARmod<-2}
    if(model=="CAR"){ARmod<-3}
    if(model=="SARAR"){ARmod<-4}
  
  W<-as(W,"dgCMatrix")
  
  WplusWt<-forceSymmetric(W+t(W))
  WtW<-forceSymmetric(t(W)%*%W)
  #W<-bdiag(W,W,W,W,W,W,W,W,W,W)
  #W<-bdiag(W,W,W,W)
  
  ns<-length(Ys)
  hilf<-dim(X)
  n<-hilf[1]
  p<-hilf[2]
  p1<-p+1
  I<-Diagonal(n,rep(1,n))
  
  indNs<-1:ns
  indNr<-(ns+1):n
  
  
A<-I-rho[1]*W

Ass <- A[indNs,indNs]
Asr <- A[indNs,indNr]
Ars <- A[indNr,indNs]
Arr <- A[indNr,indNr]


switch(ARmod,
{
  # SAR model
  # (Ainv)ss= (Ass - Asr*Arrinv*Ars)^{-1}
  #therefore   (Ainv)ss^{-1}=Ass - Asr*Arrinv*Ars
  
  
  BsX<-solve(A,X)[indNs,]
  
  
  #takes longer 
  #  Vss<- solve(Ainv.ss%*%t(Ainv.ss)+Ainv.sr%*%t(Ainv.sr))
  
  
  M<-I-rho*WplusWt+rho^2*WtW   #M=AtA
  Mss<-M[indNs,indNs]
  Msr<-M[indNs,indNr]
  Mrs<-M[indNr,indNs]
  Mrr<-M[indNr,indNr]
  
  
    Mrr.inv<-solve(Mrr)
    Msr.Mrr.inv.Mrs<-forceSymmetric(Msr %*% Mrr.inv %*% Mrs)
   
    
  Vss<- Mss - Msr.Mrr.inv.Mrs 
  
  
},
{ # SARerr model
  
  BsX<-X[indNs,]

  M<-I-rho*WplusWt+rho^2*WtW 
  Mss<-M[indNs,indNs]
  Msr<-M[indNs,indNr]
  Mrs<-M[indNr,indNs]
  Mrr<-M[indNr,indNr]
  
  if(is(Mrr,"dgCMatrix")){Mrr.inv<-solve(Mrr,sparse=TRUE)}else{Mrr.inv<- solve(Mrr,sparse=FALSE)}
    
    Msr.Mrr.inv.Mrs<-forceSymmetric(Msr %*% Mrr.inv %*% Mrs) 
  
  #cat("hello\n")
  #print(dim(Msr.Mrr.inv.Mrs))
  Vss<- Mss - Msr.Mrr.inv.Mrs 
  
}
, 
{ #CAR model
  BsX<-X[indNs,]
  
  Mss<- Ass
  Msr<-Asr
  Mrs<-Ars
  Mrr<-Arr
  Mrr.inv <- solve(Arr)
  Mss.inv <- solve(Ass)
  #Ainv.ss <- solve(Ass-Asr%*%Arr.inv%*%Ars)
  Vss <- Mss - forceSymmetric(Msr%*%Mrr.inv%*%Mrs) #solve(Ainv.ss )
  
}
,
{ # SARAR
  
  BsX<-solve(A,X)[indNs,]
  
  A2<-I - rho[2]*W
  
  M <- A2%*%A
  M<- t(M)%*%M
  
  Mss<-M[indNs,indNs]
  Msr<-M[indNs,indNr]
  Mrs<-M[indNr,indNs]
  Mrr<-M[indNr,indNr]
  
  
  Mrr.inv<-solve(Mrr)
  Msr.Mrr.inv.Mrs<-forceSymmetric(Msr %*% Mrr.inv %*% Mrs)
  
  
  Vss<- Mss - Msr.Mrr.inv.Mrs   
  
  
}

)
      
     
      BsXbeta<-BsX%*%beta
      ry<-Ys - BsXbeta  #for grad
      Vssry<- Vss%*%ry  #for grad
      rytVssry<-t(ry)%*%Vssry
     
      par.names<-c("sigma2",paste("beta",1:p),"rho")
      


  if(ARmod<=2){
  dMdrho<- -WplusWt+2*rho*WtW
  }
  if(ARmod==3){    
    #CAR
    dMdrho<-    -W
  }
  if(ARmod==4){
    
   dMdrho <- -2*rho[1]%*%t(W)%*%(A2tA2)%*%A1   
    
  }

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
    
    if(is(A,"dgCMatrix")){Ainv<- solve(A,sparse=TRUE)}else{Ainv<- solve(A,sparse=FALSE)}
    
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
  
  
  dVssdrho <- dMssdrho - dMsrdrho%*%Mrr.invMrs + MsrMrr.inv%*%(dMrrdrho%*%Mrr.invMrs - dMrsdrho)
  rytdVssdrhory<-t(ry)%*%dVssdrho%*%ry
  
  #cat("drytdrho%*%Vssry/omega\n")
  #print(dLdrho)
  
  dLdrho<- dLdrho -0.5*tr.prod(Vssinv,dVssdrho) + 0.5*rytdVssdrhory/omega 
  
  
  dLdbeta<- - BsXtVssry/omega
  
 
  dLdomega<- 0.5*ns/omega - 0.5/(omega^2)*rytVssry
  grad<-as.matrix(rBind(dLdomega,dLdbeta,dLdrho))
  rownames(grad)<-par.names

     
      
      VssinvdVssdrho<- Vssinv%*%dVssdrho
      #Mrr.invdMrrdrho<- Mrr.inv%*%dMrrdrho
      #Mrr.invdMrrdrhoMrr.inv<- Mrr.invdMrrdrho%*%Mrr.inv
      #MsrMrr.inv<-Msr%*%Mrr.inv
      #Mrr.invMrs<-Mrr.inv%*%Mrs
      #dMrrdrhoMrr.invdMrrdrho<- dMrrdrho%*%Mrr.invdMrrdrho
      
      #d2Mdrho2<- 2*WtW
      #d2Mssdrho2<-  d2Mdrho2[indNs,indNs]
      #d2Msrdrho2<-  d2Mdrho2[indNs,indNr]
      #d2Mrsdrho2<-  d2Mdrho2[indNr,indNs]
      #d2Mrrdrho2<-  d2Mdrho2[indNr,indNr]
      
      #d2Vssdrho2<- d2Mssdrho2 - d2Msrdrho2%*%Mrr.invMrs + 2*dMsrdrho%*%Mrr.invdMrrdrhoMrr.inv%*%Mrs
      #d2Vssdrho2<-d2Vssdrho2 + 2*Msr%*%Mrr.invdMrrdrhoMrr.inv%*%dMrsdrho - dMsrdrho%*%Mrr.inv%*%dMrsdrho - MsrMrr.inv%*%d2Mrsdrho2 
      #d2Vssdrho2<-d2Vssdrho2 - 2*MsrMrr.inv%*%Mrr.invdMrrdrhoMrr.inv%*%Mrr.invMrs + MsrMrr.inv%*%d2Mrrdrho2%*%Mrr.invMrs
      
      #   "omega"  "beta "   "rho"  
      #omega/omega
      V11 <- 0.5*ns/omega^2
      #beta/beta
      V22 <- t(BsX)%*%Vss%*%BsX/omega
      #first term -1/2*tr(Vssinv*dVssdrho*Vssinv*dVssdrho) - 1/2*tr(Vssin*d2Vssdrho2) +1/2*tr(d2Vssdrho2*Vssinv)/omega + drytdrho*Vss*drydrho/omega
      #rho/rho
      V33 <-  0.5*tr.prod(VssinvdVssdrho,VssinvdVssdrho) 
      
      #omega and beta
      V12 <- matrix(0,1,p)
      #omega and rho
      V13<- - 0.5*tr(VssinvdVssdrho)/omega
      #beta and rho
      V23<-matrix(0,p,1)
      
      if(ARmod==1){
        #additional term because for AR model mean depends on rho
        #rho/rho
        V33 <- V33  + drytdrho%*%Vss%*%t(drytdrho)/omega 
        #beta and rho
        V23 <- - t(drytdrho%*%Vss%*%BsX)/omega
      }
      
      
      
      Info<-rBind(cBind(V11,V12,V13),cBind(t(V12),V22,V23),cBind(t(V13),t(V23),V33))
      Info<-as.matrix(Info)
      
      
      
      
    
  
  
  #print(Info)
  #print(grad)
  #print(par.names)
  rownames(Info)<-colnames(Info)<- par.names  
  rownames(grad)<-par.names
  
  
  Cov<-solve(Info)
  Cov.omega<-Cov[1,1]
  Cov.beta<-Cov[2:p1,2:p1]
  Cov.rho<-Cov[p+2,p+2]
  

return(list(grad=grad,Cov=Cov,Cov.sigma2=Cov.omega,Cov.beta=Cov.beta,Cov.rho=Cov.rho,Info=Info))
  
}#END FUNCTION AR.EXACT.INFO





############################################
#LeSage and Page 2004
###########################################
EM.algorithm.LesagePage<-function(Ys,W,X,model,K=50,tol=1e-6,iter.max=1e2,approx.mean=TRUE){


  
ns<-length(Ys)
hilf<-dim(X)
p<-hilf[2]
n<-hilf[1]
indNs<-1:ns
indNr<-(ns+1):n
rownames(X)<-rownames(W)<-1:n

Xs<-X[indNs,]
Xr<-X[indNr,]

rlm<-lm(Ys~Xs-1)  
beta<-rlm$coefficients
omega<-summary(rlm)$sigma^2
rho<-0  
#EYu/Ys
Y<-c(Ys,Xr%*%beta)


listw<-mat2listw(forceSymmetric(1*(W>0)),style="W")


#check whether symmetric
#if(sum(abs(W-t(W)))<1e-6){meth="Matrix_J"}else{meth="eigen"}
meth<-"Matrix"

WtW<-t(W)%*%W
WplusWt<-W+t(W)


#calculate approximate values
rho.seq<-rep(rho,K+1)^(0:K)
WX.list<-list()
WX.list[[1]]<-X
for(i in 2:(K+1)){
  WX.list[[i]]<-W%*%WX.list[[i-1]]
}#end

iter<-0
diff<-1

I<-.symDiagonal(n)

while(diff>tol & iter<iter.max){
iter<-iter+1

beta0<-beta
rho0<-rho
omega0<-omega



  if(model=="SAR"){
    
    hilf <- try(lagsarlm(Y~X, listw=listw,method=meth))
    
    beta<-as.matrix(hilf$coefficients)
    rho<-hilf$rho
    omega<-hilf$s2
  
    #print(beta)
    #print(rho)
    #print(omega)
    
    Ly<-hilf$LL  
  
  
  
  #print(length(WX.list))  
  #print(rho.seq)
  
  
  if(approx.mean){
    
    
    M<- I - rho*WplusWt+rho^2*WtW
    C<- I- M/diag(M)
    rho.seq<-rep(rho,K+1)^(0:K)
    BX<-Reduce("+",get.list(WX.list,rho.seq))
    #print(BX[indNs,])
    
    Eyr<-as.numeric(BX[indNr,]%*%beta +  C[indNr,indNs]%*%(Ys-BX[indNs,]%*%beta))
    
    
    
    
    
    
  }else{
    # exact
    A<-I-rho*W
    M<- I - rho*WplusWt+rho^2*WtW
    C2<- - solve(M[indNr,indNr],M[indNr,indNs],sparse=TRUE)
      
    AinvX<-solve(A,X)
    #print(ys)
    #print(AinvX[indNs,]%*%beta)
    Eyr<-as.numeric( AinvX[indNr,]%*%beta +  C2%*%(Ys-AinvX[indNs,]%*%beta))
  }
  

  Y<-c(Ys,Eyr)
  
 # Lyr.ys<-dmvnorm(Yr, mean =, sigma = diag(p), log = TRUE)
  
  } #end if 

if(model=="SARerr"){
  
  hilf <- try(errorsarlm(Y~X, listw=listw,method=meth))
  
  beta<-as.matrix(hilf$coefficients)
  rho<-hilf$lambda
  omega<-hilf$s2
  
  Ly<-hilf$LL
  
   
  if(approx.mean){
    
    M<- I - rho*WplusWt+rho^2*WtW
    C<-I- M/diag(M)
    Eyr<-as.numeric(Xr%*%beta +  C[indNr,indNs]%*%(Ys-Xs%*%beta))
  
    }else{
    
    M<- I - rho*WplusWt+rho^2*WtW
    C2<- - solve(M[indNr,indNr],M[indNr,indNs],sparse=TRUE)
    Eyr<-as.numeric(Xr%*%beta +  C2%*%(Ys-Xs%*%beta))
  }
  
  
  Y<-c(Ys,Eyr)
  
  
  } #end if  
  
if(model=="CAR"){
  
} #end if  

if(iter%%10==0){
  cat("Iteration: ",iter,"\n")
cat("beta",beta,"\n")
cat("rho",rho,"\n")
cat("omega",omega,"\n")
cat("complete log-lik:",Ly,"\n")

diff<-abs(rho-rho0)+abs(omega-omega0)+sum(abs(beta-beta0))
cat("diff",diff,"\n")

}


if(abs(omega)> (1/tol)^2 | abs(rho)> (1/tol)^2 | sum(abs(beta))> (1/tol)^2){
  conv<-FALSE;return(list(beta=beta,rho=rho,sigma2=omega,conv=conv)) 
}#end if

}#end while



if(diff>tol | iter>=iter.max){conv<-FALSE}else{conv<-TRUE}

if(model=="SAR"){
  
 #in fact r is A*r
 rs<- (Y - rho*W%*%Y -X%*%beta)[indNs,]  # A*r=A* (Y - Ainv*X*beta)=A*Y-X*beta=Y-rho*W*Y-X*beta
 omega.LP<-t(rs)%*%rs/ns
 #print(rs)
 #print(omega.LP)
}

if(model=="SARerr"){
 rs<-Y-X%*%beta
 rs<-(rs-rho*W%*%rs)[indNs,]
 omega.LP<-t(rs)%*%rs/ns
 #print(rs)
 #print(omega.LP)
}

return(list(beta=beta,rho=rho,sigma2=omega,sigma2.LP=omega.LP,conv=conv,L=Ly,iter=iter)) 

  
}#end EM.Lesage



###########################################################################
#   EM correct
###########################################################################
EM.algorithm<-function(Ys,W,X,model="SARerr",K=50,tol=1e-6,rho0=0,method=2,iter.max=1e2,omega0=NULL,beta0=NULL,approx.mean=1,se=F){
#new EM algorithm  
  #W.sparse<- class(W)=="dgCMatrix"
  ys<-Ys
  ns<-length(ys)
  hilf<-dim(X)
  p<-hilf[2]
  n<-hilf[1]
  
  rownames(W)<-colnames(W)<-1:n
  s<-1:ns
  Xs<-X[s,]
  Xr<-X[-s,]
  
  rlm<-lm(Ys~Xs-1)  
  if(is.null(beta0)){
  
  beta0<-rlm$coefficients
  
  }
  if(is.null(omega0)){
  omega0<-summary(rlm)$sigma^2
  }
  nr<-n-ns
  
  #print(beta0)
  #print(omega0)
  
  #print(nr)
  #print(ns)
  #check symmetry of W
  
  chol.list<-list()
  chol.list$n<-n
  
  if(isSymmetric(W)){
    
  W<-forceSymmetric(W)
  
  chol.list$Wsimn<--W
  chol.list$Wsimp<-W
  chol.list$Cn <- Cholesky(chol.list$Wsimn, super = TRUE, Imult = 2) #needed for positive rho's
  chol.list$Cp <- Cholesky(chol.list$Wsimp, super = TRUE, Imult = 2) #needed for negative rho's 
  
  chol.list$sim<-TRUE
    
  }else{
      Whilf<-1*(W>0)
      rSW<-rowSums(W)
      if(isTRUE(all.equal(min(rSW[rSW!=0]),max(rSW[rSW!=0]),1)) && isSymmetric(Whilf) ){
        
        
        d<-sqrt(rowSums(Whilf))
        d[d==0]<-1
        Whilf<-Whilf/d
        Whilf<-t(Whilf)/d
        Whilf<-forceSymmetric(Whilf)
        
        chol.list$Wsimn<--Whilf
        chol.list$Wsimp<-Whilf
        chol.list$Cn <- Cholesky(chol.list$Wsimn, super = TRUE, Imult = 2) #needed for positive rho's
        chol.list$Cp <- Cholesky(chol.list$Wsimp, super = TRUE, Imult = 2) #needed for negative rho's 
        
        chol.list$sim<-TRUE
        
      }else{
        
        chol.list$Cp<-chol.list$Cn<-NULL
        chol.list$Wsimn<-NULL
        chol.list$Wsimp<-NULL
        chol.list$sim<-FALSE
        
      }#end if/else rowsum==1
      
    }#end else if symm
    
    

  
  WtW<-t(W)%*%W
  WplusWt<-W+t(W)
  
  #Taylor
  if(approx.mean==1 | method==7){
  
  
  
  A <- WplusWt[-s,-s]
  B <- WtW[-s,-s]
  #6th order Taylor approximation of Muu.inv
  A2<-A%*%A
  A3<-A%*%A2
  A4<-A%*%A3
  A5<-A%*%A4
  A6<-A%*%A5
  
  B2<-B%*%B
  B3<-B%*%B2
  B4<-B%*%B3
  B5<-B%*%B4
  B6<-B%*%B5
  AB<-A%*%B
  A2B<-A%*%AB
  A3B<-A%*%A2B
  A4B<-A%*%A3B
  A5B<-A%*%A4B
  
  AB2<-A%*%B2
  A2B2<-A%*%AB2
  A3B2<-A%*%A2B2
  A4B2<-A%*%A3B2
  AB3 <- A%*%B3
  A2B3<- A%*%AB3 
  A3B3<-A%*%A2B3
  AB4 <-A%*%B4
  A2B4<-A%*%AB4
  AB5<-A%*%B5
  
  
  
  Muu.inv.list<- list(.symDiagonal(nr),A,A2-B,A3-2*AB,A4-3*A2B+B2,A5-4*A3B+3*AB2,
                      A6-5*A4B+6*A2B2-B3,-6*A5B+10*A3B2-4*AB3,15*A4B2-10*A2B3+B4,-20*A3B3+5*AB4,15*A2B4-B5,-6*AB5,B6)
  
  
  
  Muu.inv.Mus.list<-list()
  Aus<-WplusWt[-s,s]
  Bus<-WtW[-s,s]
  #rho                               rho^0 
  Muu.inv.Mus.list[[1]]<- - Muu.inv.list[[1]]%*%Aus
  for(i in 2:13){
    Muu.inv.Mus.list[[i]]<- - Muu.inv.list[[i]]%*%Aus + Muu.inv.list[[i-1]] %*%Bus
  }#end
  # rho^14                        rho^10
  Muu.inv.Mus.list[[14]]<- Muu.inv.list[[13]]%*%Bus
    
  
   
  #true inverse
  # Muu <- .symDiagonal(nr) -rho0*WplusWt[-s,-s]+rho0^2*WtW[-s,-s]
  # Muu.inv.true<- solve(Muu)
  # summary(c(as.matrix(Muu.inv.true-Muu.inv)))
  # Mus <-  -rho0*WplusWt[-s,s]+rho0^2*WtW[-s,s] 
  # Muu.inv.Mus.true<-Muu.inv.true%*%Mus 
  # summary(c(as.matrix(Muu.inv.Mus.true-Muu.inv.Mus)))
  
  }#
  
  
  
  #EYu/Ys
  ##################################################################################
  if(model=="SAR"){
    
    Ls<-NULL
    
    #check whether symmetric
    # needed
    # ys,Eus,Vus,X,XtX,WX,XtWX,XtWtWX,WEY
    XtX<-t(X)%*%X
    
    hilf<-get.W.lists(W,X,K=K,Xonly=TRUE)$WX.list
    WX.list<-list()
    WX.list[[1]]<-X
    WX.list<-append(WX.list,hilf)
    #rho0<-0
    #rho0.seq<-rep(rho0,K+1)^(0:K)
    #XtEy,XtWEy,WEy,Ey,XtX,XtWX,XtWtWX,M0rr,omega0
    #XtX,XtWX,XtWtWX,M0rr
    
    WX<-WX.list[[2]]
    XtWX<-t(X)%*%WX
    XtWtWX<-t(WX)%*%WX
    
    
    XtWplusWt<-t(X)%*%WplusWt
    XtWtW<-t(X)%*%WtW
    
    ys<-Matrix(ys,ns,1,dimnames=list(1:ns,NULL))
    
    #these only depend on ys - don't need update
    
    
    
    
    #all terms depend on Eyr need update
    
    #I<-.symDiagonal(n)
    #M<-(I-rho0*WplusWt+rho0^2*WtW)
    
    
    #XtEy,XtWEy,WEy,Ey,XtX,XtWX,XtWtWX,M0rr,omega0
    # check: sum(abs(M/diag(M)-Diagonal(n,1/diag(M))%*%M))
    # update Yu with E(Yu|ys)
    
    
    
  
    #Eyr<-BX[-s,]%*%beta+  C[-s,s]%*%(ys-BX[s,]%*%beta)
    
    #alternative
    #Minv<-solve(M)
    # C1<- Minv[-s,s]%*%solve(Minv[s,s])
    # max(abs(C1-C[-s,s]))
    
  
    rho<-rho0
    beta<-beta0
    omega<-omega0
    
    #maxrhodiff<-0.1
    #print(s)
    
    theta.hat<-matrix(c(rho0,omega0,beta0),1,2+p)
    
    L1<-L2<-NULL
    
    diff<-1
    iter<-0
    while(diff>tol & iter<iter.max){
      
      rho0<-rho
      beta0<-beta
      omega0<-omega
      
      
      iter<-iter+1
      
      
      
      # check: sum(abs(M/diag(M)-Diagonal(n,1/diag(M))%*%M))
      # update Yu with E(Yu|ys)
      
      # XtEy,XtWEy,WEy,Ey,XtX,XtWX,XtWtWX,M0rr,omega0
      
      ##################################################################
      
        
        I<-.symDiagonal(n)
        M<-(I-rho*WplusWt+rho^2*WtW)  
        
        rho.seq<-rep(rho,K+1)^(0:K)
        
        if(approx.mean){
          if(approx.mean==1){ #6th Order Taylor
            
            rho.seq.Muu.inv<- rho^seq(0,12) #length 13
            rho.seq.Muu.inv.Mus<- rho^seq(1,14)  #length 14
            
            Muu.inv<-Reduce("+",get.list(Muu.inv.list,rho.seq.Muu.inv))  
            Muu.inv.Mus<-   Reduce("+",get.list(Muu.inv.Mus.list,rho.seq.Muu.inv.Mus))  
            
            C<- - Muu.inv.Mus
            
            BX<-Reduce("+",get.list(WX.list,rho.seq))
            Eyr<-BX[-s,]%*%beta +  C%*%(ys-BX[s,]%*%beta) 
            
          }#end
          
          if(approx.mean==2){ #LeSage Method
            
            C<-I- M/diag(M)
            BX<-Reduce("+",get.list(WX.list,rho.seq))
            Eyr<-BX[-s,]%*%beta0 +  C[-s,s]%*%(ys-BX[s,]%*%beta0)
            
            Muu.inv<-solve(M[-s,-s])
            
          }#end
          
          
          
        }else{
         
        #inside E-step
          #cat("beta",beta,"\n")
          #cat("rho",rho,"\n")
          #cat("omega",omega,"\n")
        
          
          A<-I-rho*W
          #Ainv<-solve(A)       
          #V<-Ainv%*%t(Ainv)
          #C1<- V[-s,s]%*%solve(V[s,s])  #Vus*Vss^{-1}
          # max(abs(C1-C[-s,s]))
          
          C2<- - solve(M[-s,-s],M[-s,s],sparse=TRUE)
  
          #sum(abs(C2-C1))
          
          
          AinvX<-solve(A,X)
          #print(dim(AinvX))
          #print(beta)
          #print(dim(C1))
          Eyr<-AinvX[-s,]%*%beta +  C2%*%(ys-AinvX[s,]%*%beta)
        }
        
        
        Ey<-rBind(ys,Eyr)
        rownames(Ey)<-1:n
        #Ey 506 x 1 Matrix of class "dgeMatrix"
        # Matrix(X%*%beta0)506 x 1 Matrix of class "dgeMatrix"
        
        XtEy<-t(X)%*%Ey
        WEy<-W%*%Ey
        XtWEy<-t(X)%*%WEy
                
        
        if(method==2 | method==3){
          M0rrinv<- solve(M[-s,-s])
          #M0rrinv<-sparseinv(M[-s,-s])
        
        }
        if(method==4){
         M0rrinv<-M[-s,-s]  #don't provide inv but M0rr
        }
        if(method==1 | method==5 | method==6){M0rrinv<-NULL}
      
        if(method==7){
          #6th order for trace
          rho.seq.Muu.inv<- rho^seq(0,12) #length 13
          M0rrinv<-Reduce("+",get.list(Muu.inv.list,rho.seq.Muu.inv))
        }  
        if(method==8){
          #first order for trace
          M0rrinv <- (I[-s,-s] + rho*WplusWt[-s,-s]-rho^2*WtW[-s,-s])   
        }
      
        
        res<- try(optimise(f=loglikSAM.EM,interval=c(-1+tol,1-tol),s,n,ns,XtEy,XtWEy,WEy,Ey,X,XtX,XtWX,XtWtWX,M0rrinv,omega0,chol.list,WplusWt,WtW,method=method))
        
        rho<-res$minimum
        Ly<-res$objective
        
        #print(res)
        #print(s)
        
        omega<-as.numeric(attr(res$objective,"omega"))
        beta<-as.matrix(attr(res$objective,"beta"))
        
        #theta.hat<-rbind(theta.hat,c(rho,omega,beta))
      
      if(iter%%10==0){
       
        cat("SAR Iteration: ",iter,"\n")  
         
      cat("beta",beta,"\n")
      cat("rho",rho,"\n")
      cat("omega",omega,"\n")
      cat("complete Log-Lik",Ly,"\n")
      
      
      #nromega<-nr*omega
      
  
      
      
      diff<-abs(rho-rho0)+abs(omega-omega0)+sum(abs(beta-beta0))
      cat("diff",diff,"\n")
        }
      #Sys.sleep(1)
      
      
      if(abs(omega)> (1/tol)^2 | abs(rho)> (1/tol)^2 | sum(abs(beta))> (1/tol)^2){
        conv<-FALSE;return(list(beta=beta,rho=rho,sigma2=omega,conv=conv)) 
      }#end if
      
    }#end while 
    
  } #end if  SAM
  #####################################################################
  
  
 
  ##################################################################################
  if(model=="SARerr"){
  
  
  #check whether symmetric
  # needed
  # ys,Eus,Vus,X,XtX,WX,XtWX,XtWtWX,WEY
  XtX<-t(X)%*%X
  
  WX.list<-get.W.lists(W,X,K=50,Xonly=TRUE)$WX.list
  
  WX<-WX.list[[1]]
  XtWX<-t(X)%*%WX
  XtWtWX<-t(WX)%*%WX
 
  
  XtWplusWt<-t(X)%*%WplusWt
  XtWtW<-t(X)%*%WtW
  
  ys<-Matrix(ys,ns,1,dimnames=list(1:ns,NULL))
  
  #these only depend on ys - don't need update
  

  
  rho<-rho0
  beta<-beta0
  omega<-omega0
 
  #maxrhodiff<-0.1
  #theta.hat<-matrix(c(rho0,omega0,beta0),1,2+p)
  
 
  Ls<-NULL
  
  diff<-1
  iter<-0
  while(diff>tol & iter<iter.max){
  
  rho0<-rho
  beta0<-beta
  omega0<-omega
  
    iter<-iter+1
  
  
  
  I<-.symDiagonal(n)
  M<-(I-rho*WplusWt+rho^2*WtW)
  
  # check: sum(abs(M/diag(M)-Diagonal(n,1/diag(M))%*%M))
  # update Yu with E(Yu|ys)
  
  if(approx.mean){
    
    if(approx.mean==1){ #5th Order Taylor
      
      rho.seq.Muu.inv<- rho^seq(0,12) #length 13
      rho.seq.Muu.inv.Mus<- rho^seq(1,14)  #length 14
      
      Muu.inv<-Reduce("+",get.list(Muu.inv.list,rho.seq.Muu.inv))  
      Muu.inv.Mus<-   Reduce("+",get.list(Muu.inv.Mus.list,rho.seq.Muu.inv.Mus))  
      
      C<- - Muu.inv.Mus
      
      Eyr<- Xr%*%beta+  C%*%(ys-Xs%*%beta)
      
      
      
    }#end
    
    if(approx.mean==2){ #LeSage Method
      
      
      C<-I- M/diag(M)
      Eyr<- Xr%*%beta+  C[-s,s]%*%(ys-Xs%*%beta)
      
      Muu.inv<-solve(M[-s,-s])
      
    }#end
    
  }else{
    #exact method
    #Muu.inv<-solve(M[-s,-s])
    
    C2<- - solve(M[-s,-s],M[-s,s],sparse=TRUE)
    #sum(abs(C2-C1))
    
    
    Eyr<-X[-s,]%*%beta +  C2%*%(ys-X[s,]%*%beta)
  }
  
  Ey<-rBind(ys,Eyr)
  
  r<-Ey-X%*%beta
  rtr<-t(r)%*%r
  Wr<-W%*%r
  rtWr<-t(r)%*%Wr
  rtWtWr<- t(Wr)%*%Wr
  
  XtEy<-t(X)%*%Ey
  XtWplusWtEy<- XtWplusWt%*%Ey
  XtWtWEy<-XtWtW%*%Ey
  
  
  
  #M0rr<-M[-s,-s]
  if(method==2 | method==3){
    M0rrinv<- solve(M[-s,-s])
    #M0rrinv<-sparseinv(M[-s,-s])
    
  }
  if(method==4){
    M0rrinv<-M[-s,-s]  #don't provide inv but M0rr
  }
  if(method==1 | method==5 | method==6){M0rrinv<-NULL}
  
  if(method==7){
    #6th order for trace
    rho.seq.Muu.inv<- rho^seq(0,12) #length 13
    M0rrinv<-Reduce("+",get.list(Muu.inv.list,rho.seq.Muu.inv))
  }  
  if(method==8){
    #first order for trace
    M0rrinv <- (I[-s,-s] + rho*WplusWt[-s,-s]-rho^2*WtW[-s,-s])  
    
  }
  
  #Mrrinv<-sparseinv(M[-s,-s])
  #cat("method",method,"\n")
  
  res<- try(optimise(f=loglikSEM.EM,interval=c(-1+tol,1-tol),s,n,ns,rtr,rtWr,rtWtWr,Ey,XtEy,XtWplusWtEy,XtWtWEy,omega0,M0rrinv,XtX,XtWX,XtWtWX,X,WX,W,chol.list,WplusWt,WtW,method=method))
  
  rho<-res$minimum
  Ly<-res$objective
  #cat("optimisation done\n")
  #record list of objective functions
  #hilf<-AR.exact.log.lik.fct(rho,ys,W,X,WplusWt,WtW,ARmod=2,chol.list,pr=TRUE,approx=F,REML=F)
  #cat("Ls done\n")
  # record list of true objective fct.
  #Ls<-c(Ls,hilf)
  
  
  
  omega<-as.numeric(attr(res$objective,"omega"))
  beta<-as.matrix(attr(res$objective,"beta"))
  
  #theta.hat<-rbind(theta.hat,c(rho,omega,beta))
  
  if(iter%%10==0){
  
    cat("SARerr Iteration: ",iter,"\n")  
    
  cat("beta",beta,"\n")
  cat("rho",rho,"\n")
  cat("omega",omega,"\n")
  cat("complete Log-Lik",Ly,"\n")
  
  
  #nromega<-nr*omega
  
  diff<-abs(rho-rho0)+abs(omega-omega0)+sum(abs(beta-beta0))
  cat("diff",diff,"\n")
  #Sys.sleep(1)
  
  }#end 
  
  if(abs(omega)> (1/tol)^2 | abs(rho)> (1/tol)^2 | sum(abs(beta))> (1/tol)^2){
    conv<-FALSE;return(list(beta=beta,rho=rho,sigma2=omega,conv=conv)) 
  }#end if
  
  }#end while 
  
  } #end if  SEM
  #####################################################################
  
  
  
  if(model=="CAR"){
    
  } #end if  

  if(diff>tol | iter>=iter.max){conv<-FALSE}else{conv<-TRUE}
  
  
 ####################### END EM now STANDARD ERRORS ######################## 
  #cat("list of approx. Ls's\n") 
  #cat(Ls,"\n")
  
  
  #calculate INFO matrix
  if(se==TRUE){
    
    if(model=="SAR"){
      
      A<-I - rho*W
      Ainv<- solve(A)
      M<-(I-rho*WplusWt+rho^2*WtW)
      dMdrho<- - WplusWt+2*rho*WtW
      
      
      
      AinvX<-Ainv%*%X
      BsX<-AinvX[s,]
     
      Muu<-M[-s,-s]
      Mss<-M[s,s]
      Msu<-M[s,-s]
      Mus<-t(Msu)
      Muu.inv<-solve(Muu)
      Bss<- Mss - Msu%*%Muu.inv%*%t(Msu)
      rs<-ys - BsX%*%beta
      rstBss<-t(rs)%*%Bss
      rstBssrs<-rstBss%*%rs
      BsXBssrs<-t(BsX)%*%t(rstBss)
      AinvW<-Ainv%*%W
      AinvWAinvX<-AinvW%*%AinvX
      dBsXdrho<- AinvWAinvX[s,]
      
      drsdrho<- - dBsXdrho%*%beta
      
        
      dMsudrho.Muu.inv.Mus<- dMdrho[s,-s]%*%Muu.inv%*%Mus
      Muu.inv.Muu.Muu.inv<-Muu%*%Muu.inv%*%Muu.inv
      
      
      dBssdrho<- dMdrho[s,s] - (dMsudrho.Muu.inv.Mus + t(dMsudrho.Muu.inv.Mus)) + Msu%*%Muu.inv.Muu.Muu.inv%*%Mus
        
      BsX.dBssdrho.rs <- t(BsX)%*%dBssdrho%*%rs
      BsX.Bss.drsdrho <-  t(BsX)%*%Bss%*%drsdrho 
      dBsXdrho.Bss.rs <- t(dBsXdrho)%*%t(rstBss)
      
      Bssinv<-solve(Bss)
      Bssinv.dBsdrho<-Bssinv%*%dBssdrho
      
      d2rsdrho2<- 2*(AinvW%*%AinvWAinvX%*%beta)[s,]
      
      d2Mdrho2<- 2*WtW
      
      d2Msudrho2.Muu.inv.Mus<-d2Mdrho2[s,-s]%*%Muu.inv%*%Mus
      dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus<-dMdrho[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      
      d2Bssdrho2<- d2Mdrho2[s,s] - (d2Msudrho2.Muu.inv.Mus + t(d2Msudrho2.Muu.inv.Mus))
      d2Bssdrho2<-d2Bssdrho2 + 2*(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus+t(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus))
      d2Bssdrho2<-d2Bssdrho2 - 2*dMdrho[s,-s]%*%Muu.inv%*% dMdrho[-s,s] 
      d2Bssdrho2<-d2Bssdrho2 - 2*M[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      d2Bssdrho2<-d2Bssdrho2 + Msu%*%Muu.inv%*%d2Mdrho2[-s,-s]%*%Muu.inv%*%Mus
      
      I.rho.rho <-    0.5*tr.prod(Bssinv.dBsdrho,Bssinv.dBsdrho) - 0.5*tr.prod(Bssinv,d2Bssdrho2)
      I.rho.rho <- I.rho.rho + 1/omega*( t(drsdrho)%*%Bss%*%drsdrho + t(d2rsdrho2)%*%t(rstBss) +2*t(drsdrho)%*%dBssdrho%*%rs +t(rs)%*%d2Bssdrho2%*%rs)
      
      I.rho.omega <- -1/(2*omega^2)*(2*rstBss%*%drsdrho+t(rs)%*%dBssdrho%*%rs)
      I.rho.beta  <- -(BsX.dBssdrho.rs+BsX.Bss.drsdrho+dBsXdrho.Bss.rs)/omega
      I.omega.omega <- -ns/(2*omega^2) +rstBssrs/omega^3
      I.beta.beta <- t(BsX)%*%Bss%*%BsX/omega 
      I.omega.beta <- - BsXBssrs/omega^2
      
      
    }#end if SAR
    
    if(model=="SARerr"){
      
      
      M<-(I-rho*WplusWt+rho^2*WtW)
      dMdrho<- - WplusWt+2*rho*WtW
      
    
      Muu<-M[-s,-s]
      Mss<-M[s,s]
      Msu<-M[s,-s]
      Mus<-t(Msu)
      Muu.inv<-solve(Muu)
      Bss<- Mss - Msu%*%Muu.inv%*%t(Msu)
      rs<-ys - Xs%*%beta
      rstBss<-t(rs)%*%Bss
      rstBssrs<-rstBss%*%rs
      XsBssrs<-t(Xs)%*%t(rstBss)
      
      
      
      dMsudrho.Muu.inv.Mus<- dMdrho[s,-s]%*%Muu.inv%*%Mus
      Muu.inv.Muu.Muu.inv<-Muu%*%Muu.inv%*%Muu.inv
      
      
      dBssdrho<- dMdrho[s,s] - (dMsudrho.Muu.inv.Mus + t(dMsudrho.Muu.inv.Mus)) + Msu%*%Muu.inv.Muu.Muu.inv%*%Mus
      
      Xs.dBssdrho.rs <- t(Xs)%*%dBssdrho%*%rs
      
      
      Bssinv<-solve(Bss)
      Bssinv.dBsdrho<-Bssinv%*%dBssdrho
      
      
      d2Mdrho2<- 2*WtW
      
      d2Msudrho2.Muu.inv.Mus<-d2Mdrho2[s,-s]%*%Muu.inv%*%Mus
      dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus<-dMdrho[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      
      d2Bssdrho2<- d2Mdrho2[s,s] - (d2Msudrho2.Muu.inv.Mus + t(d2Msudrho2.Muu.inv.Mus))
      d2Bssdrho2<-d2Bssdrho2 + 2*(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus+t(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus))
      d2Bssdrho2<-d2Bssdrho2 - 2*dMdrho[s,-s]%*%Muu.inv%*% dMdrho[-s,s] 
      d2Bssdrho2<-d2Bssdrho2 - 2*M[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      d2Bssdrho2<-d2Bssdrho2 + Msu%*%Muu.inv%*%d2Mdrho2[-s,-s]%*%Muu.inv%*%Mus
      
      I.rho.rho <-    0.5*tr.prod(Bssinv.dBsdrho,Bssinv.dBsdrho) - 0.5*tr.prod(Bssinv,d2Bssdrho2)
      I.rho.rho <- I.rho.rho + t(rs)%*%d2Bssdrho2%*%rs/omega
      
      I.rho.omega <- -1/(2*omega^2)*(t(rs)%*%dBssdrho%*%rs)
      I.rho.beta  <- -(Xs.dBssdrho.rs)/omega
      I.omega.omega <- -ns/(2*omega^2) +rstBssrs/omega^3
      I.beta.beta <- t(Xs)%*%Bss%*%Xs/omega 
      I.omega.beta <- - XsBssrs/omega^2 
      
    }#end if SARerr
    
    if(model=="CAR"){
      
      M<-(I-rho*W)
      dMdrho<- - W
      
      
      Muu<-M[-s,-s]
      Mss<-M[s,s]
      Msu<-M[s,-s]
      Mus<-t(Msu)
      Muu.inv<-solve(Muu)
      Bss<- Mss - Msu%*%Muu.inv%*%t(Msu)
      rs<-ys - X[s,]%*%beta
      rstBss<-t(rs)%*%Bss
      rstBssrs<-rstBss%*%rs
      XsBssrs<-t(X[s,])%*%t(rstBss)
      
      
      
      dMsudrho.Muu.inv.Mus<- dMdrho[s,-s]%*%Muu.inv%*%Mus
      Muu.inv.Muu.Muu.inv<-Muu%*%Muu.inv%*%Muu.inv
      
      
      dBssdrho<- dMdrho[s,s] - (dMsudrho.Muu.inv.Mus + t(dMsudrho.Muu.inv.Mus)) + Msu%*%Muu.inv.Muu.Muu.inv%*%Mus
      
      Xs.dBssdrho.rs <- t(BsX)%*%dBssdrho%*%rs
      
      
      Bssinv<-solve(Bss)
      Bssinv.dBsdrho<-Bssinv%*%dBssdrho
      
      
      
      
      
      d2Msudrho2.Muu.inv.Mus<-d2Mdrho2[s,-s]%*%Muu.inv%*%Mus
      dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus<-dMdrho[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      
      
      d2Bssdrho2<-  2*(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus+t(dMsudrho.Muu.inv.dMuu.drho.Muu.inv.Mus))
      d2Bssdrho2<-d2Bssdrho2 - 2*dMdrho[s,-s]%*%Muu.inv%*% dMdrho[-s,s] 
      d2Bssdrho2<-d2Bssdrho2 - 2*M[s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%dMdrho[-s,-s]%*%Muu.inv%*%Mus
      
      
      I.rho.rho <-    0.5*tr.prod(Bssinv.dBsdrho,Bssinv.dBsdrho) - 0.5*tr.prod(Bssinv,d2Bssdrho2)
      I.rho.rho <-     I.rho.rho + t(rs)%*%d2Bssdrho2%*%rs/omega
      
      I.rho.omega <- -1/(2*omega^2)*(t(rs)%*%dBssdrho%*%rs)
      I.rho.beta  <- -(Xs.dBssdrho.rs)/omega
      I.omega.omega <- -ns/(2*omega^2) +rstBssrs/omega^3
      I.beta.beta <- t(Xs)%*%Bss%*%Xs/omega 
      I.omega.beta <- - XsBssrs/omega^2 
      
    }#end if CAR
    
    #print(I.rho.rho)
    #print(I.rho.omega)
    #print(I.rho.beta)
    #print(I.beta.beta)
    #print(I.omega.omega)
    #print(I.omega.beta)
    
    
    l.theta<-p+2
      
    Info<-matrix(NA,l.theta,l.theta) 
    Info[1,1]<-as.numeric(I.rho.rho)
    Info[1,2]<-Info[2,1]<-as.numeric(I.rho.omega)
    Info[2,2]<-as.numeric(I.omega.omega)
    Info[3:l.theta,1]<-as.numeric(I.rho.beta)
    Info[1,3:l.theta]<-as.numeric(I.rho.beta)
    Info[3:l.theta,2]<-as.numeric(I.omega.beta)
    Info[2,3:l.theta]<-as.numeric(I.omega.beta)
    Info[3:l.theta,3:l.theta]<-as.matrix(I.beta.beta)
    
  
  Cov<-as.matrix(solve(Info))
  #print(Cov)
  
  rho.se<-sqrt(Cov[1,1])
  omega.se<-sqrt(Cov[2,2])
  beta.se<-sqrt(diag(Cov[3:l.theta,3:l.theta]))
  
  
  #rho<-theta.hat[L+1,1]
  #omega<-theta.hat[L+1,2]
  #beta<-theta.hat[L+1,3:l.theta]
  
  return(list(beta=beta,rho=rho,sigma2=omega,conv=conv,Cov=Cov,rho.se=rho.se,
              beta.se=beta.se,sigma2.se=omega.se,L=Ly,iter=iter,Ls=Ls))  
  
  }#end if(se=T){
    
    
  
  #if successful
  return(list(beta=beta,rho=rho,sigma2=omega,conv=conv,L=Ly,iter=iter,Ls=Ls))  
  
  
}#end EM
##################################################################################




#################################################################################
loglikSEM.EM<-function(rho,s,n,ns,rtr,rtWr,rtWtWr,Ey,XtEy,XtWplusWtEy,XtWtWEy,omega0,M0rrinv,XtX,XtWX,XtWtWX,X,WX,W,chol.list,WplusWt=W+t(W),WtW=t(W)%*%W,method=1){
#Y contains filled data
  
# M  = I -rho*(W+Wt)+rho^2*WtW 
  #A  <- -rho*W
  #diag(A)<-1
  #print(rho)
  #print(XtX)
  #print(XtWX)
  #print(XtWtWX)
  
  XtMX<- XtX - rho*(XtWX+t(XtWX)) + rho^2*XtWtWX  
  
  
  XtMEy <- XtEy - rho*XtWplusWtEy + rho^2*XtWtWEy
  #print(XtMEy)
   #For SEM and CAR
  beta<-solve(XtMX)%*%XtMEy
   
  #print(beta)
  
  rtMr<- rtr - 2*rho*rtWr + rho^2*rtWtWr 
  
  T5<-0
  n1<-n
  
  switch( method,
{
  #  default method=1: assume at end of M-step that  rho=rho0 and omega=omega0
  omega<- rtMr/ns
  T4<-n
}
,
  { 
  #  method=2: correct method - i.e. calculate trace as rho/rho0 and omega/omega0 differ
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
  }
  ,
  {
    #REML version
   p<-length(beta)  
   n1<-n-p
   #method=3: if we use   omega<- rtMr/ns, then last term does not equal "n" unless rho=rho0 
   Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
    
   trace.term<-  tr.prod(Mrr,M0rrinv)
   
   omega<- (rtMr+omega0*trace.term)/n1
   
   T4<-(rtMr+omega0*trace.term)/omega
   
   T5<-log(det(XtMX))
  
   #print(T5)
   #print(T4)
   
  }
,
{
  # method=4: use approximation of trace  
  #first order approx of trace
  Mrr <- .symDiagonal(n-ns) - rho*WplusWt[-s,-s] + rho^2*WtW[-s,-s]
  trace.term<- sum(diag(Mrr)/diag(M0rrinv)) #just consider diagonals
  #print(trace.term)
  omega<- (rtMr+omega0*trace.term)/n    
  T4<-n
}  ,
  {#method=5: LeSage and Pace 2004
    rs<-Ey-X%*%beta
    rs<-(rs-rho*W%*%rs)[s,]
    omega<-t(rs)%*%rs/ns
      T4<-n
  }
,

{ #method=6: assuming rho=rho'
  omega<- (rtMr+omega0*(n-ns))/n    
  T4<-n  
  
}
,
{ 
  #  method=7: correct except that M0rr is 6th order approx
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
}
,
{ 
  #  method=8:correct except that M0rr is 1st order approx
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
}


  )
  
  if(chol.list$sim){
  #  symmetric
   logdetArho<-ldetA(rho,chol.list) 
       
  }else{
    #not symmetric
    M<- .symDiagonal(n) -rho*WplusWt + rho^2*WtW;
    C<-Cholesky(-rho*WplusWt+rho^2*WtW,super=TRUE,Imult=1);logdetArho<- determinant(C)$modulus
  }#end 

  #cat("calculated\n")
  #print(2*logdetArho)
  #cat("true\n")
  #print(log(det(M)))
  #print(rho)
  #   terms                    T1               T2           T3               T4
  
  # last term n is (rtM + romega0*trace.term)/omega=n, as omega= (rtM + romega0*trace.term)/n
  
  #det.log.M is missing
  L<- as.numeric(0.5*(n*log(2*pi) - 2*logdetArho +   n1*log(omega) + T4   + T5))            
  attr(L, "omega") <-   omega  
  attr(L, "beta") <-   beta  
  return(L)
}#end function
#############################
#  END  SEM
#############################

#############################
#  BEGIN  SAM
#############################
loglikSAM.EM<-function(rho,s,n,ns,XtEy,XtWEy,WEy,Ey,X,XtX,XtWX,XtWtWX,M0rrinv,omega0,chol.list,WplusWt=W+t(W),WtW=t(W)%*%W,method=1){
  #Y contains filled data
  
  # M  = I -rho*(W+Wt)+rho^2*WtW 
  #A  <- -rho*W
  #diag(A)<-1
  #print(rho)
  #print(XtX)
  #print(XtWX)
  #print(XtWtWX)

  
  
  XtAEy <- XtEy - rho*XtWEy
  #print(XtMEy)
  #For SEM and CAR
  beta<-solve(XtX)%*%XtAEy
  
  #print(beta)
  Xbeta<-X%*%beta
  # r_A = A*Y-Xbeta= (I-rhoW)Y-Xbeta=Y-rho*W*Y -Xbeta 
  r<- Ey -rho*WEy -Xbeta
  
  rtMr<- t(r)%*%r
  
  T5<-0
  n1<-n
  
  switch( method,
{
  #  default method=1: assume at end of M-step that  rho=rho0 and omega=omega0
  omega<- rtMr/ns
  T4<-n
}
,
{ 
  #  method=2: correct method - i.e. calculate trace as rho/rho0 and omega/omega0 differ
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
}
,
{
  #ethod=3: REML version
  XtMX<- XtX  -  rho*(XtWX+t(XtWX))+rho^2*XtWtWX
  
  p<-length(beta)  
  n1<-n-p
  #method=3: if we use   omega<- rtMr/ns, then last term does not equal "n" unless rho=rho0 
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<-  tr.prod(M0rrinv,Mrr)
  
  omega<- (rtMr+omega0*trace.term)/n1
  
  T4<-(rtMr+omega0*trace.term)/omega
  
  T5<-log(det(XtMX))
  
  #print(T5)
  #print(T4)
  
}
,
{
# method=4: use approximation of trace  
#first order approx of trace
Mrr <- .symDiagonal(n-ns) - rho*WplusWt[-s,-s] + rho^2*WtW[-s,-s]
trace.term<-sum(diag(Mrr)/diag(M0rrinv)) #just consider diagonals
#print(trace.term)
omega<- (rtMr+omega0*trace.term)/n    
T4<-n
}
,
{#method=5: LeSage and Pace 2004
  omega<-t(r[s,])%*%r[s,]/ns
  T4<-n
},


{ #method=6: assuming rho=rho'
  omega<- (rtMr+omega0*(n-ns))/n    
  T4<-n  
  
}
,
{ 
  #  method=7: correct except that M0rr is 6th order approx
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
}
,
{ 
  #  method=8:correct except that M0rr is 1st order approx
  Mrr<- .symDiagonal(n-ns) - rho*WplusWt[-s,-s]+rho^2*WtW[-s,-s]
  
  trace.term<- tr.prod(M0rrinv,Mrr) # tr.prod(Mrr,solve(M0rr))
  
  #print(trace.term)
  
  omega<- (rtMr+omega0*trace.term)/n
  
  T4<-n
}
  )

if(chol.list$sim){
  #  symmetric
  logdetArho<-ldetA(rho,chol.list) 
  
}else{
  #not symmetric
  M<- .symDiagonal(n) -rho*WplusWt + rho^2*WtW;
  C<-Cholesky(-rho*WplusWt+rho^2*WtW,super=TRUE,Imult=1);logdetArho<- determinant(C)$modulus
}#end 

#cat("calculated\n")
#print(2*logdetArho)
#cat("true\n")
#print(log(det(M)))
#print(rho)
#   terms                    T1               T2           T3               T4

# last term n is (rtM + romega0*trace.term)/omega=n, as omega= (rtM + romega0*trace.term)/n

#det.log.M is missing
L<- as.numeric(0.5*(n*log(2*pi) - 2*logdetArho +   n1*log(omega) + T4   + T5))            
attr(L, "omega") <-   omega  
attr(L, "beta") <-   beta  
return(L)
}#end function
# END SAM /SAR
#################################



calculate.det.sparse.symm<-function(W,eps=1e-4,lower=0,upper=1){
#calculates the determinant of A=I-rho*W for a vector of rho-values
# rho from lower to upper with step-size eps
# and from -upper to - lower with step-size
  # by default for all rho from -1 to +1
  
  # check symmetry
  if(sum(abs((USCounties+t(USCounties))/2-USCounties))>eps){
    cat("not symmetric\n")
    return(NULL)
  }
  IM <- .symDiagonal(n)
  nWC <- -W
  pWC <-  W
  
  rho<-seq(lower+eps,upper-eps,eps)
  #rho.n<- -rho.p
  #rho<-rho.p
  
  #system.time(MJp <- sapply(rho.p, function(x)
  #  determinant(IM - x * USCounties, logarithm = TRUE)$modulus))
  #system.time(MJn <- sapply(rho.n, function(x)
  #  determinant(IM - x * USCounties, logarithm = TRUE)$modulus))
  
  Cn <- Cholesky(nWC, super = TRUE, Imult = 2) #needed for positive rho's
  Cp <- Cholesky(pWC, super = TRUE, Imult = 2) #needed for negative rho's 
  MJprho <- n * log(rho) + Matrix:::ldetL2up(Cn, nWC, 1/rho)
  MJnrho <- n * log(rho) + Matrix:::ldetL2up(Cp, pWC, 1/rho) 
  
  return(rho=c(-rho,0,rho),log.det.A=c(MJnrho,1,MJnrho))
  
}#end calculate.det.sparse
  
calculate.det.sparse.unsymm<-function(WplusWt,WtW,eps=1e-4){
  
  
  IM <- .symDiagonal(n)
  nWC <- -W
  pWC <-  W
  
  rho<-seq(eps,1-eps,eps)
  #rho.n<- -rho.p
  #rho<-rho.p
  
  #system.time(MJp <- sapply(rho.p, function(x)
  #  determinant(IM - x * USCounties, logarithm = TRUE)$modulus))
  #system.time(MJn <- sapply(rho.n, function(x)
  #  determinant(IM - x * USCounties, logarithm = TRUE)$modulus))
  
  Cn <- Cholesky(nWC, super = TRUE, Imult = 2) #needed for positive rho's
  Cp <- Cholesky(pWC, super = TRUE, Imult = 2) #needed for negative rho's 
  MJprho <- n * log(rho) + Matrix:::ldetL2up(Cn, nWC, 1/rho)
  MJnrho <- n * log(rho) + Matrix:::ldetL2up(Cp, pWC, 1/rho) 
  
  return(rho=c(-rho[length(rho):1],0,rho),log.det.A=c(MJnrho[length(rho):1],0,MJnrho))
  
}#end calculate.det.sparse


###############################################################
#maximising marginal log-likelihood with approximation
###################################################################
AR.pop.approx.sparse<-function(Ys,W,X,model="SAR",se=TRUE,tol=1e-6,K=10,Wlists=NULL,s=1:length(Ys),rho.limits=c(-1,1)){
  
  #ARmod="SAR"
  #ARmod="SARerr"
  #ARmod="CAR"  
  if(model=="SAR"){ARmod<-1}
  if(model=="SARerr"){ARmod<-2}
  if(model=="CAR"){ARmod<-3}
  
  #Obtain W.lists have to have 3 elements: WX.list,WWt.list,WplusWt.list
  if(is.null(Wlists)){
  cat("Calculating W-lists started! \n")
  Wlists<-get.W.lists(W,X,K=K)
  cat("Calculating W-lists finished! \n")
  }
  p<-dim(X)[2]
  p1<-p+1
  
  
  ns<-length(Ys)
  
  WX.list.sample <- make.sample.list(Wlists$WX.list,s,1:p)
  WWt.list.sample <- make.sample.list(Wlists$WWt.list,s,s)
  WplusWt.list.sample <-  make.sample.list(Wlists$WplusWt.list,s,s)
 
  
  
  n<-dim(X)[1]
  
  Xs<-as.matrix(X[s,])
  
  
  #check
  
  
  
  eps<-tol
  
  counter<<-0
  
  
system.time({  
  switch(ARmod,
{
  #SAM=SAR
  res<- try(optimize(AR.pop.approx,interval=c(rho.limits[1]+eps,rho.limits[2]-eps),WX.list.sample,Ys,Xs=Xs,ARmod=1,grad=FALSE,Info=FALSE,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=K,tol =tol^2))
  
}
,
{
  #SEM=SARerr
  res<-try(optimize(AR.pop.approx,interval=c(rho.limits[1]+eps,rho.limits[2]-eps),WX.list.sample=NULL,Ys,Xs=Xs,ARmod=2,grad=FALSE,Info=FALSE,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=K,tol=tol^2))
}
,
{ 
  #CAR
  res<-try(optimize(AR.pop.approx,interval=c(rho.limits[1]+eps,rho.limits[2]-eps),WX.list.sample=NULL,Ys,Xs=Xs,ARmod=3,grad=FALSE,Info=FALSE,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=K,tol=tol^2))
  
})

}) #end system.time
cat(counter," Iterations needed.\n")
CleanEnvir(counter)

rho<-res$minimum




L<-AR.pop.approx(x=rho,WX.list.sample,Ys,Xs=Xs,ARmod=ARmod,grad=se,Info=se,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=K)

# print("hello")

omega<-attr(L, "omega")
beta<-attr(L, "beta") 
if(se){
  grad<-attr(L,"gradient")
  Info<-attr(L,"Info")
  
  # return(list(grad=grad,Cov=Cov,Cov.sigma2=Cov.omega,Cov.beta=Cov.beta,Cov.rho=Cov.rho,Info=Info))
  
  
  Cov<-try(solve(Info))
  if(!inherits(Cov,"try-error")){
  Cov.omega<-Cov[1,1]
  Cov.beta<-Cov[2:p1,2:p1]
  Cov.rho<-Cov[p+2,p+2]
  if(p>1){beta.se<-sqrt(diag(Cov.beta))}else{beta.se<-sqrt(Cov.beta)}
  rho.se<-sqrt(Cov.rho)
  sigma2.se<-sqrt(Cov.omega)
  }else{
    sigma2.se<-rho.se<-beta.se<-Cov.rho<-Cov.beta<-Cov.omega<-NULL     
  }
}
L<- c(-L)




return(list(L=L,beta=beta,rho=rho,sigma2=omega,grad=grad,Info=Info,Cov=Cov,Cov.sigma2=Cov.omega,Cov.rho=Cov.rho,Cov.beta=Cov.beta,
            beta.se=beta.se,rho.se=rho.se,sigma2.se=sigma2.se,Wlists=Wlists))

}#end function

AR.approx.Info<-function(Ys,W,X,rho,model="SAR",tol=1e-6,K=10,Wlists=NULL,s=1:length(Ys),rho.limits=c(-1,1)){
  
  
  #ARmod="SAR"
  #ARmod="SARerr"
  #ARmod="CAR"  
  if(model=="SAR"){ARmod<-1}
  if(model=="SARerr"){ARmod<-2}
  if(model=="CAR"){ARmod<-3}
  
  #Obtain W.lists have to have 3 elements: WX.list,WWt.list,WplusWt.list
  if(is.null(Wlists)){
    cat("Calculating W-lists started! \n")
    Wlists<-get.W.lists(W,X,K=K)
    cat("Calculating W-lists finished! \n")
  }
  p<-dim(X)[2]
  p1<-p+1
  
  ns<-length(Ys)
  
  WX.list.sample <- make.sample.list(Wlists$WX.list,s,1:p)
  WWt.list.sample <- make.sample.list(Wlists$WWt.list,s,s)
  WplusWt.list.sample <-  make.sample.list(Wlists$WplusWt.list,s,s)
  
  
  
  n<-dim(X)[1]
  
  Xs<-as.matrix(X[s,])
  
  
  #check
  
  
L<-AR.pop.approx(x=rho,WX.list.sample,Ys,Xs=Xs,ARmod=ARmod,grad=TRUE,Info=TRUE,pr=TRUE,WWt.list.sample,WplusWt.list.sample,K=K)


omega<-attr(L, "omega")
beta<-attr(L, "beta") 


  grad<-attr(L,"gradient")
  Info<-attr(L,"Info")
  
  # return(list(grad=grad,Cov=Cov,Cov.sigma2=Cov.omega,Cov.beta=Cov.beta,Cov.rho=Cov.rho,Info=Info))
  
  
  Cov<-try(solve(Info))
  if(!inherits(Cov,"try-error")){
    Cov.omega<-Cov[1,1]
    Cov.beta<-Cov[2:p1,2:p1]
    Cov.rho<-Cov[p+2,p+2]
    beta.se<-sqrt(diag(Cov.beta))
    rho.se<-sqrt(Cov.rho)
    sigma2.se<-sqrt(Cov.omega)
  }else{
    Cov.omega<- Cov.beta<-Cov.rho<-beta.se<-rho.se<-sigma2.se<-NULL
    
  } 

L<- c(-L)


return(list(L=L,beta=beta,rho=rho,sigma2=omega,grad=grad,Info=Info,Cov=Cov,Cov.sigma2=Cov.omega,Cov.rho=Cov.rho,Cov.beta=Cov.beta,
            beta.se=beta.se,rho.se=rho.se,sigma2.se=sigma2.se,Wlists=Wlists))

  
}#end AR.approx.Info

##################################################################################
##################################################################################
#                        measurement error model
##################################################################################
##################################################################################
AR.error<-function(Y,W,X,model="SAR",se=TRUE,REML=FALSE,tol=1e-4,rho0=NULL,omega.eps0=NULL,omega.Y0=NULL,rho.limits=c(-1,1)){
  
  
  if(model=="SAR"){ARmod<-1}
  if(model=="SARerr"){ARmod<-2}
  if(model=="CAR"){ARmod<-3}
  
 
  n<-dim(W)[1]
  
  p<-dim(X)[2]
  
  
  rownames(W)<-colnames(W)<-1:n
  rownames(X)<-1:n
  W<- as(W, "CsparseMatrix")
  WplusWt<-forceSymmetric(W+t(W))
  WWt<-forceSymmetric(W%*%t(W))
  WY<-W%*%Y  
  WX<-W%*%X
  I<-Diagonal(n,rep(1,n))
  
  
  if(!is.null(omega.Y0) & !is.null(omega.eps0)){
  theta0<-omega.Y0/omega.eps0
  }else{
   theta0<-1  
  }
  if(!is.null(rho0)){
    rho0<-0
  }
  
  
  par<-c(rho0,log(theta0))
  
  # source("ARmodels.sparse.r")
  #AR.error.log.lik.fct(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=TRUE,REML=T)
  
  # meth<- "Nelder-Mead" 
  if(0){
  cat("Simulated Annealing\n")  
  st0<-system.time(
    res0 <- try(optim(par, fn=AR.error.log.lik.fct, gr = NULL,Y=Y,W=W,X=X,WplusWt=WplusWt,WWt=WWt,WY=WY,WX=WX,ARmod=ARmod,pr=TRUE,REML=REML,
                     method = "SANN",
                     lower = c(-1+tol,-10), upper = c(1-tol,10),
                     control = list(), hessian = FALSE))
    #control = list(factr=tol^2,pgtol=tol), hessian = FALSE))
    
    )
  par<-res0$par
  print(par)
  print(st0)
  
  
  }

  meth<-"L-BFGS-B"
  counter<<-0 
  cat("Optimisation L-BFGS-B\n")  
  st<-system.time(
  res <- try(optim(par, fn=AR.error.log.lik.fct, gr = NULL,Y=Y,W=W,X=X,WplusWt=WplusWt,WWt=WWt,WY=WY,WX=WX,ARmod=ARmod,pr=TRUE,REML=REML,
        method = meth,
        lower = c(-1+tol,-10), upper = c(1-tol,10),
        control = list(), hessian = FALSE))
        #control = list(factr=tol^2,pgtol=tol), hessian = FALSE))
  )
  #print(counter)
  L<-AR.error.log.lik.fct(res$par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=TRUE,REML=REML)
  cat("Number Iterations:",counter-1,"\n")
  cat("Time for fitting",st,"\n")
  
  beta<-as.matrix(attr(L,"beta"))
  rho<-attr(L,"rho")
  omega.Y<- as.numeric(attr(L,"omega_Y"))
  omega.eps<- as.numeric(attr(L,"omega_eps"))
  Lpos<- -c(L)
  
  log.theta<-res$par[2]
  
  
  cat("Fitting finished\n")
  
  
  cat("rho:",rho,"\n")
  cat("omega.Y:",omega.Y,"\n")
  cat("omega.eps:",omega.eps,"\n")
  
  
  
  if(se){
    
    cat("Calculation of Standard Errors started\n")
    #calculate first derivative
    #print(log.lik.fct(par1,W,Ys,X,WplusWt,WtW,ARmod=ARmod,pr=F))
    
    g<-try({
    f.rho<-function(rho,omega.Y,omega.eps,beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML){
    par<-c(rho,log(omega.Y)-log(omega.eps),omega.eps,beta)
      h<-AR.error.log.lik.fct(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=FALSE,REML=REML)
    #return(-as.matrix(attr(h,"Ctb")-attr(h,"CtC")%*%beta)/omega.eps)
    return(-(as.matrix(attr(h,"Ctb")-attr(h,"CtC")%*%beta))/omega.eps)
    }
    # Ctb-CtC%*%beta=(AX)t*K.theta.inv*AY- t(AX)%*%K.theta.inv%*%AX%*%beta = XtMr, derivative is -(1/omega_eps)*Xt%*% d(M*r)/drho  
    
    
    #correct
    f.rho.omegas<-function(par,beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML){
      rho<-par[1]
      omega.eps<-par[2]
      omega.Y<-par[3]
      par<-c(rho,log(omega.Y)-log(omega.eps),omega.eps,beta)
      h<-AR.error.log.lik.fct(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=FALSE,REML=REML)
      return(h)
    }
    f.all<-function(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML){
      p<-dim(X)[2]
      beta<-par[1:p]
      rho<-par[p+1]
      omega.Y<-par[p+3]
      omega.eps<-par[p+2]
      par<-c(rho,log(omega.Y)-log(omega.eps),omega.eps,beta)
      h<-AR.error.log.lik.fct(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=FALSE,REML=REML)
      return(h)
    }
    
    
    
    par.names<-c(paste("beta",1:p),"rho","omega.eps","omega.Y")
    
    #calculate 2nd derivative
    
    if(REML){n1<-n-p}else{n1<-n}
    
    Info.beta.beta<- as.matrix(attr(L,"CtC"))/omega.eps
    Info.beta.omega.Y   <- matrix(0,p,1)
    Info.beta.omega.eps <- matrix(0,p,1)
    if(ARmod==1){
      # f.rho(rho,omega.Y,omega.eps,beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML)
      # Info.beta.rho<- jacobian(f.rho,rho,method="Richardson",method.args=list(),omega.Y,omega.eps,beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML) 
      Info.beta.rho<-  f.rho(rho,omega.Y,omega.eps,beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML)
    }else{
      Info.beta.rho<-  matrix(0,p,1)
    }  
    
    
    Info.rho.omegas<-hessian(f.rho.omegas,c(rho,omega.eps,omega.Y),method="Richardson", method.args=list(),beta,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML)
    
    # AR.error.log.lik.fct(res$par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML)
    
    grad<-jacobian(func=AR.error.log.lik.fct,x=res$par,method="Richardson", method.args=list(),Y=Y,W=W,X=X,WplusWt=WplusWt,WWt=WWt,WY=WY,WX=WX,ARmod=ARmod,REML=REML)
    
      
    Info1<-cbind(Info.beta.beta,Info.beta.rho,Info.beta.omega.eps,Info.beta.omega.Y)
    Info2<-cbind(t(Info.beta.rho),Info.rho.omegas[1,1],Info.rho.omegas[1,2],Info.rho.omegas[1,3])
    Info3<-cbind(t(Info.beta.omega.eps),Info.rho.omegas[2,1],Info.rho.omegas[2,2],Info.rho.omegas[2,3]) 
    Info4<-cbind(t(Info.beta.omega.Y),Info.rho.omegas[3,1],Info.rho.omegas[3,2],Info.rho.omegas[3,3]) 
      
    Info<-rbind(Info1,Info2,Info3,Info4)
    
    
    #compare with 
    if(0){
    Info1<-hessian(f.all,c(beta,rho,omega.eps,omega.Y),method="Richardson", method.args=list(),Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,REML=REML)
    
    cat("Info\n")
    print(Info)
    cat("Info1\n")
    print(Info1)
    cat("dLdomegas\n")
    print(Info.rho.omegas)
    }
    
    Cov<-try(solve(Info))
    if(inherits(Cov,"try-error")){
      Cov<-try(ginv(Info))
    }
    
    if(!inherits(Cov,"try-error")){
    
    colnames(Info)<-colnames(Cov)<- par.names
    rownames(Info)<-rownames(Cov)<- par.names
    beta.se<-sqrt(diag(Cov[1:p,1:p]))
    rho.se<-sqrt(Cov[p+1,p+1])
    omega.Y.se<-sqrt(Cov[p+3,p+3])
    omega.eps.se<-sqrt(Cov[p+2,p+2])
    
    
    }
    })
    if(!inherits(g,"try-error")){
    
    
      cat("Calculation of Standard Errors finished\n")
    return(list(L=Lpos,grad=grad,omega.Y=omega.Y,omega.eps=omega.eps,beta=beta,rho=rho,
                Cov=Cov,Info=Info,beta.se=beta.se,rho.se=rho.se,omega.Y.se=omega.Y.se,omega.eps.se=omega.eps.se,time=st))
    
    }#end if(!inherits(g,"try-error")){
           
  }#if se
  
  return(list(L=Lpos,omega.Y=omega.Y,omega.eps=omega.eps,beta=beta,rho=rho,
              beta.se=NA,rho.se=NA,omega.Y.se=NA,omega.eps.se=NA,time=st))
  
}#end function AR.error




#############################################################
#             beginning of log-lik function
#############################################################
AR.error.log.lik.fct<-function(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=1,pr=TRUE,REML=T,n=dim(X)[1],p=dim(X)[2]){
 
#optimimise.theta<-TRUE  
  
  
  
  
  
  
  if(pr){
   rho<-par[1]
   
   #if(!optimise.theta){
   log.theta<-par[2]  #omega/omega_epsilon
   theta<-exp(log.theta)
   #}
   #cat("rho",rho,"theta",theta,"\n")
   
   }else{
    rho<-par[1]
    log.theta<-par[2]
    theta<-exp(log.theta) #theta
    omega<-par[3] #omega.eps
    beta<-par[4:(3+p)]     #beta
  }
    I<-Diagonal(n)
   
    #print(dim(Y))
    #print(dim(WplusWt))
    #print(dim(WWt))
    #print(dim(I))
    #print(rho)
    K<-I-rho*WplusWt+rho^2*WWt  #A*A' =(I-rho*W)*(I-rho*W)' I -rho*(W+W') +rho^2*(W*W')
    K<-forceSymmetric(K)
  
    if(ARmod==2){AX<-X-rho*WX}else{AX<-X}
    AY<- Y-rho*WY
    
 
  if(1){
    #faster
    #system.time({  
    chol.K <- Cholesky(K, super = TRUE, Imult = 0)
    #chol.K.theta <- Matrix:::ldetL2up(chol.K, K, theta)
    
    
    
    ldet.K<- 2*determinant(chol.K)$modulus #compare log(det(K))
  
  
      
    #f.theta<-function(theta){  
    chol.K.theta<- update(chol.K, K, mult=theta)
    ldet.K.theta<-2*determinant(chol.K.theta)$modulus
    logdetV<- ldet.K.theta - ldet.K  #compare with log(det(I+theta*solve((I-rho*t(W))%*%(I-rho*W))))
    
    b<- solve(chol.K.theta,AY,system = "A")
    C<- solve(chol.K.theta,AX,system = "A")
    #}) 
    }else{
  
    
    #system.time({
    chol.K.theta<-chol(K+theta*I,pivot=TRUE)
    ldet.K<-2*ldetA(rho,chol.list)
    logdetV<- 2*sum(log(diag(chol.K.theta))) - ldet.K
    b<- solve(t(chol.K.theta),AY)
    C<- solve(t(chol.K.theta),AX)
    #})
    }#end if else
    
    
    #compare CtC with XtMX
    # A<-I-rho*W;AtA<-t(A)%*%A; t(solve(A)%*%X)%*%solve(I+theta*solve(AtA))%*%(solve(A)%*%X)
    # t(AX)%*%solve(K+theta*I)%*%AX
    CtC<- t(AX)%*%C
    Ctb<- t(AX)%*%b  # t(AX)%*%solve(K+theta*I)%*%AY # compare with Ctb
  
    # compare with 
  
    btb<- t(AY)%*%b  # t(AY)%*%solve(K+theta*I)%*%AY # compare with btb
  
    
    # r <- expand(chol.K)
    # L.P <- with(r, crossprod(L,P))  ## == L'P
    # PLLP <- crossprod(L.P)          ## == (L'P)' L'P == P'LL'P  = XX = M M'
  
    if(pr){
      #calculate beta
    # V<-  I+theta*solve((I-rho*t(W))%*%(I-rho*W))
    # beta <- solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%Y  #SARerr
    beta<-solve(CtC)%*%Ctb  # solve(t(X)%*%solve(K+theta*I)%*%X)%*%t(X)%*%solve(K+theta*I)%*%AY
    #
    
    rtMr<-as.numeric(btb - t(Ctb)%*%beta)  
    # rtMr <- t(Y-X%*%beta)%*%solve(V)%*%(Y-X%*%beta)
    }else{
      #don't calculate beta
    rtMr<-as.numeric(btb - 2*t(Ctb)%*%beta + t(beta)%*%CtC%*%beta)   
    }
    
    # as.numeric(btb - 2* t(Ctb)%*%beta + t(beta)%*%CtC%*%beta )  
    # t(AY-X%*%beta)%*%solve(K+theta*I)%*%(AY-X%*%beta)
    XtMX<-forceSymmetric(CtC) 
    
    #XtMX <- t(X)%*%solve(V)%*%X
  
    if(pr){
      #calculate omega
    if(REML){
      n1<-n-p
      add<- 0.5*2*sum(log(diag(chol(XtMX))));omega<-rtMr/n1
    }else{
      n1<-n
      add<- 0;omega<-rtMr/n1
    }
    }else{
      #don't calculate omega
      if(REML){
        n1<-n-p
        add<- 0.5*2*sum(log(diag(chol(XtMX))));
      }else{
        n1<-n
        add<- 0;
      }  
      
      
    }
    
    #return(0.5*n*log(omega)   + 0.5*logdetV + 0.5*rtMr/omega + add)
    #}
    
    #}#end if optimise


#optimise(f.theta,interval=c(0+tol,1e4))
         
#print(L)
L<- c(0.5*n*log(2*pi) + 0.5*n1*log(omega)   + 0.5*logdetV + 0.5*rtMr/omega + add)   


attr(L, "omega_eps") <-   omega
attr(L, "beta") <- as.matrix(beta)
attr(L, "omega_Y")  <- theta*omega
attr(L,"rho")<-rho
attr(L,"CtC")<-CtC
attr(L,"Ctb")<-Ctb
attr(L,"rtMr")<-rtMr
#attr(L,"C")<-C
#attr(L,"b")<-b


#cat("hello")
if(exists("counter")){counter<<-counter+1;if(counter%%1==0){cat("Iteration:",counter," L:",L," rho:",rho," log.theta:",log.theta,"\n")}}


#print(L)
return(L)

#############################################################
} #             end of log-lik function of AR exact
#############################################################


