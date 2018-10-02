library(Matrix)
library(numDeriv)
#source("ARmodels.sparse.add.fcts.r")


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
  #meth<-"BFGS"
  counter<<-0 
  cat("Optimisation L-BFGS-B\n")  
  st<-system.time(
  res <- try(optim(par, fn=AR.error.log.lik.fct, gr = NULL,Y=Y,W=W,X=X,WplusWt=WplusWt,WWt=WWt,WY=WY,WX=WX,ARmod=ARmod,pr=TRUE,REML=REML,
        method = meth,
        lower = c(-1+tol,-10), upper = c(1-tol,10),
        #control = list(), hessian = FALSE))
        control = list(factr=1e7,pgtol=0), hessian = FALSE))
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
} #             end of log-lik function of AR.error.log.lik.fct
#############################################################


