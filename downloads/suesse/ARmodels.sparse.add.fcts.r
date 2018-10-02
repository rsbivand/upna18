f_YinvX <- function(cholY,X) {
  X %>%
    forwardsolve(t(cholY),.) %>%
    backsolve(cholY,.) 
}

CleanEnvir <- function(pattern = "tmp"){
  rm(list = ls(envir=globalenv())[
    grep("tmp", ls(envir=globalenv()))], envir = globalenv())
}

combine.lists<-function(...){
  mapply(abind,...,along=1)
}#end combine lists


sparseDist <- function(m, k) {
  m <- t(m)
  n <- ncol(m)
  d <- vapply( seq_len(n-1L), function(i) { 
    d<-colSums((m[, seq(i+1L, n), drop=FALSE]-m[,i])^2)
    o<-sort.list(d, na.last=NA, method='quick')[seq_len(k)]
    c(sqrt(d[o]), o+i) 
  }, numeric(2*k)
  )
  dimnames(d) <- list(c(paste('d', seq_len(k), sep=''),
                        paste('i', seq_len(k), sep='')), colnames(m)[-n])
  d
}

get.W<-function(dist,nbs,k=10){
  n<-dim(dist)[2]+1
  number.neighbours.inv<-W<-Matrix(0,n,n)
  for(i in 1:(n-1)){
    is.nna<-!is.na(nbs[,i])
    ind<-nbs[is.nna,i][dist[is.nna,i]<k]
    if(length(ind>0)){
      W[i,ind]<-W[ind,i]<-1
    }
  }
  #normalise
  number.neighbours<-rowSums(W)
  diag(number.neighbours.inv)<-1/number.neighbours
  #dist<-apply(dist,1,sort)
  return(list(W=W,W.row=number.neighbours.inv%*%W,number.neighbours=number.neighbours))
}


calculate.W<-function(x,y,k){
  n<-length(x)
  dist<-number.neighbours.inv<-W<-Matrix(0,n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      dist[i,j]<-sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
      if(dist[i,j]<k){W[i,j]<-W[j,i]<-1}
    }#end
  }#end for i
  number.neighbours<-rowSums(W)
  diag(number.neighbours.inv)<-1/number.neighbours
  #dist<-apply(dist,1,sort)
  return(list(dist=dist,W=W,W.row=number.neighbours.inv%*%W,number.neighbours=number.neighbours)) 
}#end function

tr.prod<-function(A,B){
  sum(A*t(B))  
}

tr<-function(A){sum(diag(A))}

get.W.lists<-function(W,X=NULL,K,Xonly=FALSE){
  
  if(Xonly){
    
    WX.list<-list()
    WX.list[[1]]<-W%*%X  
    for(k in 1:(K-1)){
      WX.list[[k+1]]<-W%*%WX.list[[k]]     
    }
    
    return(list(WX.list=WX.list))  
    
    
  }else{  
    
    
    WplusWt.list<-list()
    WX.list<-list()
    W.list<-list()
    
    
    if(K>=1){
      W1W1t<-W%*%t(W)
      W1<-W  
      WplusWt.list[[1]]<-W+t(W)
      if(!is.null(X)){
        WX.list[[1]]<-W%*%X
        WX.list[[2]]<-W%*%WX.list[[1]]
      }
      W.list[[1]]<-W
      
      cat("K=1\n")
      
    }
    
    if(K>=2){
      W2<-W%*%W
      W.list[[2]]<-W2
      W1W2t<-(W%*%t(W2))
      W2W2t<-W%*%W1W2t 
      WplusWt.list[[2]]<-W2+t(W2)
      if(!is.null(X)){
        WX.list[[3]]<-W%*%WX.list[[2]]
        WX.list[[4]]<-W%*%WX.list[[3]]
      }
      cat("K=2\n")
    }
    if(K>=3){
      W3<-W%*%W2 
      W.list[[3]]<-W3
      W1W3t<-W%*%t(W3)
      W2W3t<-W%*%W1W3t
      W3W3t<-W%*%W2W3t  
      WplusWt.list[[3]]<-W3+t(W3) 
      if(!is.null(X)){
        WX.list[[5]]<-W%*%WX.list[[3]]
        WX.list[[6]]<-W%*%WX.list[[4]]
      }
      cat("K=3\n")
    }
    
    if(K>=4){
      W4<-W%*%W3
      W.list[[4]]<-W4
      W1W4t<-W%*%t(W4)
      W2W4t<-W%*%W1W4t
      W3W4t<-W%*%W2W4t
      W4W4t<-W%*%W3W4t   
      W2W3t+t(W2W3t)+W1W4t+t(W1W4t)
      WplusWt.list[[4]]<-W4+t(W4)
      if(!is.null(X)){
        WX.list[[7]]<-W%*%WX.list[[6]]
        WX.list[[8]]<-W%*%WX.list[[7]]
      }
      cat("K=4\n")
    }
    
    if(K>=5){
      W5<-W%*%W4
      W.list[[5]]<-W5
      W1W5t<-W%*%t(W5)  
      W2W5t<-W%*%W1W5t
      W3W5t<-W%*%W2W5t
      W4W5t<-W%*%W3W5t
      W5W5t<-W%*%W4W5t
      WplusWt.list[[5]]<-W5+t(W5)
      if(!is.null(X)){
        WX.list[[9]]<-W%*%WX.list[[8]]
        WX.list[[10]]<-W%*%WX.list[[9]]
      }
      cat("K=5\n")
    }
    
    if(K==10){
      W6<-W%*%W5
      W7<-W%*%W6
      W8<-W%*%W7
      W9<-W%*%W8
      W10<-W%*%W9
      
      W1W6t<-W%*%t(W6)  
      W2W6t<-W%*%W1W6t
      W3W6t<-W%*%W2W6t
      W4W6t<-W%*%W3W6t
      W5W6t<-W%*%W4W6t
      W6W6t<-W%*%W5W6t
      
      cat("K=6\n")
      
      WplusWt.list[[6]]<-W6+t(W6)
      
      W1W7t<-W%*%t(W7)  
      W2W7t<-W%*%W1W7t
      W3W7t<-W%*%W2W7t
      W4W7t<-W%*%W3W7t
      W5W7t<-W%*%W4W7t
      W6W7t<-W%*%W5W7t
      W7W7t<-W%*%W6W7t
      WplusWt.list[[7]]<-W7+t(W7)
      
      cat("K=7\n")
      
      W1W8t<-W%*%t(W8)  
      W2W8t<-W%*%W1W8t
      W3W8t<-W%*%W2W8t
      W4W8t<-W%*%W3W8t
      W5W8t<-W%*%W4W8t
      W6W8t<-W%*%W5W8t
      W7W8t<-W%*%W6W8t
      W8W8t<-W%*%W7W8t 
      WplusWt.list[[8]]<-W8+t(W8)
      
      cat("K=8\n")
      
      W1W9t<-W%*%t(W9)  
      W2W9t<-W%*%W1W9t
      W3W9t<-W%*%W2W9t
      W4W9t<-W%*%W3W9t
      W5W9t<-W%*%W4W9t
      W6W9t<-W%*%W5W9t
      W7W9t<-W%*%W6W9t
      W8W9t<-W%*%W7W9t 
      W9W9t<-W%*%W8W9t 
      WplusWt.list[[9]]<-W9+t(W9)
      
      cat("K=9\n")
      
      W1W10t<-W%*%t(W10)  
      W2W10t<-W%*%W1W10t
      W3W10t<-W%*%W2W10t
      W4W10t<-W%*%W3W10t
      W5W10t<-W%*%W4W10t
      W6W10t<-W%*%W5W10t
      W7W10t<-W%*%W6W10t
      W8W10t<-W%*%W7W10t 
      W9W10t<-W%*%W8W10t 
      W10W10t<-W%*%W9W10t 
      WplusWt.list[[10]]<-W10+t(W10)
      
      if(!is.null(X)){
        WX.list[[11]]<-W%*%WX.list[[10]]
        WX.list[[12]]<-W%*%WX.list[[11]]
        WX.list[[13]]<-W%*%WX.list[[12]]
        WX.list[[14]]<-W%*%WX.list[[13]]
        WX.list[[15]]<-W%*%WX.list[[14]]
        WX.list[[16]]<-W%*%WX.list[[15]]
        WX.list[[17]]<-W%*%WX.list[[16]]
        WX.list[[18]]<-W%*%WX.list[[17]]
        WX.list[[19]]<-W%*%WX.list[[18]]
        WX.list[[20]]<-W%*%WX.list[[19]]
      }
      
      cat("K=10\n")
      
    }
    
    
    switch(K,
           #2
           WWt.list<-list(W1W1t)  #K=1
           ,
           #2  3  4
           WWt.list<-list(W1W1t,W1W2t+t(W1W2t),W2W2t) #K=2
           ,
           #2   3              4                               5            6
           WWt.list<-list(W1W1t,W1W2t+t(W1W2t),W2W2t+W1W3t+t(W1W3t),W2W3t+t(W2W3t),W3W3t)  #K=3
           ,
           #2   3              4                               5                         6                       7              8
           WWt.list<-list(W1W1t,W1W2t+t(W1W2t),W2W2t+W1W3t+t(W1W3t),W2W3t+t(W2W3t)+W1W4t+t(W1W4t),W3W3t+W2W4t+t(W2W4t),W3W4t+t(W3W4t),W4W4t)  #K=4
           ,
           #2   3              4                               5                         6                                  7                               8                  9                10
           WWt.list<-list(W1W1t,W1W2t+t(W1W2t),W2W2t+W1W3t+t(W1W3t),W2W3t+t(W2W3t)+W1W4t+t(W1W4t),W3W3t+W2W4t+t(W2W4t)+W1W5t+t(W1W5t),W3W4t+t(W3W4t)+W2W5t+t(W2W5t),W4W4t+W3W5t+t(W3W5t),W4W5t+t(W4W5t),W5W5t)  #K=5
           ,
{
  #K=6    
}
,
{
  #K=7  
}
,
{
  #K=8  
}
,
{
  #K=9  
}
,
{
  #K=10 
  #   2     #K=3            # 4                  # 5                          #6                                    #7                                             #8                                                 
  WWt.list<-list(W1W1t,W1W2t+t(W1W2t),W2W2t+W1W3t+t(W1W3t),W2W3t+t(W2W3t)+W1W4t+t(W1W4t),W3W3t+W2W4t+t(W2W4t)+W1W5t+t(W1W5t),W3W4t+t(W3W4t)+W2W5t+t(W2W5t)+W1W6t+t(W1W6t),W4W4t+W3W5t+t(W3W5t)+W2W6t+t(W2W6t)+W1W7t+t(W1W7t))    
  #i+j=9
  WWt.list[[8]]<- W4W5t+t(W4W5t)+W3W6t+t(W3W6t)+W2W7t+t(W2W7t)+W1W8t+t(W1W8t)
  #i+j=10 
  WWt.list[[9]]<- W5W5t+W4W6t+t(W4W6t)+W3W7t+t(W3W7t)+W2W8t+t(W2W8t)+W1W9t+t(W1W9t)
  #i+j=11 
  WWt.list[[10]]<- W5W6t+t(W5W6t)+W4W7t+t(W4W7t)+W3W8t+t(W3W8t)+W2W9t+t(W2W9t)+W1W10t+t(W1W10t)
  #i+j=12  
  WWt.list[[11]]<- W6W6t+W5W7t+t(W5W7t)+W4W8t+t(W4W8t)+W3W9t+t(W3W9t)+W2W10t+t(W2W10t) 
  #i+j=13  
  WWt.list[[12]]<- W6W7t+t(W6W7t)+W5W8t+t(W5W8t)+W4W9t+t(W4W9t)+W3W10t+t(W3W10t) 
  #i+j=14  
  WWt.list[[13]]<- W7W7t+W6W8t+t(W6W8t)+W5W9t+t(W5W9t)+W4W10t+t(W4W10t)
  #i+j=15  
  WWt.list[[14]]<- W7W8t+t(W7W8t)+W6W9t+t(W6W9t)+W5W10t+t(W5W10t)
  #i+j=16 
  WWt.list[[15]]<- W8W8t+W7W9t+t(W7W9t)+W6W10t+t(W6W10t)
  #i+j=17 
  WWt.list[[16]]<- W8W9t+t(W8W9t)+W7W10t+t(W7W10t)
  #i+j=18 
  WWt.list[[17]]<- W9W9t+W8W10t+t(W8W10t)
  #i+j=19  
  WWt.list[[18]]<- W9W10t+t(W9W10t)
  #i+j=20 
  WWt.list[[19]]<- W10W10t
}
    )

return(list(WplusWt.list=WplusWt.list,WWt.list=WWt.list,WX.list=WX.list,W.list=W.list))
  }
}


make.sample.list<-function(List,s1,s2){
  new.list<-list()  
  for(i in 1:length(List)){
    new.list[[i]]<-List[[i]][s1,s2]
  }
  return(new.list)
}



ldetA<-function(rho,chol.list,tol=1e-6){ 
#ldetA(0.5,chol.list)
if(rho>tol){
  #based on I-rho*W = rho*(-W + 1/rho *I) 
  logdetArho <- chol.list$n * log(rho) + Matrix:::ldetL2up(chol.list$Cn, chol.list$Wsimn, 1/rho) #positive rho
}else if(rho< -tol){
  #based on I-rho*W = -rho*( W + 1/(-rho) *I) 
  logdetArho <- chol.list$n * log(-rho) + Matrix:::ldetL2up(chol.list$Cp, chol.list$Wsimp, -1/rho) #negative rho
}else{ logdetArho<- 0 }
return(logdetArho)
}#end function



output<-function(mat,row.names=NULL,digits=3,nsmall=2){
hilf<-dim(mat)
m<-hilf[1]
n<-hilf[2]

if(is.null(row.names)){row.names<-rownames(mat)}
output1<-NULL
for(i in 1:m){
  output1<-c(output1,row.names[i])
for(j in 1:n){
if(!is.na(mat[i,j])){  
 output1<-c(output1," & ",format(mat[i,j] ,nsmall=nsmall,digits=digits))
}else{
  output1<-c(output1," & -- ") 
}
 
}  
output1<-c(output1,"\\\\\n")
  
}#end   
cat(output1,"\n")
}#end

output1<-function(mat,row.names=NULL,digits=3,nsmall=2){
  hilf<-dim(mat)
  m<-hilf[1]
  n<-hilf[2]
  
  if(is.null(row.names)){row.names<-rownames(mat)}
  output1<-NULL
  for(i in 1:m){
    output1<-c(output1,row.names[i])
    for(j in 1:n){
      if(!is.na(mat[i,j])){  
        if(j%%2==1){
        output1<-c(output1," & $",format(mat[i,j] ,nsmall=nsmall,digits=digits),"$")
        }
        if(j%%2==0){
        output1<-c(output1,"  ($",format(mat[i,j] ,nsmall=nsmall,digits=digits),"$)")
        } 
        
        }else{
        output1<-c(output1," & -- ") 
      }
      
    }  
    output1<-c(output1,"\\\\\n")
    
  }#end   
  cat(output1,"\n")
}#end


expit<-function(x){ex<-exp(x);return(ex/(1+ex))}


####################################################################################
# plotting
plot.AR.error<-function(Y,W,X,ARmod=1,REML=T,log.theta.lim=c(-10,10),tol=1e-3){
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
  
  log.theta.lim<-c(-8,-2)
  rho.lim<-c(0.95,1)
  tol<-2e-3
  log.theta<-seq(log.theta.lim[1],log.theta.lim[2],by=(log.theta.lim[2]-log.theta.lim[1])*tol/2) #from -10 to +10
  rho<-seq(rho.lim[1]+tol,rho.lim[2]-tol,by=(rho.lim[2]-rho.lim[1])*tol/2) # from -1 to +1
  
  z<-matrix(NA,length(log.theta),length(rho))
  
  for(i in 1:length(log.theta)){
   for(j in 1:length(rho)){
      par<-c(rho[j],log.theta[i])
     z[i,j]<- AR.error.log.lik.fct(par,Y,W,X,WplusWt,WWt,WY,WX,ARmod=ARmod,pr=TRUE,REML=REML,n=n,p=p)
   if(j%%10==0){
     cat("rho=",rho[j],"  log.theta=",log.theta[i]," L=",z[i,j],"\n")
     
   }# end if
     }#end for j
  }#end for i
  
  save(z,file="z1.RData")
  
  return(list(z=z,log.theta=log.theta.ind,rho=rho.ind))
  
}#end function  
  