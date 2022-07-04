mrsir.model<- function(x,y,slice){

  n=nrow(x)
  p=ncol(x)
  q=ncol(y)

  nx=stz.x(x) #standardize of X

##cov. matrix of X and its square root of inverse
    xcov=cov.x(as.matrix(nx))
    sqxcov=solve(chol(xcov))

    km=array(0,dim=c(p,p,q)) #kernel matrix
    lam.bar=rep(0,q) #ratio of sum of nonzero eigenvalues to total eigenvalues
    lam.m=matrix(0,nrow=p,ncol=q) #store eigenvalues for each response
    wk.m=matrix(0,nrow=p,ncol=p)

##mrSIR algorithm
    for(i in 1:q){
      out<- sir.model(nx,y[,i],slice)
      km[,,i]=out$scov
      lam.bar[i]=out$lambar
      lam.m[,i]=out$evl
    }
    if (sum(lam.bar)>0) {lam.bar=lam.bar/sum(lam.bar)} 
    else {lam.bar=lam.m[1,]/sum(lam.m[1,])}
#print(lam.bar)   #detecting for piecewise significance

    for(i in 1:q){
      wk.m=wk.m+lam.bar[i]*km[,,i]
    }

    ker.mat=t(sqxcov)%*%wk.m%*%sqxcov
    eg=svd(ker.mat)
    evl=eg$d
    evt=sqxcov%*%eg$u 

  return(list(evt=evt,evl=evl))
}

