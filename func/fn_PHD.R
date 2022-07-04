phd.model<- function(x,y){
  
  n=nrow(x)
  p=ncol(x)

##cov. matrix of X and its square root of inverse
    sxx<-covw.x(x,rep(1,n))
    sqxcov=solve(chol(sxx))
    syxx=covw.x(x,y-mean(y))

    ker.mat=t(sqxcov)%*%syxx%*%sqxcov
    eg=svd(ker.mat)
    evl=eg$d
    evt=sqxcov%*%eg$u          #matrix of eigenvectors

    idx=order(abs(evl),decreasing = TRUE)
    evl<-abs(evl[idx])
    evt=evt[,idx]

    pvalues=rep(1,p)
    sq.ev=sum(evl^2)
    var.y=sd(y)^2

    for(k in 0:(p-1)){
        t=n*sq.ev/(2*var.y)  #test statistic
        pvalues[k+1]=pchisq(t,df=(0.5*(p-k)*(p-k+1)),lower.tail = FALSE) #p-value
        sq.ev=sq.ev-evl[k+1]^2
      }

    dim.x=length(pvalues[pvalues<0.05])
    if(dim.x>0) lam.bar<-sum(evl[1:dim.x])/sum(evl)
    if(dim.x==0) lam.bar<-0

  return(list(kermat=ker.mat,lambar=lam.bar,syxx=syxx,evt=evt,evl=evl,pvalues=pvalues))
}

