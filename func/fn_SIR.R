sir.model<- function(x,y,slice){

  n=nrow(x)
  p=ncol(x)

  c.slice=rep(floor(n/slice),slice) #No. of cases in slices 
  if((n%%slice)!=0){
    for(j in slice-(n%%slice)+1:slice)c.slice[j]=c.slice[j]+1
  }
  x.mean=apply(x,2,mean)
  
##cov. matrix of X and its square root of inverse
    xcov=cov.x(as.matrix(x))
    sqxcov=solve(chol(xcov))
##  sqxcov=expm::sqrtm(solve(xcov))

    data=as.data.frame(cbind(y,x))  #putting y and x.pr together
    data=data[order(data[,1]),]     #sort data with order of y

##compute slice mean for each slice (H by p)
    s.mean=matrix(0,ncol=p,nrow=slice)
    for(j in 1:slice){ 
      t=sum(c.slice[1:j])         #determine the last position of each slice
      h=data[(t-c.slice[j]+1):t,] #catch the slice from raw data
      s.mean[j,]=apply(h[,2:(p+1)],2,mean) #compute slice mean
    }

##compute the slice cov. matrix
    scov=matrix(0,ncol=p,nrow=p)
    for(j in 1:slice) {
      s.mean[j,]=s.mean[j,]-x.mean 
      scov=scov+(c.slice[j]/n)*s.mean[j,]%*%t(s.mean[j,])
    }

##Chi-squared test for no. of sig. dir.
##H0 : no. of sig. dir. = k
##H1 : no. of sig. dir. > k

    ker.mat=t(sqxcov)%*%scov%*%sqxcov
    eg=svd(ker.mat)
    evl=eg$d
    evt=sqxcov%*%eg$u          #matrix of eigenvectors
    pvalues=rep(1,min(p,(slice-1)))
    for(k in 0:(p-1)){
      if((slice-k-1)>0){       #control the test if no. of slices < dim. of X
        t=n*sum(eg$d[(k+1):p]) #test statistic 
        pvalues[k+1]=pchisq(t,df=((p-k)*(slice-k-1)),lower.tail = FALSE) #p-value
        }
      }

    dim.x=length(pvalues[pvalues<0.05])
    if(dim.x>0) lam.bar<-sum(eg$d[1:dim.x])/sum(eg$d)
    if(dim.x==0) lam.bar<-eg$d[1]

  return(list(kermat=ker.mat,lambar=lam.bar,scov=scov,evt=evt,evl=evl,pvalues=pvalues))
}

