mrphd.model<- function(x,y){

  n=nrow(x)
  p=ncol(x)
  q=ncol(y)

  z<-stz.x(x) #standardize of X
  ny=stz.y(y) #standardize of Y

    pc<- prcomp(~.,data=as.data.frame(ny))
    pc.evt=pc$rotation
    pc.evl=(pc$sdev)^2
    pc.y=ny%*%pc.evt
    wt=pc.evl/sum(pc.evl)

##cov. matrix of Z and its square root of inverse
    szz<-covw.x(z,rep(1,n))
    sqzcov=solve(chol(szz))

    hess.m=matrix(0,nrow=p,ncol=p)
    sigma=matrix(0,nrow=p,ncol=p)

    for(i in 1:ncol(pc.y)){
        res=lm(pc.y[,i]~z)$residuals
        srzz<-covw.x(z,res)
        ker.mat=t(sqzcov)%*%srzz%*%sqzcov
        eg=svd(ker.mat)
        evl=eg$d
        evt=sqzcov%*%eg$u       #matrix of eigenvectors
        idz=order(abs(evl),decreasing = TRUE)
        evl<-abs(evl[idz])
        evt=evt[,idz]
        hess.m=evt%*%diag(evl)%*%t(evt)
        sigma=sigma+wt[i]*hess.m
     }
    if(sum(wt)!=0){
        sig.mat=t(sqzcov)%*%sigma%*%sqzcov
        eg.1=svd(sig.mat)
        evl=eg.1$d
        evt=sqzcov%*%eg.1$u
        idz=order(abs(evl),decreasing = TRUE)
        evl<-abs(evl[idz])
        evt=evt[,idz]
     }
  return(list(evt=evt,evl=evl))
}

#----------------------------------------------------------
mrphd.permutation.test<- function(x,y,evalues,evv,no.perm){
#no.perm=no of permutations,say 100.

  n=nrow(x)
  p=ncol(x)
  q=ncol(y)
  permut.pval=rep(0,p)

    perm.c=0
    total.evl=sum(evalues)
    for(i in 1:no.perm){
        v2x=as.matrix(x)%*%evv
        v2x=v2x[sample(n),]
        zz=mrphd.model(v2x,y)
        if(sum(zz$evl)>total.evl) perm.c=perm.c+1
    }
    permut.pval[1]=perm.c/no.perm

    for(i in 1:(p-1)){
       perm.c=0
       for(j in 1:no.perm){
          v1x=as.matrix(x)%*%evv[,1:i]
          v2x=as.matrix(x)%*%evv[,(i+1):p]
          v2x=v2x[sample(n),]
          zz=mrphd.model(cbind(v1x,v2x),y)
          if(sum(zz$evl[(i+1):p])>sum(evalues[(i+1):p])) perm.c=perm.c+1
       }
       permut.pval[i+1]=perm.c/no.perm
    }
  return(permut.pval)
}

#----------------------------------------------------------

#example.
#a=mrphd.model(x,y)
#b=mrphd.permutation.test(x,y,a$evl,a$evt,100)


