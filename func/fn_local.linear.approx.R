local.linear.approx<- function(x,y,h){ #for single direction. eg.,h_x=0.5 or h_y=3

   nobs=nrow(x)
   p=ncol(x)

     sxx<-covw.x(x,rep(1,nobs))
     beta=rep(1/sqrt(p),p)   #initial est. for beta 
     C=matrix(0,nrow=nobs,ncol=p)
     D<-array(0,dim=c(p,p,nobs))
     F=matrix(0,nrow=nobs,ncol=p)
     a=rep(0,nobs)
     d=rep(0,nobs)

     for(j in 1:nobs){
        kw<- mkerl.wt(x,x[j,],h)
        cf<- CFj(x,x[j,],y,kw)
        C[j,]<-cf$Cj
        F[j,]<-cf$Fj
        E<-kw%*%y
        D[,,j]<-Dj(x,x[j,],kw)
        m=rbind(cbind(1,t(C[j,])%*%beta),cbind(t(C[j,])%*%beta,t(beta)%*%D[,,j]%*%beta))
        m.inv<-ginv(m)
        ad=m.inv%*%rbind(E,t(beta)%*%F[j,])
        a[j]=ad[1]
        d[j]=ad[2]
     }
     if(length(d[d==0])<nobs){
        dD=0
        dFC=0
        for(j in 1:nobs){
           dD=dD+d[j]^2*D[,,j]
           dFC=dFC+d[j]*(F[j,]-a[j]*C[j,])
        }
        m1.inv<-ginv(dD)
        beta=m1.inv%*%dFC
     }
     loop.c=0  #loop count
     dis.beta=10
     repeat{
        loop.c=loop.c+1
#  print("MAVE count:")
#  print(loop.c)
        beta1=beta
        for(j in 1:nobs){
           kw<- kerl.wt(x,x[j,],beta,h)
           cf<- CFj(x,x[j,],y,kw)
           C[j,]<-cf$Cj
           F[j,]<-cf$Fj
           E<-kw%*%y
           D[,,j]<-Dj(x,x[j,],kw)
           m=rbind(cbind(1,t(C[j,])%*%beta),cbind(t(C[j,])%*%beta,t(beta)%*%D[,,j]%*%beta))
           m.inv<-ginv(m)
           ad=m.inv%*%rbind(E,t(beta)%*%F[j,])
           a[j]=ad[1]
           d[j]=ad[2]
        }
        dD=0
        dFC=0
        for(j in 1:nobs){
           dD=dD+d[j]^2*D[,,j]
           dFC=dFC+d[j]*(F[j,]-a[j]*C[j,])
        }
        m1.inv<-ginv(dD)
        beta=m1.inv%*%dFC
        dis.beta=norm(beta1-beta)
        if (dis.beta<0.001 || loop.c>=30) break #breaking rule
     }
     beta=beta/norm(beta)  #normalize beta
     phi=as.matrix(sxx)%*%as.matrix(beta)

#     print("MAVE estimation:")
#     print(beta)

  return(list(beta=beta,phi=phi,g=a))
}

