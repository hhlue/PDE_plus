##---------------------------------------------------------------
mkerl.wt<- function(x,xj,h){  #multidimensional kernel weight
   nobs=nrow(x)
   ww=rep(0,nobs)
      for(i in 1:nobs){
          q.inp=norm(x[i,]-xj)
          tt=q.inp/h
          if (q.inp<=h) {ww[i]=15*(1-tt^2)^2/16}
          if (q.inp>h)  {ww[i]=0}
#print(ww[i])
      }
   k.wt=ww/sum(ww)
   return(k.wt)  }

biwt.kerl<- function(x1,x2,h){  #biweight kernel
   dt=abs(x1-x2)
   tt=dt/h
     if (dt<=h) bk=15*(1-tt^2)^2/16
     if (dt>h)  bk=0
   return(bk)  }

kerl.wt<- function(x,xj,beta,h){  #kernel weight, given beta
   nobs=nrow(x)
   mx.b=x%*%beta
   xj.b=xj%*%beta
   ww=rep(0,nobs)
      for(i in 1:nobs){
          ww[i]=biwt.kerl(mx.b[i],xj.b,h)
      }
   kwt=ww/sum(ww)
   return(kwt)  }

Dj<- function(x,xi,w){  #Dj=Sigma_i[wij(xi-xj)(xi-xj)']
   nobs=nrow(x)  
   p=ncol(x)
   m.c<-x-rep(1,nobs)%o%xi
   mc.w<-m.c*rep(w,each=1,times=p)
   covw=t(as.matrix(mc.w))%*%as.matrix(m.c)
   return(covw)  }

CFj<- function(x,xi,y,w){  #Cj=Sigma_i[wij(xi-xj)] or Fj=Sigma_i[wij(xi-xj)yi]
   nobs=nrow(x)  
   p=ncol(x)
   m.c<-x-rep(1,nobs)%o%xi
   mc.w<-m.c*rep(w,each=1,times=p)
   mc.wy<-mc.w*rep(y,each=1,times=p)
   Cj=t(mc.w)%*%rep(1,nobs)
   Fj=t(mc.wy)%*%rep(1,nobs)
   return(list(Cj=Cj,Fj=Fj))  }

norm<- function(x){
   norm.x<-sqrt(sum(x^2))
   return(norm.x)  }

covw.x<- function(x,w){  #weighted covariance matrix of x (x: nxp)
   x.c<-apply(x,2,mean)
   nobs=nrow(x)  
   p=ncol(x)
   m.c<-x-rep(1,nobs)%o%x.c
   mc.w<-m.c*rep(w,each=1,times=p)
   covw=(as.matrix(t(mc.w))%*%as.matrix(m.c))/(nobs-1)
   return(covw)  }

cov.x<- function(x){   #cov. matrix of x (x: n by p)
   x.c<-apply(x,2,mean)
   nobs=length(x[,1])
   m.c<-x-rep(1,nobs)%o%x.c
   cov=(as.matrix(t(m.c))%*%as.matrix(m.c))/nobs
   return(cov)  }

stz.x<- function(x){   #cov. matrix of x (x: n by p)
   nobs=length(x[,1])
   x.c<-apply(x,2,mean)
   m.c<-x-rep(1,nobs)%o%x.c
   xcov=cov.x(as.matrix(x))
   sqxcov=solve(chol(xcov))
   stzx=as.matrix(m.c)%*%sqxcov
   return(stzx)  }

stz.y<- function(y){   #individual stdz for matrix of y (y: n by q)
   n=nrow(y)
   q=ncol(y)
   y.c<-apply(y,2,mean)
   yy=matrix(0,nrow=n,ncol=q)
      for(i in 1:q){
          yy[,i]=(y[,i]-rep(y.c[i],n))/sd(y[,i])
      }
   return(yy)  }

cos.angle<- function(v,beta){   #(v: length(beta) by no.iter)
   no.iter=ncol(v)
   cos.v=rep(0,no.iter)
      for(i in 1:no.iter){
          cos.v[i]=abs(beta%*%v[,i]/(norm(beta)*norm(v[,i])))
      }
#   mean.cos=mean(cos.v)
#   min.cos=min(cos.v)
#   sd.cos=sd(cos.v)
#   print(list(mean.cos,min.cos,sd.cos))
   return(cos.v)  }

rsqr.evt<- function(v,theta1,theta2,x){  #R^2
#     x=stz.x(x)
     t1.x=x%*%theta1
     t2.x=x%*%theta2
     vx=x%*%v
     tx=cbind(t1.x,t2.x)
     v.r2=summary(lm(vx~tx))$r.squared
#print(v.r2)
   return(v.r2)  }

align.basis<- function(f,ph1,ph2){
    ph1f=as.numeric(t(ph1)%*%f)
    ph2f=as.numeric(t(ph2)%*%f)
    ph11=as.numeric(t(ph1)%*%ph1)
    ph12=as.numeric(t(ph1)%*%ph2)
    ph22=as.numeric(t(ph2)%*%ph2)
    a1=(ph1f*ph22-ph2f*ph12)/(ph11*ph22-ph12*ph12)
    a2=(ph1f*ph12-ph2f*ph11)/(ph12*ph12-ph11*ph22)
    v1=a1*ph1+a2*ph2
    r.v1f=t(v1)%*%v1/t(f)%*%f
  return(r.v1f)  }

loop.align.basis<- function(f,m.phi){
   no.loop=ncol(m.phi[,,1])
   r.v1f=rep(0,no.loop)
      for(i in 1:no.loop){
          r.v1f[i]=align.basis(f,m.phi[,i,1],m.phi[,i,2])
      }
   return(r.v1f)  }

##---------------------------------------------------------------
##  example.
##  xx=stz.x(x)
##  apply(xx,2,mean)
##  diag(cov.x(as.matrix(xx)))
##---------------------------------------------------------------
#correlation
#cor.xy<- function(x,y) {
#   corr=cov(x,y)/(sd(x)*sd(y)) 
#   return(corr)  }
##---------------------------------------------------------------

