linefill=function(y,tidx,tfull)
{
 spl <- splinefun(tidx,y,method="natural")  
 spl(tfull)								 
}

getZhat=function(xt,pxl,ut,theta,phi,g,sg,omega.hat,avgt,KofNN,saveComps=F)
{
   no=nrow(xt)
   x.t=xt%*%theta
   g.1=g[,1];g.2=g[,2]
   ph1=phi[,1];ph2=phi[,2]
   sg1=sg[,1];sg2=sg[,2]
   ss=seq(1:10);lss=length(ss)  
   
   i=KofNN
      G1=rep(NA,no);G2=rep(NA,no)
      
      for(j in 1:no){
          sc.1=order(abs(pxl[,1]-x.t[j,1]))[1:ss[i]]
          G1[j]=mean(g.1[sc.1]*sg1[sc.1])
          sc.2=order(abs(pxl[,2]-x.t[j,2]))[1:ss[i]]
          G2[j]=mean(g.2[sc.2]*sg2[sc.2])
         }
      y.hat=cbind(G1,G2)%*%rbind(ph1,ph2)+omega.hat+avgt  
   
   if(!saveComps) return(y.hat)   else
   {
      #return(list(zhat=y.hat, xt=x.t, fs=cbind(G1,G2), w=rbind(ph1,ph2)))
	  v1=pxl[,1];v2=pxl[,2]
      v1v2=cbind(sort(v1),sort(v2))
      fs=cbind(g.1[order(v1)]*sg1[order(v1)],g.2[order(v2)]*sg2[order(v2)]) 	  
	  return(list(zhat=y.hat, variate=v1v2, fs=fs, w=rbind(ph1,ph2)))        
   }	  
}

stz.x<- function(x){   #cov. matrix of x (x: n by p)
   nobs=length(x[,1])
   x.c<-apply(x,2,mean)
   m.c<-x-rep(1,nobs)%o%x.c
   xcov=cov.x(as.matrix(x))
   sqxcov=solve(chol(xcov))
   stzx=as.matrix(m.c)%*%sqxcov
   return(stzx)  }

cov.x<- function(x){   #cov. matrix of x (x: n by p)
   x.c<-apply(x,2,mean)
   nobs=length(x[,1])
   m.c<-x-rep(1,nobs)%o%x.c
   cov=(as.matrix(t(m.c))%*%as.matrix(m.c))/nobs
   return(cov)  }
   
