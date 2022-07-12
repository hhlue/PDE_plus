PDEplus=function(xl,yl,xt,sl,st,tidx,tfull,hx,hy, method="both", saveComp=F,KofNN=3)
{  
   if(method=="both")   
    out=pde.curve2(xl,yl,10,hy,hx,100)
   if(method=="phd")   
    out=pde.curve2.mrphd(xl,yl,hy,hx)
   if(method=="sir")   
    out=pde.curve2.mrsir(xl,yl,3,hy,hx,100) 
    	
	nl=nrow(xl)
	itmax=1;iter=0
	Rtiny=( mad(unlist(yl))/abs(median(unlist(yl))) )< 4.5
	
    res.omega=matrix(NA,nl,length(tidx))
    for(j in 1:nl){res.omega[j,]=lm(yl[j,]~out$phi-1)$residuals }
    tavg=colMeans(res.omega)*0   #different from autoFRK
    
    
	fit=stKrig(res.omega-rep(1,nl)%o%tavg,sl,st,Rtiny,tidx,tfull,test.apart=T)
    omega=fit$yhat1	   
    nt=nrow(xt)
	
  while (iter<itmax){
   res=yl-omega

   if(method=="both")   
    out=pde.curve2(xl,res,10,hy,hx,100)
   if(method=="phd")   
    out=pde.curve2.mrphd(xl,res,hy,hx)
   if(method=="sir")   
    out=pde.curve2.mrsir(xl,res,3,hy,hx,100)	
    

    res.omega=matrix(NA,nl,length(tidx))
    for(j in 1:nl){res.omega[j,]=lm(yl[j,]~out$phi-1)$residuals }
    tavg=colMeans(res.omega)*0  #different from autoFRK

    fit=stKrig(res.omega-rep(1,nl)%o%tavg,sl,st,Rtiny,tidx,tfull,test.apart=T)
    omega=fit$yhat1    
    iter=iter+1
  }

    avgt=rep(1,nt)%o%rep(0,length(tfull))
    
    omega.hat=fit$yhat0						
     pxl=xl%*%out$theta 
     n1=nrow(xl)
     g1=out$g[,1];g2=out$g[,2]
     phi1=out$phi[,1];phi2=out$phi[,2]
     phi1=phi1/norm(phi1);phi2=phi2/norm(phi2) 
	 phit=cbind(linefill(phi1,tidx,tfull), linefill(phi2,tidx,tfull)) 
     sg1=rep(NA,n1);sg2=rep(NA,n1)
     for(k in 1:n1){phi=cbind(g1[k]*phi1,g2[k]*phi2)
                    coef=lm(yl[k,]~phi-1)$coefficients 
                    sg1[k]=coef[1];sg2[k]=coef[2] }
     sg=cbind(sg1,sg2)
	 
   zhat=getZhat(xt,pxl,ut=0,out$theta,phit,out$g,sg,omega.hat,avgt,KofNN, saveComps=saveComp)    
   return(zhat)   
}  