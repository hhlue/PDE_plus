PDE=function(xl,yl,xt,tidx,tfull,hx,hy,KofNN=3, 
    method="both", saveComp=F)
{  
   if(method=="both")   
    out=pde.curve2(xl,yl,10,hy,hx,100)
   if(method=="phd")   
    out=pde.curve2.mrphd(xl,yl,hy,hx)
   if(method=="sir")   
    out=pde.curve2.mrsir(xl,yl,3,hy,hx,100)
   	
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
   tavg=rep(0,length(tfull))
   avgt=rep(1,nrow(xt))%o%tavg
   omega.hat=matrix(0,nrow(xt),length(tfull))
   zhat=getZhat(xt,pxl,ut=0,out$theta,phit,out$g,sg,omega.hat,avgt,KofNN, saveComps=saveComp)    
   return(zhat)   
}  