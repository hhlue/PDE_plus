pde.curve1<- function(x,y,slice,h1,h2){#eg. out=pde.curve1(x,y,10,3,0.5)

  n=nrow(x)
  p=ncol(x)
  q=ncol(y)
  nx=stz.x(x)  #standardize of X

     b=mrsir.model(x,y,slice)
     b1=mrphd.model(x,y)
     if(b$evl[1]>=b1$evl[1]) {theta.1=b$evt[,1];mr_evl=b$evl }
     if(b$evl[1]< b1$evl[1]) {theta.1=b1$evt[,1];mr_evl=b1$evl }
#no.perm=100
#perm_pval=mrphd.permutation.test(x,y,b$evl,b$evt,no.perm)
     u_t=apply(y,2,mean) #est. mean curve

     theta.1x=nx%*%theta.1
     d1=local.linear.approx(y,theta.1x,h1)
     phi.1=d1$phi

     coef=rep(0,n)
     dis.phi.0=10
     loop.c=0  #loop count
     repeat{
        loop.c=loop.c+1
        phi_1=phi.1
        for(j in 1:n){coef[j]=lm((y[j,]-u_t)~phi.1-1)$coef }
        a1=local.linear.approx(y,coef,h1)
        phi.1=a1$phi
        dis.phi=norm(phi.1-phi_1)
        if (dis.phi<=0.001 || dis.phi==dis.phi.0 || loop.c>=30) break #breaking rule
        dis.phi.0=dis.phi
     }
     dd.y=y%*%phi.1
     loop.c=0  #loop count
     repeat{
        loop.c=loop.c+1
        theta_1=theta.1
        d1=local.linear.approx.p1(x,dd.y,theta.1,h2)
        theta.1=d1$beta
        gj.1=d1$g
        dis.theta.1=norm(theta.1-theta_1)
        if (dis.theta.1<=0.001 || loop.c>=30) break #breaking rule
     }

return(list(theta1=theta.1,phi1=phi.1,mr.evl=mr_evl, #perm.pval=perm_pval,
g1=gj.1,ut=u_t))
}

#-----------------------------------------------------
#    muhat=apply(y-coef%*%t(varphi),2,mean) 
#-----------------------------------------------------

