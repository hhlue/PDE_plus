pde.curve2.mrsir<- function(x,y,slice,h1,h2,no.perm){#eg. a=pde.curve2.mrsir(x,y,5,3,0.5,100)

  n=nrow(x)
  p=ncol(x)
  q=ncol(y)
  nx=stz.x(x)  #standardize of X

    b=mrsir.model(x,y,slice)
    mrsir_evl=b$evl
    theta.1=b$evt[,1]
    theta.2=b$evt[,2]
    u_t=apply(y,2,mean) #est. mean curve

      theta.1x=nx%*%theta.1
      theta.2x=nx%*%theta.2
      d1=local.linear.approx(y,theta.1x,h1)
      d2=local.linear.approx(y,theta.2x,h1)
      phi.1=d1$phi
      phi.2=d2$phi
      phi.2=phi.2-(dot(phi.2,phi.1)/sum(phi.1^2))*phi.1

     coef=matrix(0,nrow=n,ncol=2)
     loop.c=0  #loop count
     repeat{
        loop.c=loop.c+1
        phi_1=phi.1
        phi_2=phi.2
        phi=cbind(phi.1,phi.2)
        for(j in 1:n){coef[j,]=lm((y[j,]-u_t)~phi-1)$coef }
        a1=local.linear.approx(y,coef[,1],h1)
        a2=local.linear.approx(y,coef[,2],h1)
        phi.1=a1$phi
        phi.2=a2$phi
        phi.2=phi.2-(dot(phi.2,phi.1)/sum(phi.1^2))*phi.1
        z1=align.basis(phi_1,phi.1,phi.2)
        z2=align.basis(phi_2,phi.1,phi.2)
        if (loop.c>=30 || (z1>=0.95 & z2>=0.95)) break #breaking rule
     }
     dd.y=y%*%cbind(phi.1,phi.2)
     loop.c=0  #loop count
     repeat{
        loop.c=loop.c+1
        theta_1=theta.1
        a1=local.linear.approx.p1(x,dd.y[,1],theta.1,h2)
        theta.1=a1$beta
        gj.1=a1$g
        dis.theta.1=norm(theta.1-theta_1)
        if (dis.theta.1<0.001 || loop.c>=30) break #breaking rule
     }
     loop.c=0  #loop count
     repeat{
        loop.c=loop.c+1
        theta_2=theta.2
        a2=local.linear.approx.p1(x,dd.y[,2],theta.2,h2)
        theta.2=a2$beta
        theta.2=theta.2-(dot(theta.2,theta.1)/sum(theta.1^2))*theta.1
        gj.2=a2$g
        dis.theta.2=norm(theta.2-theta_2)
        if (dis.theta.2<0.001 || loop.c>=30) break #breaking rule
     }
     theta.1=theta.1/norm(theta.1);theta.2=theta.2/norm(theta.2)
     theta=cbind(theta.1,theta.2)
     phi=cbind(phi.1,phi.2)
     gj=cbind(gj.1,gj.2)

return(list(theta=theta,phi=phi,evl=mrsir_evl,g=gj,ut=u_t))
}

#-----------------------------------------------------
#    muhat=apply(y-coef%*%t(varphi),2,mean) 
#-----------------------------------------------------

