generator1=function(vxi=1,n=100,TT=20){

 ## n: no. of spatial locations  
 ## TT: time points 
 ## Output: nxTT matrix 
 
   mu=0
   f=list( f1=function(s) cos(((s[,1]+0.5)^2+(s[,2]+0.5)^2)*0.5*pi),
           f2=function(s) sin(((s[,1]-0.5)^2+(s[,2]-0.5)^2)*0.5*pi)
   )

   w=list( w1=function(TT) (0.5*1:TT-5)^2, 
           w2=function(TT) 5*sin((1:TT)/TT*2*pi)
   )   

   int.case=seq(-1,1,l=250)
   s=rand.case(int.case,n)  
   tj=seq(0,1,l=TT)
   K=length(f)
   wt=matrix(NA,TT,K)
   for(k in 1:K) wt[,k]=w[[k]](TT)
   y=z=matrix(NA,n,TT)
   mu=matrix(mu,n,TT)
   Fk=matrix(NA,n,K)   
   for(k in 1:K) Fk[,k]=f[[k]](s)
   uModel=gstat(formula=z~1, locations=~x+y,
                dummy=T, #True for unconditional simulation   
                beta=0,  #average value over the field                
                model=vgm(psill=1,range=2,nugget=0,model='Exp'),  
				        nmax=40)   				
   coord=data.frame(s); colnames(coord)=c("x","y")				
   for(tt in 1:TT)
   {
    z[,tt]=mu[,tt]+Fk%*%wt[tt,]+
	  predict(uModel, newdata=coord, nsim=1)$sim1*1
    y[,tt]=z[,tt]+0.5*rnorm(n,0,sd=sqrt(vxi)) 
   }	

 return(list(yy=y,ss=s,zz=z))
}
#--------------------------------------------------------------


rand.case<- function(int.case,cc){
   len=length(int.case)
   b=array(0,dim=c(len,2,len))
   for (j in 1:len){
     for (i in 1:len){
         b[i,,j]=c(int.case[j],int.case[i])
     }
   }
   c=rbind(b[,,1],b[,,2])
     for (i in 3:len){
         c=rbind(c,b[,,i])
     }
   case=sample(c(1:nrow(c)),cc)
   cc.x=c[case,]
 return(cc.x)}