#--------------------------------------------------------------
#load required packages and additional functions
#--------------------------------------------------------------
library(gstat)  
library(sp)
library(spacetime)
library(mc2d)
library(dr)
library(MASS)
library(pracma)
library(tsDyn)
library(mvtnorm)

source("func/fn_SIR.R")
source("func/fn_sir.program.R")
source("func/fn_PHD.R")
source("func/fn_mrsir.R")
source("func/fn_mrphd.R")
source("func/fn_local.linear.approx.R")
source("func/fn_local.linear.approx.p1.R")
source("func/fn_pde.curve1.R")
source("func/fn_pde.curve2.R")
source("func/fn_pde.curve2.mrphd.R")
source("func/fn_pde.curve2.mrsir.R")
source("func/fn_pde.curve2.mrsir.R")
source("func/stKrig.R")
source("func/func-extra.R")
source("func/PDE.R")
source("func/PDEplus.R")


#--------------------------------------------------------------
#simulate a learning set and a testing set
#--------------------------------------------------------------   
#simulation program for Example 1 is in "generator1.R"  
source("generator1.R")
set.seed(0) 
m1=generator1(n=100,TT=20)
r.tt=sample(2:19,4)  
r.case=sample(1:100,20)  
tidx=c(1:20)[-r.tt]   #observed time in learing
tfull=1:20	          #complete time for prediction
 
xx=cbind(m1$ss[,1],m1$ss[,2],m1$ss[,1]^2,m1$ss[,2]^2)
x=xx 
#x=stz.x(as.matrix(xx)) 
#use the above line for standardization when locations not within [-1,1]x[-1,1] or [0,1]x[0,1]
 
#learning set 
yl=m1$yy[-r.case,-r.tt]  
sl=m1$ss[-r.case,]  #learning set of locations
xl=x[-r.case,]      #learning set of spatial covariates

#testing set
yt=m1$yy[r.case,]   
st=m1$ss[r.case,]   #testing set of locations 
xt=x[r.case,]       #testing set of spatial covariates    


#--------------------------------------------------------------
#get predicted value from PDE or PDE+
#--------------------------------------------------------------   
hy=3; hx=0.5        #bandwidth setting
zhat_from_pde=PDE(xl,yl,xt,tidx,tfull,hx,hy,method="both")
zhat_pdeplus=PDEplus(xl,yl,xt,sl,st,tidx,tfull,hx,hy,method="both")  
#predicted value for z in size of ntxTT
#at nt testing locations and TT time points

#
#note for method in PDE() or PDEplus()
#
#method="phd": all directions using pe-mrPHD
#method="sir": all directions using mrSIR 
#method="both": directions using mrSIR and pe-mrPHD


#--------------------------------------------------------------
#visualize components from PDE+
#--------------------------------------------------------------   
cmps=PDEplus(xl,yl,xt,sl,st,tidx,tfull,hx,hy,method="both",saveComp=TRUE)  

x11()
par(mfrow=c(2,2))
plot(tfull,cmps$w[1,],xlab="time", ylab="basis for w1")
plot(tfull,cmps$w[2,],xlab="time", ylab="basis for w2")
plot(cmps$variate[,1],cmps$fs[,1],type="l", xlab="theta1'x", ylab="first scaled coefficient")
plot(cmps$variate[,2],cmps$fs[,2],type="l", xlab="theta2'x", ylab="second scaled coefficient")
