library(gstat)
library(sp)
library(spacetime)


stKrig=function(yl,sl,st,Rtiny=FALSE,tidx,tfull, test.apart=F, 
    covM="productSum")
{ 
 TT=length(tfull) #TT=ncol(yl)
 Tl=ncol(yl)
 DataFit=data.frame(y=c(yl),
   s1=rep(sl[,1],Tl),
   s2=rep(sl[,2],Tl))
 
 now=Sys.time()
 Time=now+3600*(tidx)
 Time0=now+3600*(tfull)

 stfdf = STFDF(sp=SpatialPoints(sl), Time, DataFit)

 sp=SpatialPoints(st)
 si=SpatialPoints(sl)
 gridded(sp)=FALSE
 st0= STF(sp=sp, time=Time0)
 if(!test.apart) st1= STF(sp=si, time=Time0) else
 st1= STF(sp=si, time=Time)

 fx=y~1  
 if(!Rtiny)
 suppressWarnings(empVgm<-variogramST(fx,data=stfdf,tlags=0:8, progress=F,
            cutoff=1.3, width = 1.3/20))  else
 suppressWarnings( empVgm<-variogramST(fx,data=stfdf,tlags=0:8, progress=F,
            cutoff=max(dist(sl),na.rm=T)/3 ) )			
 vs=mean(apply(yl,1,var))
 vt=mean(apply(yl,2,var))			
 
 vModel <- vgmST(covM,
                 #method = "Nelder-Mead",
                 space=vgm(vs,"Exp", 1.5,0),
                 time =vgm(vt,"Exp", 2,0),	
				 joint=vgm(vs+vt,"Exp", 12, 0),	stAni=2	,		 
                 nugget=0.5, k=2,temporalUnit="hours")

 vfit <- try(Fit.StVariogram(empVgm, vModel))
 if(any(class(vfit)=="try-error"))
  {
  vModel <- vgmST(covM,
                 #method = "Nelder-Mead",
                 space=vgm(vs,"Exp", 1.5,0),
                 time =vgm(vt,"Exp", 2,0),
				 joint=vgm(vs+vt,"Exp", 5, 0), stAni=1,
				 nugget=1.5, k=2,temporalUnit="hours")
  vfit <- try(Fit.StVariogram(empVgm, vModel))
 }
 if(any(class(vfit)=="try-error"))
 {
  vModel <- vgmST(covM,
                 #method = "Nelder-Mead",
                 space=vgm(vs,"Exp", 1.5,0),
                 time =vgm(vt,"Exp", 2,0),
				 joint=vgm(vs+vt,"Exp", 1, 0), stAni=0.5,	
                 nugget=1.5,  k=2, temporalUnit="hours")
  vfit <- try(Fit.StVariogram(empVgm, vModel))
 }
 if(any(class(vfit)=="try-error"))
 {
   vModel <- vgmST(covM,
                   #method = "Nelder-Mead",
                   space=vgm(vs,"Exp", 3,0),
                   time =vgm(vt,"Exp", 4,0),
				   joint=vgm(vs+vt,"Exp", 1, 0), stAni=1,	
                   nugget=20,  k=2, temporalUnit="hours")
   vfit <- Fit.StVariogram(empVgm, vModel)
 }

 
 vModel=vfit; vModel$temporalUnit="hours"

 krig0=krigeST(y~1, data=stfdf, newdata=st0, 
				  modelList=vModel) 			  
 
 if(!test.apart)
 {
  yhat=matrix(unlist(krig0@data),ncol=TT)				  
  yhat1=yhat[1:nrow(sl),]
  yhat0=yhat[(nrow(sl)+1):nrow(st),]
  return(list(yhat1=yhat1,yhat0=yhat0))		
 }

 if(test.apart)
 { 
  krig1=krigeST(y~1, data=stfdf, newdata=st1, 
				  modelList=vModel) 			  
  yhat1=matrix(unlist(krig1@data),ncol=Tl)				  
  yhat0=matrix(unlist(krig0@data),ncol=TT)				  				  
  return(list(yhat1=yhat1,yhat0=yhat0))				  
 } 
}


Fit.StVariogram=function (object, model, ..., method = "L-BFGS-B", lower, upper, 
                          fit.method = 6, stAni = NA, wles) 
{
  if (!inherits(object, "StVariogram")) 
    stop("\"object\" must be of class \"StVariogram\"")
  if (!inherits(model, "StVariogramModel")) 
    stop("\"model\" must be of class \"StVariogramModel\".")
  sunit <- attr(object$spacelag, "units")
  tunit <- attr(object$timelag, "units")
  tu.obj = attr(model, "temporal unit")
  if (!is.null(tu.obj)) 
    stopifnot(identical(tunit, tu.obj))
  object <- na.omit(object)
  ret <- model
  if (!missing(wles)) {
    if (wles) 
      fit.method = 1
    else fit.method = 6
  }
  if (fit.method == 0) {
    attr(ret, "optim.output") <- "no fit"
    attr(ret, "MSE") <- mean((object$gamma - variogramSurface(model, 
                                                              data.frame(spacelag = object$dist, timelag = object$timelag))$gamma)^2)
    attr(ret, "spatial unit") <- sunit
    attr(ret, "temporal unit") <- tunit
    return(ret)
  }
  if ((fit.method == 7 | fit.method == 11) & is.null(model$stAni) & 
      is.na(stAni)) {
    message("[An uninformed spatio-temporal anisotropy value of '1 (spatial unit)/(temporal unit)' is automatically selected. Consider providing a sensible estimate for stAni or using a different fit.method.]")
    stAni <- 1
  }
  weightingFun <- switch(fit.method, function(obj, ...) obj$np, 
                         function(obj, gamma, ...) obj$np/gamma^2, function(obj, 
                                                                            ...) obj$np, function(obj, gamma, ...) obj$np/gamma^2, 
                         function(obj, ...) stop("fit.method = 5 (REML), is not yet implemented"), 
                         function(obj, ...) 1, function(obj, curStAni, ...) if (is.na(stAni)) obj$np/(obj$dist^2 + 
                                                                                                        (curStAni * obj$timelag)^2) else obj$np/(obj$dist^2 + 
                                                                                                                                                   (stAni * obj$timelag)^2), function(obj, ...) {
                                                                                                                                                     dist <- obj$dist
                                                                                                                                                     dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                                                                                                                                                     obj$np/dist^2
                                                                                                                                                   }, function(obj, ...) {
                                                                                                                                                     dist <- obj$timelag
                                                                                                                                                     dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                                                                                                                                                     obj$np/dist^2
                                                                                                                                                   }, function(obj, gamma, ...) 1/gamma^2, function(obj, 
                                                                                                                                                                                                    curStAni, ...) {
                                                                                                                                                     if (is.na(stAni)) 1/(obj$dist^2 + (curStAni * obj$timelag)^2) else 1/(obj$dist^2 + 
                                                                                                                                                                                                                             (stAni * obj$timelag)^2)
                                                                                                                                                   }, function(obj, ...) {
                                                                                                                                                     dist <- obj$dist
                                                                                                                                                     dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                                                                                                                                                     1/(obj$dist^2)
                                                                                                                                                   }, function(obj, ...) {
                                                                                                                                                     dist <- obj$timelag
                                                                                                                                                     dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                                                                                                                                                     1/(obj$timelag^2)
                                                                                                                                                   })
  if (is.null(weightingFun)) 
    stop(paste("fit.method =", fit.method, "is not implementend"))
  fitFun = function(par, trace = FALSE, ...) {
    curModel <- insertPar(par, model)
    gammaMod <- variogramSurface(curModel, data.frame(spacelag = object$dist, 
                                                      timelag = object$timelag))$gamma
    resSq <- (object$gamma - gammaMod)^2
    resSq <- resSq * weightingFun(object, gamma = gammaMod, 
                                  curStAni = curModel$stAni)
    if (trace) 
      print(c(par, MSE = mean(resSq)))
    mean(resSq)
  }
  if (missing(lower)) {
    min.s <- min(object$dist[object$dist > 0]) * 0.05
    min.t <- min(object$dist[object$timelag > 0]) * 0.05
    pos <- sqrt(.Machine$double.eps)
    lower <- switch(strsplit(model$stModel, "_")[[1]][1], 
                    separable = c(min.s, 0, min.t, 0, 0), productSum = c(0, 
                                                                         min.s, 0, 0, min.t, 0, pos), productSumOld = c(0, 
                                                                                                                        min.s, 0, 0, min.t, 0, 0), sumMetric = c(0, min.s, 
                                                                                                                                                                 0, 0, min.t, 0, 0, pos, 0, pos), simpleSumMetric = c(0, 
                                                                                                                                                                                                                      min.s, 0, min.t, 0, pos, 0, 0, pos), metric = c(0, 
                                                                                                                                                                                                                                                                      pos, 0, pos), stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  }
  if (missing(upper)) 
    upper <- switch(strsplit(model$stModel, "_")[[1]][1], 
                    separable = c(Inf, 1, Inf, 1, Inf), productSum = Inf, 
                    productSumOld = Inf, sumMetric = Inf, simpleSumMetric = Inf, 
                    metric = Inf, stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  pars.fit <- optim(extractPar(model), fitFun, ..., method = method, 
                    lower = lower+pos, upper = upper)
  ret <- insertPar(pars.fit$par, model)
  attr(ret, "optim.output") <- pars.fit
  attr(ret, "MSE") <- mean((object$gamma - variogramSurface(insertPar(pars.fit$par, 
                                                                      model), data.frame(spacelag = object$dist, timelag = object$timelag))$gamma)^2)
  attr(ret, "spatial unit") <- sunit
  attr(ret, "temporal unit") <- tunit
  return(ret)
}

pos <- sqrt(.Machine$double.eps)

insertPar=function (par, model) 
{
  switch(strsplit(model$stModel, "_")[[1]][1], separable = insertParSeparable(par, 
                                                                              model), productSum = gstat:::insertParProdSum(par, model), productSumOld = insertParProdSumOld(par, 
                                                                                                                                                                     model), sumMetric = insertParSumMetric(par, model), simpleSumMetric = gstat:::insertParSimpleSumMetric(par, 
                                                                                                                                                                                                                                                                    model), metric = insertParMetric(par, model), stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

insertParProdSum=function (par, model) 
{
  vgmST("productSum", space = vgm(par[1], as.character(model$space$model[2]), 
                                  par[2], par[3], kappa = model$space$kappa[2]), time = vgm(par[4], 
                                                                                            as.character(model$time$model[2]), par[5], par[6], kappa = model$time$kappa[2]), 
        k = par[7])
}


