
# This code is provided 'as is'. 
# Original code by Devin Johnson
# Jon Horne added code for Synoptic model 5-05-2010; jhorne@uidaho.edu 
################################################################################

## Synoptic Exponential Power Model #######################################################
seple = function(track,AvailList,locAvailFile, start.val = NULL)
{
require(MASS)
habmat = AvailList[[1]] #Just take the first one for now...

#Parameters are: meanx, meany,ln(scale), ln(shape), ThetaW1, ThetaW2, ...)
ubounds = c(rep(Inf,2),Inf,Inf, rep(Inf,ncol(habmat)-2))
lbounds = c(rep(-Inf,2),-Inf,-4.6, rep(-Inf,ncol(habmat)-2))

stime = Sys.time()
#If problems during optimization occur, try a different optimization method

#mle.sep = optim(start.val, sepLogLik, method = "Nelder-Mead",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=500, reltol = .00001))

mle.sep = optim(start.val, sepLogLik, method = "L-BFGS-B",upper = ubounds, lower = lbounds, hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=100))

etime = Sys.time()
Hessian = T
  if(Hessian){
    covmat = 2*ginv(mle.sep$hessian)
    se = sqrt(diag(covmat))
  } else {
    covmat = NULL
    se = rep(NA,length(mle.sep$par))
  }       
  z = mle.sep$par/se
  p.val = 2*(1-pnorm(abs(z)))
 
  parTable = cbind(Est = mle.sep$par, SE = se, Lower=mle.sep$par-1.96*se,
              Upper=mle.sep$par+1.96*se, Z=z, P = p.val)
  dimnames(parTable)[[1]] = c('mu.x','mu.y','ln.scale','ln.shape', colnames(habmat)[-c(1:2)])

#----------------------------------------------------------------------------------------
#Calculate probability of use distribution
UseDistList = list()
for (k in 1:length(AvailList)){
starttime <- proc.time()[3]

habmat = AvailList[[k]]

paramSEP = mle.sep$par
  meanx = paramSEP[1]
  meany = paramSEP[2]
  scale = exp(paramSEP[3])
  shape = exp(paramSEP[4])

# Calculate value of Gamma function evaluated at (1+shape/2)
coef = c(76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.00120858003, -0.536382 * 10^-5)
STP = 2.50662827465
FpF = 5.5
x = shape-1
TMP = x+FpF
TMP = (x+.5)*log(TMP)-TMP
SER = 1
for (i in 1:6) 
{
x = x+1
SER = SER+coef[i]/x
}
Gamma = exp(TMP+log(STP*SER))

#Calculate volume under non-normalized use function
maxy = max(habmat[,2])
maxx = max(habmat[,1])
tempX = habmat[,1]-maxx
tempX2 = 1/tempX
maxX2 = ginv(min(tempX2))+maxx
dimx = maxx-maxX2
tempY = habmat[,2]-maxy
tempY2 = 1/tempY
maxY2 = ginv(min(tempY2))+maxy
dimy = maxy-maxY2
cellsize = as.numeric(dimx*dimy)
MapDist = sqrt((habmat[,1]-meanx)^2+(habmat[,2]-meany)^2)
MapShapeFunc = exp(-(MapDist/scale)^(2/shape))
Normalize = 2/(shape*2*pi*(scale^2)*Gamma)
Map.g.a = MapShapeFunc*Normalize
# Selection Function in original paper ***do not use--not tested***
#  W = 1
#  if (length(paramSEP)==4){
#   wMap = 1
#  } else {
#   for (j in 5:length(paramSEP))
#   {
#   Wj = 1+habmat[,(j-2)]*paramSEP[j]
#   W = W*Wj
#   }
#   wMap = W
#  } #end else
# Exponential Selection Function
 if (length(paramSEP)==4){
  wMap = 1
 } else {
  if (length(paramSEP)==5){
   wMap = exp(habmat[,3]*paramSEP[5])
  } else {
   wMap = exp(habmat[,3:ncol(habmat)]%*%paramSEP[5:length(paramSEP)])
  } 
 }
MapNonNorm.g.u = Map.g.a*wMap
MapVolume = sum(MapNonNorm.g.u*cellsize)

#Calculate probability of use grid and cummulative probability grid
habmat = cbind(habmat,0)
colnames(habmat)[ncol(habmat)]= "Prob"
habmat[,ncol(habmat)] = MapNonNorm.g.u*cellsize/MapVolume
habmat =habmat[order(habmat[,ncol(habmat)],decreasing=T),]
habmat = cbind(habmat,0)
colnames(habmat)[ncol(habmat)]= "CumProb"
habmat[,ncol(habmat)]=cumsum(habmat[,(ncol(habmat)-1)])

UseDistList[k] = list(habmat)

endtime <- proc.time()[3]
runtime <- endtime - starttime
print("seple Calculate probability of use distribution loop time:")
print(runtime)
} #end AvailList loop


K = length(start.val)
AIC = mle.sep$value + 2*K
AICc = AIC+(2*K*(K+1)/(length(track)-K-1))

#Return a list of results
    list(parTable=(parTable),
          covmat = 2*ginv(mle.sep$hessian),
          Neg2xLikelihood = mle.sep$value,
	    AICc = AICc,
          evalTime = difftime(etime,stime),
          convergence=(mle.sep$convergence==0), UseDistList)
}

# =======================================================================================
# Likelihood Function
sepLogLik = function(paramSEP, track,AvailList)
{
  meanx = paramSEP[1]
  meany = paramSEP[2]
  scale = exp(paramSEP[3])
  shape = exp(paramSEP[4])

if (shape<0.01){  #For Nelder-Mead optimization which doesn't allow bounds on parameters
shape = 0.01}

# Calculate value of Gamma function evaluated at (1+shape/2)
coef = c(76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.00120858003, -0.536382 * 10^-5)
STP = 2.50662827465
FpF = 5.5
x = shape-1
TMP = x+FpF
TMP = (x+.5)*log(TMP)-TMP
SER = 1
for (i in 1:6)
{
x = x+1
SER = SER+coef[i]/x
}
Gamma = exp(TMP+log(STP*SER))

#Calculate volume under non-normalized use function
habmat = AvailList[[1]] # Just use the first one for now... no covariate values are used
maxy = max(habmat[,2])
maxx = max(habmat[,1])
tempX = habmat[,1]-maxx
tempX2 = 1/tempX
maxX2 = ginv(min(tempX2))+maxx
dimx = maxx-maxX2
tempY = habmat[,2]-maxy
tempY2 = 1/tempY
maxY2 = ginv(min(tempY2))+maxy
dimy = maxy-maxY2
cellsize = as.numeric(dimx*dimy)
MapDist = sqrt((habmat[,1]-meanx)^2+(habmat[,2]-meany)^2)
MapShapeFunc = exp(-(MapDist/scale)^(2/shape))
Normalize = 2/(shape*2*pi*(scale^2)*Gamma)
Map.g.a = MapShapeFunc*Normalize

SumLogLik = 0
LogLoc.g.u=array(0,nrow(track))

for (i in 1:nrow(track)){
  #get appropriate availability map
  availname = locAvailFile[i]
  habmat = AvailList[[availname]]
# Selection Function in original synoptic paper 
#  W = 1
#  if (length(paramSEP)==4){
#   wMap = 1
#  } else {
#   for (j in 5:length(paramSEP))
#   {
#   Wj = 1+habmat[,(j-2)]*paramSEP[j]
#   W = W*Wj
#   }
#   wMap = W
#  } 
# Exponential Selection Function
 if (length(paramSEP)==4){
  wMap = 1
 } else {
  if (length(paramSEP)==5){
   wMap = exp(habmat[,3]*paramSEP[5])
  } else {
   wMap = exp(habmat[,3:ncol(habmat)]%*%paramSEP[5:length(paramSEP)])
  }
 }

MapNonNorm.g.u = Map.g.a*wMap

MapVolume = sum(MapNonNorm.g.u*cellsize)

#Calculate Log-likelihood
# Selection function in original synoptic paper 
#  W = 1
#  if (length(paramSEP)==4){
#   wLoc = 1
#  } else {
#  for (j in 5:length(paramSEP))
#   {
#    Wj = 1+track[i,(j-1)]*paramSEP[j]
#    W = W*Wj
#   }
#  wLoc = W
#  }

# Exponential selection function
 if (length(paramSEP)==4){
  wLoc = 1
 } else {
  if (length(paramSEP)==5){
   wLoc = exp(track[i,-(1:3)]*paramSEP[5])
  } else {
   wLoc = exp(track[i,-(1:3)]%*%paramSEP[5:length(paramSEP)])
  }   
 }

LocDist = sqrt((track[i,1]-meanx)^2+(track[i,2]-meany)^2)
LocShapeFunc = exp(-(LocDist/scale)^(2/shape))
Loc.g.a = LocShapeFunc*Normalize
density = Loc.g.a*wLoc/MapVolume
if (is.na(density)){density=10^-320}
if (density==0){density=10^-320}
LogLoc.g.u[i] = log(density)
SumLogLik = SumLogLik+LogLoc.g.u[i]
} #end loop through locations
-2*SumLogLik
}

## Synoptic Bivariate Normal Model #######################################################

sbvnle = function(track,AvailList,locAvailFile,start.val = NULL)
{
require(MASS)

habmat = AvailList[[1]] #Just take the first one for now...

#Parameters are: meanx, meany,ln(stddev x), ln(stddev y), correlation, ThetaW1, ...)
ubounds = c(rep(Inf,2),rep(Inf,2),0.999, rep(Inf,ncol(habmat)-2))
lbounds = c(rep(-Inf,2),rep(-Inf,2),-0.999, rep(-Inf,ncol(habmat)-2))

stime = Sys.time()
#If problems during optimization occur, try a different optimization method

#mle.sep = optim(start.val, sepLogLik, method = "Nelder-Mead",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=500, reltol = .00001))

mle.sbvn = optim(start.val, sbvnLogLik, method = "L-BFGS-B", lower = lbounds, upper=ubounds,hessian = TRUE, track = track,AvailList=AvailList,control=list(maxit=100))
etime = Sys.time()

Hessian = T
  if(Hessian){
    covmat = 2*ginv(mle.sbvn$hessian)
    se = sqrt(diag(covmat))
  } else {
    covmat = NULL
    se = rep(NA,length(mle.sbvn$par))
  }       
  z = mle.sbvn$par/se
  p.val = 2*(1-pnorm(abs(z)))
  parTable = cbind(Est = mle.sbvn$par, SE = se, Lower=mle.sbvn$par-1.96*se,
              Upper=mle.sbvn$par+1.96*se, Z=z, P = p.val)
  dimnames(parTable)[[1]] = c('mu.x','mu.y','ln.stddev.x','ln.stddev.y', 'correlation',  colnames(habmat)[-c(1:2)])
#----------------------------------------------------------------------------------------
#Calculate probability of use distribution
UseDistList = list()
for (k in 1:length(AvailList)){
starttime <- proc.time()[3]
habmat = AvailList[[k]]

paramSBVN = mle.sbvn$par
  meanx = paramSBVN[1]
  meany = paramSBVN[2]
  sdx = exp(paramSBVN[3])
  sdy = exp(paramSBVN[4])
  corr = (paramSBVN[5])

#Calculate volume under non-normalized use function
maxy = max(habmat[,2])
maxx = max(habmat[,1])
tempX = habmat[,1]-maxx
tempX2 = 1/tempX
maxX2 = ginv(min(tempX2))+maxx
dimx = maxx-maxX2
tempY = habmat[,2]-maxy
tempY2 = 1/tempY
maxY2 = ginv(min(tempY2))+maxy
dimy = maxy-maxY2
cellsize = as.numeric(dimx*dimy)
Map.g.a.1 = 1/((2*pi)*sdx*sdy*sqrt(1-corr^2))
Map.g.a.2 = -1/(2*(1-corr^2))
Map.g.a.3 = (((habmat[,1]-meanx)/sdx)^2)+(((habmat[,2]-meany)/sdy)^2)-(2*corr*((habmat[,1]-meanx)/sdx)*((habmat[,2]-meany)/sdy))
Map.g.a = Map.g.a.1*exp(Map.g.a.2*Map.g.a.3) 

# Selection Function in original paper ***do not use--not tested***
#  W = 1
#  for (j in 6:length(paramSBVN))
#  {
#  Wj = 1+habmat[,(j-2)]*paramSBVN[j]
#  W = W*Wj
#  }
#  wMap = W
# Exponential Selection Function
 if (length(paramSBVN)==5){
  wMap = 1
 } else {
  if(length(paramSBVN)==6){
   wMap = exp(habmat[,3]*paramSBVN[6])
  } else {
   wMap = exp(habmat[,3:ncol(habmat)]%*%paramSBVN[6:length(paramSBVN)])
  }
 }
MapNonNorm.g.u = Map.g.a*wMap
MapVolume = sum(MapNonNorm.g.u*cellsize)

#Calculate probability of use grid and cummulative probability grid
habmat = cbind(habmat,0)
colnames(habmat)[ncol(habmat)]= "Prob"
habmat[,ncol(habmat)] = MapNonNorm.g.u*cellsize/MapVolume
habmat =habmat[order(habmat[,ncol(habmat)],decreasing=T),]
habmat = cbind(habmat,0)
colnames(habmat)[ncol(habmat)]= "CumProb"
habmat[,ncol(habmat)]=cumsum(habmat[,(ncol(habmat)-1)])

UseDistList[k] = list(habmat)

endtime <- proc.time()[3]
runtime <- endtime - starttime
print("seple Calculate probability of use distribution loop time:")
print(runtime)
} #edn AvailList loop
K = length(start.val)
AIC = mle.sbvn$value + 2*K
AICc = AIC+(2*K*(K+1)/(length(track)-K-1))

#Return a list of results
    list(parTable=(parTable),
          covmat = 2*ginv(mle.sbvn$hessian),
          Neg2xLikelihood = mle.sbvn$value,
          AICc = AICc,
          evalTime = difftime(etime,stime),
          convergence=(mle.sbvn$convergence==0), UseDistList)
}

# =======================================================================================
# Likelihood Function
sbvnLogLik = function(paramSBVN, track,AvailList)
{
  meanx = paramSBVN[1]
  meany = paramSBVN[2]
  sdx = exp(paramSBVN[3])
  sdy = exp(paramSBVN[4])
  corr = (paramSBVN[5])

#Calculate volume under non-normalized use function
habmat = AvailList[[1]]  #Just use the first one for now... no covariate values are used
maxy = max(habmat[,2])
maxx = max(habmat[,1])
tempX = habmat[,1]-maxx
tempX2 = 1/tempX
maxX2 = ginv(min(tempX2))+maxx
dimx = maxx-maxX2
tempY = habmat[,2]-maxy
tempY2 = 1/tempY
maxY2 = ginv(min(tempY2))+maxy
dimy = maxy-maxY2
cellsize = as.numeric(dimx*dimy)
Map.g.a.1 = 1/((2*pi)*sdx*sdy*sqrt(1-corr^2))
Map.g.a.2 = -1/(2*(1-corr^2))
Map.g.a.3 = (((habmat[,1]-meanx)/sdx)^2)+(((habmat[,2]-meany)/sdy)^2)-(2*corr*((habmat[,1]-meanx)/sdx)*((habmat[,2]-meany)/sdy))
Map.g.a = Map.g.a.1*exp(Map.g.a.2*Map.g.a.3) 

SumLogLik = 0
LogLoc.g.u=array(0,nrow(track))

for (i in 1:nrow(track)){
  #get appropriate availability map
  availname = locAvailFile[i]
  habmat = AvailList[[availname]]
# Selection Function in original paper  ***do not use--not tested***
#  W = 1
#  for (j in 6:length(paramSBVN))
#  {
#  Wj = 1+habmat[,(j-2)]*paramSBVN[j]
#  W = W*Wj
#  }
#  wMap = W
# Exponential Selection Function
 if (length(paramSBVN)==5){
  wMap = 1
 } else {
  if(length(paramSBVN)==6){
   wMap = exp(habmat[,3]*paramSBVN[6])
  } else {
   wMap = exp(habmat[,3:ncol(habmat)]%*%paramSBVN[6:length(paramSBVN)])
  }
 }
MapNonNorm.g.u = Map.g.a*wMap
MapVolume = sum(MapNonNorm.g.u*cellsize)

#Calculate Log-likelihood
# Selection function in original paper  ***do not use--not tested***
#  W = 1
#  for (j in 6:length(paramSBVN))
#  {
#  Wj = 1+track[,(j-1)]*paramSBVN[j]
#  W = W*Wj
#  }
#  wLoc = W
# Exponential selection function
 if (length(paramSBVN)==5){
  wLoc = 1
 } else {
  if(length(paramSBVN)==6){
   wLoc = exp(track[i,-(1:3)]*paramSBVN[6])
  } else {
   wLoc = exp(track[i,-(1:3)]%*%paramSBVN[6:length(paramSBVN)])
  }
 }

Loc.g.a.1 = 1/((2*pi)*sdx*sdy*sqrt(1-corr^2))
Loc.g.a.2 = -1/(2*(1-corr^2))
Loc.g.a.3 = (((track[i,1]-meanx)/sdx)^2)+(((track[i,2]-meany)/sdy)^2)-(2*corr*((track[i,1]-meanx)/sdx)*((track[i,2]-meany)/sdy))
Loc.g.a = Loc.g.a.1*exp(Loc.g.a.2*Loc.g.a.3) 
density = Loc.g.a*wLoc/MapVolume
if(is.na(density)){density=10^-320}
if(density==0){density=10^-320}
LogLoc.g.u[i]=log(density)
SumLogLik = SumLogLik+LogLoc.g.u[i]
} #end loop through locations
-2*SumLogLik
}
