
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
#print("TS nsyn: Map.g.a")
#print(Map.g.a)
#print("TS nsyn: wMap")
#print(wMap)

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
#print("seple Calculate probability of use distribution loop time:")
#print(runtime)
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

## Synoptic Brownian Bridge Model #######################################################
synbbfit = function(track,start.val = NULL,k)
{
require(MASS)

#habmat = AvailList[[1]] #For finding number of covariates
#Parameters are: bbsd, ThetaW1, ... 
ubounds = c(Inf,rep(Inf,ncol(track)-4))
print("TS269")
lbounds = c(-Inf,rep(-Inf,ncol(track)-4))
print(start.val)
print("TS272")
print(start.val)
print("TS274")
print(track)
print("TS276")

stime = Sys.time()
#If problems during optimization occur, try a different optimization method

#mle.sep = optim(start.val, sepLogLik, method = "Nelder-Mead",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=500, reltol = .00001))
print("TS278")
print(ubounds)
print("TS280")
print(lbounds)
print("TS281")
mle.synbb = optim(start.val, synbbLogLik, method = "L-BFGS-B", lower = lbounds, upper=ubounds,hessian = TRUE, track = track,k=k,control=list(maxit=100))
print("TS282")
#mle.synbb = optim(start.val, sbvnLogLik, method = "BFGS", lower = lbounds, upper=ubounds,hessian = TRUE, track = track,AvailList=AvailList,control=list(maxit=100))
#mle.sep = optim(start.val, sepLogLik, method = "SANN",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=26, reltol = .00001))

etime = Sys.time()

Hessian = T
  if(Hessian){
    covmat = 2*ginv(mle.synbb$hessian)
    se = sqrt(diag(covmat))
  } else {
    covmat = NULL
    se = rep(NA,length(mle.synbb$par))
  }       
  z = mle.synbb$par/se
  p.val = 2*(1-pnorm(abs(z)))
  parTable = cbind(Est = mle.synbb$par, SE = se, Lower=mle.synbb$par-1.96*se,
              Upper=mle.synbb$par+1.96*se, Z=z, P = p.val)
  dimnames(parTable)[[1]] = c('mu.x','mu.y','ln.stddev.x','ln.stddev.y', 'correlation',  colnames(habmat)[-c(1:2)])
#----------------------------------------------------------------------------------------
#Calculate probability of use distribution
UseDistList = list()
for (k in 1:length(AvailList)){
starttime <- proc.time()[3]
habmat = AvailList[[k]]

paramSBVN = mle.synbb$par
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
#print("seple Calculate probability of use distribution loop time:")
#print(runtime)
} #end AvailList loop
K = length(start.val)
AIC = mle.synbb$value + 2*K
AICc = AIC+(2*K*(K+1)/(length(track)-K-1))

#Return a list of results
    list(parTable=(parTable),
          covmat = 2*ginv(mle.synbb$hessian),
          Neg2xLikelihood = mle.synbb$value,
          AICc = AICc,
          evalTime = difftime(etime,stime),
          convergence=(mle.synbb$convergence==0), UseDistList)
}

# =======================================================================================
# Likelihood Function
synbbLogLik = function(paramSYNBB, track, k)
{
  sdbb = exp(paramSYNBB[1])


#Calculate volume under non-normalized use function

#taking the cellsize from the first cellgrid, they should be all the same
cellsize <- cellgrid[1][[1]][[1]]$ma.cellsize 
print(cellsize)
print(cellgrid[1])
q()

SumLogLik = 0
LogLoc.g.u=array(0,nrow(track)) #

extlen <- 0
if(syn.only) {
	extlen <- nrow(track)
}
else {
	extlen <- length(cellgrid)
}
#for (i in 1:nrow(track)){
for (i in length(cellgrid)) {
  habmat = cellgrid[i][[1]][[3]] #clip by current triplicate
  Map.g.a = #jon will send me stuff
  if(i %% 2 == 0) {
	# Exponential Selection Function
	if (length(paramSYNBB)==1){
	wMap = 1
	} else {
		if(length(paramSYNBB)==2){
		wMap = exp(habmat[,3]*paramSYNBB[2])
	} else {
		wMap = exp(habmat[,3:ncol(habmat)]%*%paramSYNBB[3:length(paramSYNBB)])
  	}
 	}
	MapNonNorm.g.u = Map.g.a*wMap
	MapVolume = sum(MapNonNorm.g.u*cellsize)
	#need to check:
	# when there are no covariates, that map volume =~ 1 (> .97), ie we didn't out more than ~3%

	#Calculate Log-likelihood
	# Exponential selection function
	if (length(paramSYNBB)==1){
		wLoc = 1
 	} else {
  		if(length(paramSYNBB)==2){
   		wLoc = exp(track[i,-(1:3)]*paramSYNBB[3])
  	} else {
		wLoc = exp(track[i,-(1:3)]%*%paramSYNBB[3:length(paramSYNBB)])
  	 }
 	}

	density <- get.bb.density(sdbb,MapVolume, wLoc) #we need to get this out of this LogLik func

	if(is.na(density)){density=10^-320}
	if(density==0){density=10^-320}
	LogLoc.g.u[i]=log(density)
	SumLogLik = SumLogLik+LogLoc.g.u[i]
	}
  } #end loop through locations
-2*SumLogLik
}

#get the area of the first square that's vertices are made by the first points
#in the data. The data should be spaced evenly for all points
get.ma.gridsize <- function(ma) {
	xdif <- ma$x[2] - ma$x[1]
	ydif <- 0
	for(i in 1:length(ma$y) - 1) {
		ydif <- ma$y[i] - ma$y[i + 1]
		if(length(ydif) != 0 && ydif != 0) {return(xdif*ydif)} #return the area
	}
	return(xdif*ydif) #return the area
}

get.bb.var <- function() {
	#old code vs new code guide
	#SortedBBArray(x, 0) = x val
	#SortedBBArray(x, 1) = y val	
	#SortedBBArray(0, x) = first trip 
	#SortedBBArray(1, x) = second trip, etc
	#SortedBBArray(x, 2) = time	

	SumDistanceSqByN <- 0
	count <- length(cellgrid)
	lapply(1:length(cellgrid), function (i) { #foreach trip
		#TimeI = SortedBBArray(I, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 2 of trip
		TimeI <- cellgrid[i][[1]][[3]]$time[2] - cellgrid[i][[1]][[3]]$time[1]

		#TotalTimeI = SortedBBArray(I + 1, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 3
		TotalTimeI <-  cellgrid[i][[1]][[3]]$time[3] +  cellgrid[i][[1]][[3]]$time[1]

		#BBMeanXi = (TimeI / TotalTimeI) * (SortedBBArray(I + 1, 0) - SortedBBArray(I - 1, 0)) + SortedBBArray(I - 1, 0) #(time) * mean + 1st 
		BBMeanXi <- (TimeI / TotalTimeI) * (cellgrid[i][[1]][[3]]$x[3] - cellgrid[i][[1]][[3]]$x[1]) + cellgrid[i][[1]][[3]]$x[1]

		#BBMeanYi = (TimeI / TotalTimeI) * ((SortedBBArray(I + 1, 1) - SortedBBArray(I - 1, 1)) + SortedBBArray(I - 1, 1) #
		BBMeanYi <- (TimeI / TotalTimeI) * (cellgrid[i][[1]][[3]]$y[3] - cellgrid[i][[1]][[3]]$y[1]) + cellgrid[i][[1]][[3]]$y[1]
		
		#DistanceSqByN = ((SortedBBArray(I, 0) - BBMeanXi) ^ 2 + (SortedBBArray(I, 1) - BBMeanYi) ^ 2) / (count - 1) 
		DistanceSqByN <- ((cellgrid[i][[1]][[3]]$x[2] - BBMeanXi) ^ 2 + (cellgrid[i][[1]][[3]]$y[2] - BBMeanYi) ^ 2) / (count - 1)		
		
		BBDistanceSqByN = DistanceSqByN * (TotalTimeI / (TimeI * (TotalTimeI - TimeI)))
		SumDistanceSqByN <<- SumDistanceSqByN + BBDistanceSqByN
	 }
	)
	BBVariance = SumDistanceSqByN
	return(BBVariance)
}

get.bb.density <- function(sdbb, MapVolume, wLoc) {
	bb.density <- 0
	count <- length(cellgrid)
	lapply(1:length(cellgrid), function (i) { #foreach trip
		#TimeI = SortedBBArray(I, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 2 of trip
		TimeI <- cellgrid[i][[1]][[3]]$time[2] - cellgrid[i][[1]][[3]]$time[1]

		#TotalTimeI = SortedBBArray(I + 1, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 3
		TotalTimeI <-  cellgrid[i][[1]][[3]]$time[3] +  cellgrid[i][[1]][[3]]$time[1]

		#BBMeanXi = (TimeI / TotalTimeI) * (SortedBBArray(I + 1, 0) - SortedBBArray(I - 1, 0)) + SortedBBArray(I - 1, 0) #(time) * mean + 1st 
		BBMeanXi <- (TimeI / TotalTimeI) * (cellgrid[i][[1]][[3]]$x[3] - cellgrid[i][[1]][[3]]$x[1]) + cellgrid[i][[1]][[3]]$x[1]

		#BBMeanYi = (TimeI / TotalTimeI) * ((SortedBBArray(I + 1, 1) - SortedBBArray(I - 1, 1)) + SortedBBArray(I - 1, 1) #
		BBMeanYi <- (TimeI / TotalTimeI) * (cellgrid[i][[1]][[3]]$y[3] - cellgrid[i][[1]][[3]]$y[1]) + cellgrid[i][[1]][[3]]$y[1]
		
		#DistanceSqByN = ((SortedBBArray(I, 0) - BBMeanXi) ^ 2 + (SortedBBArray(I, 1) - BBMeanYi) ^ 2) / (count - 1) 
		DistanceSq <- ((cellgrid[i][[1]][[3]]$x[2] - BBMeanXi) ^ 2 + (cellgrid[i][[1]][[3]]$y[2] - BBMeanYi) ^ 2)
		
		Alpha <- TimeI / TotalTimeI

		Stand <- TotalTimeI * Alpha * (1 - Alpha) * sdbb

		Loc.g.a <- 1/(2 * pi * Stand) * exp(-0.5 * (DistanceSq/Stand)) #1 over all?

		bb.density <<- Loc.g.a * wLoc / MapVolume #this probably doesn't go here, density is not just from the last iteration
	 }
	)
	return(bb.density)
}
