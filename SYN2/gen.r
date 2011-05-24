
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
 
  partable = cbind(Est = mle.sep$par, SE = se, Lower=mle.sep$par-1.96*se,
              Upper=mle.sep$par+1.96*se, Z=z, P = p.val)
  dimnames(partable)[[1]] = c('mu.x','mu.y','ln.scale','ln.shape', colnames(habmat)[-c(1:2)])

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
    list(partable=(partable),
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

#records optim results
OPT.RESULTS <- c()

## Synoptic Brownian Bridge Model #######################################################
synbbfit = function(track,start.val = NULL,k)
{
require(MASS)

#habmat = AvailList[[1]] #For finding number of covariates
#Parameters are: bbsd, ThetaW1, ... 
ubounds = c(Inf,rep(Inf,ncol(track)-4))
lbounds = c(-Inf,rep(-Inf,ncol(track)-4))
print("TS272")
print(start.val)

stime = Sys.time()
#If problems during optimization occur, try a different optimization method

#mle.sep = optim(start.val, sepLogLik, method = "Nelder-Mead",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=500, reltol = .00001))
print("TS277")
print(start.val)

#methods = "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"
#mle.synbb = optim(start.val, synbbLogLik, method = "L-BFGS-B", lower = lbounds, upper=ubounds,hessian = TRUE, track = track,k=k,control=list(maxit=100))
#mle.synbb = optim(start.val, synbbLogLik, method = "L-BFGS-B",hessian = TRUE, track = track,k=k,control=list(maxit=100))

OPT.RESULTS <- c()
maxit <- 20
#mle.synbb = optim(start.val, synbbLogLik, method = "BFGS",hessian = TRUE, track = track,k=k,control=list(maxit=100))
#mle.synbb = optim(start.val, synbbLogLik, method = "BFGS",hessian = TRUE, track = track,k=k,maxit=maxit,control=list(maxit=maxit))
mle.synbb = optim(start.val, synbbLogLik, method = "CG",hessian = TRUE, track = track,k=k,maxit=maxit,control=list(maxit=maxit,reltol = .01))
print("TS289")
print(OPT.RESULTS)


#mle.synbb = optim(start.val, sbvnLogLik, method = "BFGS", lower = lbounds, upper=ubounds,hessian = TRUE, track = track,AvailList=AvailList,control=list(maxit=100))
#mle.sep = optim(start.val, sepLogLik, method = "SANN",hessian = TRUE, track =track,AvailList=AvailList,control=list(maxit=26, reltol = .00001))

print("mle.synbb")
print(mle.synbb)

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
  partable = cbind(Est = mle.synbb$par, SE = se, Lower=mle.synbb$par-1.96*se,
              Upper=mle.synbb$par+1.96*se, Z=z, P = p.val)
  dimnames(partable)[[1]] = c('ln.sdbb', colnames(track)[-c(1:4)]) #fix for use on colnames instead
#----------------------------------------------------------------------------------------
#Calculate probability of use distribution

K = length(start.val)
AIC = mle.synbb$value + 2*K
AICc = AIC+(2*K*(K+1)/(length(track)-K-1))

#Return a list of results
    list(partable=(partable),
          covmat = 2*ginv(mle.synbb$hessian),
          Neg2xLikelihood = mle.synbb$value,
          AICc = AICc,
          evalTime = difftime(etime,stime),
          #convergence=(mle.synbb$convergence==0), UseDistList)
          convergence=(mle.synbb$convergence==0))

}

# =======================================================================================
# Likelihood Function
synbbLogLik = function(paramSYNBB, track, k,maxit,opt.results)
{
sdbb = exp(paramSYNBB[1])

#Calculate volume under non-normalized use function

#taking the cellsize from the first cellgrid, they should be all the same
cellsize <- cellgrid[1][[1]][[1]]$ma.cellsize 
print(paste("cellsize: ",cellsize))

SumLogLik = 0
SumLoc.g.u=array(0,nrow(track)) #i

#for time est.
print(paste("(1 - p) ~=",proc.time()[3]))
stime <- proc.time()[3]
get.sum.log.liks <- function(i) {
	lstime <- proc.time()[3]
	habmat = cellgrid[i][[1]][[2]][[k]] #clip by current triplicate
		
	loc.id <- i * 2
	#now done in cellgrid!	
	'

	StartX1 <- track[(loc.id - 1),1]
	StartY1 <- track[(loc.id - 1),2]
	LocX2 <- track[(loc.id),1]
	LocY2 <- track[(loc.id),2]
	EndX3 <- track[(loc.id + 1),1]
	EndY3 <- track[(loc.id + 1),2]
	StartSTD1 <- track[(loc.id - 1),4]
	EndSTD3 <- track[(loc.id + 1),4]
	StartTime1 <- track[(loc.id - 1),3]
	EndTime3 <- track[(loc.id + 1),3]
	TotalTime13 <- EndTime3 - StartTime1
	Time2 <- track[(loc.id),3] - track[(loc.id - 1),3]

	Alpha <- Time2 / TotalTime13
	MeanXTime2 <- StartX1 + ((Alpha) * (EndX3 - StartX1))
	MeanYTime2 <- StartY1 + ((Alpha) * (EndY3 - StartY1))
	VarTime2 <- TotalTime13 * Alpha * (1 - Alpha) * sdbb ^ 2 + (((1 - Alpha) ^ 2) * (StartSTD1 ^ 2)) + ((Alpha ^ 2) * (EndSTD3 ^ 2)) #sdbb gets changed in this maximization routine

	Map.g.a <<- matrix(0,nrow=nrow(habmat),ncol=3)
	for(j in 1:nrow(habmat)) {
		SqDist <- ((habmat[j,1] - MeanXTime2) ^ 2) + ((habmat[j,2] - MeanYTime2) ^ 2)
		PDFTime2 <- (1 / (2 * pi * VarTime2)) * exp(-0.5 * (SqDist / VarTime2))
		Map.g.a[j,3] <<- PDFTime2
		Map.g.a[j,1] <<- habmat[j,1]
		Map.g.a[j,2] <<- habmat[j,2]
	}
	'
	Map.g.a <<- cellgrid[i][[1]][[4]]
	LocX2 <- cellgrid[i][[1]][[1]]$LocX2
	MeanXTime2 <- cellgrid[i][[1]][[1]]$MeanXTime2
	LocY2 <- cellgrid[i][[1]][[1]]$LocY2
	MeanYTime2 <- cellgrid[i][[1]][[1]]$MeanYTime2
	TotalTime13 <- cellgrid[i][[1]][[1]]$TotalTime13
	Alpha <- cellgrid[i][[1]][[1]]$Alpha
	StartSTD1 <- cellgrid[i][[1]][[1]]$StartSTD1
	EndSTD3 <- cellgrid[i][[1]][[1]]$EndSTD3
	#VarTime2 <- cellgrid[i][[1]][[1]]$VarTime2
	VarTime2 <- TotalTime13 * Alpha * (1 - Alpha) * sdbb ^ 2 + (((1 - Alpha) ^ 2) * (StartSTD1 ^ 2)) + ((Alpha ^ 2) * (EndSTD3 ^ 2)) #sdbb gets changed in this maximization routine
	
	tsstime <- proc.time()[3]
	for(j in 1:nrow(habmat)) {
		SqDist <- ((habmat[j,1] - MeanXTime2) ^ 2) + ((habmat[j,2] - MeanYTime2) ^ 2)
		PDFTime2 <- (1 / (2 * pi * VarTime2)) * exp(-0.5 * (SqDist / VarTime2))
		Map.g.a[j,3] <<- PDFTime2
		Map.g.a[j,1] <<- habmat[j,1]
		Map.g.a[j,2] <<- habmat[j,2]
	}

	if(length(paramSYNBB) == 1) {
		wMap <- 1
		wLoc <- 1
		
	} else {
		if(length(paramSYNBB) == 2) {
			wMap <- exp(habmat[,3] * paramSYNBB[2])
			wLoc <- exp(track[loc.id,5] * paramSYNBB[2])
		} else {
			wMap <- exp(sum(habmat[,3:ncol(habmat)] * paramSYNBB[2:length(paramSYNBB)]))
			wLoc <- exp(sum(track[loc.id,5:ncol(track)] * paramSYNBB[2:length(paramSYNBB)]))

		}
	}
	
	MapNonNorm.g.u <- Map.g.a[,3] * wMap
	MapVolume <- sum(MapNonNorm.g.u * cellsize)

	#print(paste("trip:",i))
	#print("MapVolume")
	#print(MapVolume)
	#print("wMap")
	#print(wMap)
	#print("wLoc")
	#print(wLoc)

	#calculate log likelihood at middle location
	
	SqDist <- ((LocX2 - MeanXTime2) ^ 2) + ((LocY2 - MeanYTime2) ^ 2)
	Loc.g.a <- (1 / (2 * pi * VarTime2)) * exp(-0.5 * (SqDist / VarTime2))
	Loc.g.u <- (Loc.g.a * wLoc) / MapVolume
	#print("TS430")
	#print((Loc.g.a * wLoc) / MapVolume)

	Dist <- ((LocX2 - MeanXTime2) ^ 2) + ((LocY2 - MeanYTime2) ^ 2)
        Loc.g.a <- (1 / (2 * pi * VarTime2)) * exp(-0.5 * (SqDist / VarTime2))
        Loc.g.u <- (Loc.g.a * wLoc) / MapVolume

	if(is.na(Loc.g.u)) { 
		print(paste("Loc.g.u was na on triplicate",i,"model",k))
		print(Loc.g.a)
		print(wLoc)
		return(0)
	} #Loc.g.u <- 10^-320} #~ need to check if we ever get 

	letime <- proc.time()[3]
	lruntime <- letime - lstime
	lruntime <- lruntime * length(cellgrid)
	print(paste("p ~=",lruntime))

	return(Loc.g.u)
	#SumLogLik <- SumLogLik + SumLoc.g.u[i]
	#print("SumLogLik")
	#print(SumLogLik)
	#print("Loc.g.u")
	#print(log(Loc.g.u))
  } #end loop through triplicates

if(parallel) {
	print(paste("Getting SumLogLik in parallel"))
	clusterExport(c1,"cellgrid")
	LocLikVector <- unlist(parLapply(c1,1:length(cellgrid),get.sum.log.liks))
} else {
	LocLikVector <- unlist(lapply(1:length(cellgrid),get.sum.log.liks))
}
SumLocLik <- sum(LocLikVector)
PointLoc.g.u <- LocLikVector / SumLocLik
SumLogLik <- sum(log(PointLoc.g.u))

etime <- proc.time()[3]
runtime <- etime - stime
print(paste("Model",k," estimated runtime - ",runtime * maxit," seconds,",(runtime * maxit) / 60, "minutes"))
OPT.RESULTS <- OPT.RESULTS#,-2*SumLogLik)

print(paste("loglikelihood for model",k))
print(SumLogLik)
-2*SumLogLik
}

#get the area of the first square thats vertices are made by the first points
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

get.bb.var <- function(al) {
	#old code vs new code guide
	#SortedBBArray(x, 0) = x val
	#SortedBBArray(x, 1) = y val	
	#SortedBBArray(0, x) = first trip 
	#SortedBBArray(1, x) = second trip, etc
	#SortedBBArray(x, 2) = time	

	SumDistanceSqByN <- 0
	count <- length(al) / 3
	lapply(1:length(al), function (i) { #foreach trip
	 if(i %% 2 == 0) {
		#TimeI = SortedBBArray(I, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 2 of trip
		TimeI <- al$time[i] - al$time[i - 1]

		#TotalTimeI = SortedBBArray(I + 1, 2) - SortedBBArray(I - 1, 2) #time diff between 1 and 3
		TotalTimeI <-  al$time[i + 1] +  al$time[i - 1]

		#BBMeanXi = (TimeI / TotalTimeI) * (SortedBBArray(I + 1, 0) - SortedBBArray(I - 1, 0)) + SortedBBArray(I - 1, 0) #(time) * mean + 1st 
		BBMeanXi <- (TimeI / TotalTimeI) * (al$x[i + 1] - al$x[i - 1]) + al$x[i - 1]

		#BBMeanYi = (TimeI / TotalTimeI) * ((SortedBBArray(I + 1, 1) - SortedBBArray(I - 1, 1)) + SortedBBArray(I - 1, 1) #
		BBMeanYi <- (TimeI / TotalTimeI) * (al$y[i + 1] - al$y[i - 1]) + al$y[i - 1]
		
		#DistanceSqByN = ((SortedBBArray(I, 0) - BBMeanXi) ^ 2 + (SortedBBArray(I, 1) - BBMeanYi) ^ 2) / (count - 1) 
		DistanceSqByN <- ((al$x[i] - BBMeanXi) ^ 2 + (al$y[i] - BBMeanYi) ^ 2) / (count - 1)		
		
		BBDistanceSqByN = DistanceSqByN * (TotalTimeI / (TimeI * (TotalTimeI - TimeI)))
	  	SumDistanceSqByN <<- SumDistanceSqByN + BBDistanceSqByN
	  } #end if i %% 2 == 0
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
