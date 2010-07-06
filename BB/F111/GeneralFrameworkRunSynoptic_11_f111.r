# General Framework for Analyzing Animal Space Use
# By: Devin Johnson, email: devin.johnson@noaa.gov; Synoptic model added by Jon Horne jhorne@uidaho.edu

# read in files:
#(1) set working directory with input files
workdir = choose.dir(caption="Select folder with input files")
setwd(workdir)

#(2) choose file containing functions: GeneralFramework_SynopticModels.R
source(file = "GeneralFrameWorkFitSynopticModels_11.r")

#(3) choose locations file:
locationsfile = "f111_summer_buffer_locs.txt"
 
origfilename = strsplit(locationsfile,"\\.")[[1]][1]
Track = as.matrix(read.table(file= locationsfile,head=T,sep=''))
 #separate text refering to availability file from track data
 locAvailFile=Track[,ncol(Track)]
 Track = apply(Track[,1:(ncol(Track)-1)],2,as.numeric)

#(4) choose a grid file for final probability of use distribution: 

UseDist = as.matrix(read.table(file="f111_msum02.txt",head=T,sep=''))
# added by Polly to remove first column in avail files which is cougar ID
UseDist = apply(UseDist[,2:(ncol(UseDist))],2,as.numeric)

#(5) specify which variables will be used: use '1' to specify a variable that will be
# used and '0' for variables that will not be used

ModelsList = list()
ModelsList[1] = list(rep(0, ncol(Track)-3))  # No covariates==> null model of space use

#Modify these to do different models...  '1' means use the covariate '0' leave out

ModelsList[2] = list(c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0))
ModelsList[3] = list(c(0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0))
ModelsList[4] = list(c(0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0))
ModelsList[5] = list(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0))
ModelsList[6] = list(c(0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0))
ModelsList[7] = list(c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0))
ModelsList[8] = list(c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0))
ModelsList[9] = list(c(0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0))
ModelsList[10] = list(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
ModelsList[11] = list(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0))
ModelsList[12] = list(c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0))
ModelsList[13] = list(c(0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0))
ModelsList[14] = list(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
ModelsList[15] = list(c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0))
ModelsList[16] = list(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
ModelsList[17] = list(c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0))
ModelsList[18] = list(c(1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0))
ModelsList[19] = list(c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0))
ModelsList[20] = list(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0))
ModelsList[21] = list(c(0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0))
ModelsList[22] = list(c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0))
ModelsList[23] = list(c(0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0))
ModelsList[24] = list(c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0))
ModelsList[25] = list(c(0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0))
ModelsList[26] = list(c(0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0))
ModelsList[27] = list(c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0))
ModelsList[28] = list(c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0))
ModelsList[29] = list(c(1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0))
ModelsList[30] = list(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0))
ModelsList[31] = list(c(1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0))
ModelsList[32] = list(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

#============================================================================
#For temporally changing habitat values, create a list of availability grids
AvailFileNames = table(locAvailFile)
AvailList=list()
for (i in 1:length(AvailFileNames)){
  #these two lines to test
  filename = paste(workdir,"\\",names(AvailFileNames)[i],sep="")
  habmat = as.matrix(read.table(file=filename,head=T,sep=''))
# added by Polly to remove first column in avail files which is cougar ID
  habmat = apply(habmat[,2:(ncol(habmat))],2,as.numeric)
  AvailList[i] = list(habmat)
}
names(AvailList)=names(AvailFileNames)
#====================================================================================
#Standardize Covariates to range 0 to 1; this helps with likelihood calculations
mins = matrix(NA,length(AvailList),ncol(AvailList[[1]]))
maxs=mins
for (i in 1:length(AvailList)){
  mins[i,]=apply(AvailList[[i]],2,min)
  maxs[i,]=apply(AvailList[[i]],2,max)
 }
colmin = apply(mins,2,min)
colmax = apply(maxs,2,max)
names(colmax) = c("mu.x", "mu.y", colnames(UseDist)[-c(1:2)])
CoVarMinMax = matrix(0,2,ncol(UseDist))
CoVarMinMax[1,]=colmin
CoVarMinMax[2,]=colmax
colnames(CoVarMinMax)=names(colmax)
rownames(CoVarMinMax)=c("min", "max")
#Write output minimum and maximum values for each input variable
outputfile = paste(origfilename,"_CoVar_MinMax.txt", sep = "")
data.frame(CoVarMinMax)
write.table (CoVarMinMax, file = outputfile)

for (i in 1:length(AvailList)){
 #Standardize Covariates in Availability grids
 for (col in 3:ncol(AvailList[[1]])){
  AvailList[[i]][,col]=(AvailList[[i]][,col]-colmin[col])/(colmax[col]-colmin[col])
 } # end column loop
} # end list loop

 #Standardize Covariates in Use distribution and location data
 for (col in 3:ncol(AvailList[[1]])){
  UseDist[,col]=(UseDist[,col]-colmin[col])/(colmax[col]-colmin[col])
  Track[,col+1]=(Track[,col+1]-colmin[col])/(colmax[col]-colmin[col])
 } # end column loop
#===================================================================================
#===================================================================================
# Loop through candidate models; Exponential power null model

PrevExpParamEsts = array(0,(ncol(UseDist)+2))
names(PrevExpParamEsts) = c("mu.x", "mu.y", "ln.scale", "ln.shape", colnames(UseDist)[-c(1:2)])

for (k in 1:length(ModelsList)){
    #delete columns (i.e., variables) in Availability grids, UseDist, and Track not used 
      Nvariables = sum(ModelsList[[k]]) #number of variables in current model
	CurrentAList=list()
      CurrentTrack=Track[,1:3] 					#keep x, y and time
	CurrentUseDist=UseDist[,1:2]					#keep x and y
	for (i in 1:length(AvailList)){
	 CurrentAList[[i]]=AvailList[[i]][,1:2]			#keep x and y
	} #end availability list loop
      names(CurrentAList)=names(AvailFileNames)     	
	cc=3
	for (col in 3:ncol(UseDist)){
	  if (ModelsList[[k]][col-2]==1){
          CurrentUseDist=cbind(CurrentUseDist,UseDist[,col])
          CurrentTrack=cbind(CurrentTrack, Track[,col+1])
          colnames(CurrentUseDist)[cc]=colnames(UseDist)[col]
          colnames(CurrentTrack)[cc+1]=colnames(UseDist)[col]
	    for (i in 1:length(AvailList)){
	       CurrentAList[[i]]=cbind(CurrentAList[[i]],AvailList[[i]][,col])
             colnames(CurrentAList[[i]])[cc]=colnames(UseDist)[col]	
	    } #end availability list loop  
		cc = cc+1                  
	  }# end if statement
	} #end variable column loop
    #-----------------------------------------------------------------------------------
### Synoptic with exponential power null model
# Get initial parameter values
ThetaW = c(rep(0, ncol(CurrentTrack)-3))	#Initial RSF coeff. set to 0; no selection 
if (k==1){
mu = as.numeric(apply(Track[,1:2],2,mean))
n = length(Track[,1])
scale = (sum(sqrt((Track[,1]-(rep(mu[1],n)))^2+(Track[,2]-(rep(mu[2],n)))^2)))/n
shape = 1
} else { #==> use estimated parameters of previous models (if they exist) for initial values
 mu = SEP.fit$parTable[1:2,1]
 scale = exp(SEP.fit$parTable[3,1])
 shape = exp(SEP.fit$parTable[4,1])
 names(ThetaW) = colnames(CurrentUseDist)[3:ncol(CurrentUseDist)]
 for (i in 1:length(PrevExpParamEsts)){
  for (j in 1:length(ThetaW)){
   if(names(PrevExpParamEsts)[i]==names(ThetaW)[j]){
    ThetaW[j]=PrevExpParamEsts[i]
   } #end if 
  } #end ThetaW loop
 } #end PrevExpParamEsts loop
} #end if

lnscale = log(scale)
lnshape = log(shape)

paramSEP = c(mu, lnscale, lnshape, ThetaW)
paramSEP = paramSEP[1:(4+length(ThetaW))]

SEP.fit = seple(CurrentTrack,CurrentAList,locAvailFile,CurrentUseDist,start.val=paramSEP)

PrevExpParamEsts[rownames(SEP.fit$parTable)]=SEP.fit$parTable[rownames(SEP.fit$parTable),1]

#Transform back parameter estimates for shape and scale parameters
 UnTransSEP.fit = SEP.fit
 UnTransSEP.fit$parTable[3:4,1:4]=exp(SEP.fit$parTable[3:4,1:4])
 rownames(UnTransSEP.fit$parTable)[3]="scale"
 rownames(UnTransSEP.fit$parTable)[4]="shape"


#Write output probability file to working directory
outputfile = paste(origfilename,"_ExpPower_Model",k,"_Prob.txt", sep = "")
write.table (UnTransSEP.fit[[7]], file = outputfile, col.names = TRUE, row.names = F, sep = ",")

#Write output table to working directory
outputfile = paste(origfilename,"_ExpPower_Model",k,"_Out.txt",sep = "")
tmp.wid = getOption("width")  # save current width
options(width=10000)
sink(outputfile)              # redirect output to file
print(UnTransSEP.fit[1:6])    # print the object
sink()                        # cancel redirection
options(width=tmp.wid)        # restore linewidth

} #end loop through candidate models list

#===================================================================================
#===================================================================================
# Loop through candidate models; Bivariate Normal null model

PrevBVNParamEsts = array(0,(ncol(UseDist)+3))
names(PrevBVNParamEsts) = c("mu.x", "mu.y", "sd.x", "sd.y", "corr.x.y",colnames(UseDist)[-c(1:2)])

for (k in 1:length(ModelsList)){
    #delete columns (i.e., variables) in Availability grids, UseDist, and Track not used 
      Nvariables = sum(ModelsList[[k]]) #number of variables in current model
	CurrentAList=list()

      CurrentTrack=Track[,1:3] 					#keep x, y and time
	CurrentUseDist=UseDist[,1:2]					#keep x and y
	for (i in 1:length(AvailList)){
	 CurrentAList[[i]]=AvailList[[i]][,1:2]			#keep x and y
	} #end availability list loop
      names(CurrentAList)=names(AvailFileNames)     	
	cc=3
	for (col in 3:ncol(UseDist)){
	  if (ModelsList[[k]][col-2]==1){
          CurrentUseDist=cbind(CurrentUseDist,UseDist[,col])
          CurrentTrack=cbind(CurrentTrack, Track[,col+1])
          colnames(CurrentUseDist)[cc]=colnames(UseDist)[col]
          colnames(CurrentTrack)[cc+1]=colnames(UseDist)[col]
	    for (i in 1:length(AvailList)){
	       CurrentAList[[i]]=cbind(CurrentAList[[i]],AvailList[[i]][,col])
             colnames(CurrentAList[[i]])[cc]=colnames(UseDist)[col]	
	    } #end availability list loop  
		cc = cc+1                  
	  }# end if statement
	} #end variable column loop
    #-----------------------------------------------------------------------------------
### Synoptic with bivariate normal
# Get initial parameter values
ThetaW = c(rep(0, ncol(CurrentTrack)-3))	#Initial RSF coeff. set to 0; no selection 
if (k==1){
mu = as.numeric(apply(Track[,1:2],2,mean))
sdx = sd(CurrentTrack[,1])
sdy = sd(CurrentTrack[,2])
corrXY = cor(CurrentTrack[,1],CurrentTrack[,2])

} else { #==> use estimated parameters of previous models (if they exist) for initial values
 mu = SBVN.fit$parTable[1:2,1]
 sdx = exp(SBVN.fit$parTable[3,1])
 sdy = exp(SBVN.fit$parTable[4,1])
 corrXY = (SBVN.fit$parTable[5,1])
 names(ThetaW) = colnames(CurrentUseDist)[3:ncol(CurrentUseDist)]
 for (i in 1:length(PrevBVNParamEsts)){
  for (j in 1:length(ThetaW)){
   if(names(PrevBVNParamEsts)[i]==names(ThetaW)[j]){
    ThetaW[j]=PrevBVNParamEsts[i]
   } #end if 
  } #end ThetaW loop
 } #end PrevBVNParamEsts loop
} #end if

lnsdx = log(sdx)
lnsdy = log(sdy)

paramSBVN = c(mu, lnsdx, lnsdy, corrXY, ThetaW)
paramSBVN = paramSBVN[1:(5+length(ThetaW))]

SBVN.fit = sbvnle(CurrentTrack,CurrentAList,locAvailFile,CurrentUseDist, start.val=paramSBVN)

PrevBVNParamEsts[rownames(SBVN.fit$parTable)]=SBVN.fit$parTable[rownames(SBVN.fit$parTable),1]

#Transform back parameter estimates for sdx and sdy
 UnTransSBVN.fit = SBVN.fit
 UnTransSBVN.fit$parTable[3:4,1:4]=exp(SBVN.fit$parTable[3:4,1:4])
 rownames(UnTransSBVN.fit$parTable)[3]="sdx"
 rownames(UnTransSBVN.fit$parTable)[4]="sxy"

#Write output probability file to working directory
outputfile = paste("filename,"_BVN_Model",k,"_Prob.txt", sep = "")
write.table (UnTransSBVN.fit[[7]], file = outputfile, col.names = TRUE, row.names = F, sep = ",")
#Write output table to working directory
outputfile = paste("BVN_Model",k,"_Out.txt",sep = "")
tmp.wid = getOption("width")  # save current width
options(width=10000)
sink(outputfile)              # redirect output to file
print(UnTransSBVN.fit[1:6])    # print the object
sink()                        # cancel redirection
options(width=tmp.wid)        # restore linewidth

} #end loop through candidate models list
