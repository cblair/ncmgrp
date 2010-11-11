# Synoptic Model for Analyzing Animal Space Use
# By: Jon Horne; jhorne@uidaho.edu

locAvailFile <- c()

syn <- function(al, ma) {
	#(3) choose locations file:
	#locationsfile = paste(cluster, "_all_locations.txt", sep="") 
 
	#origfilename = strsplit(locationsfile,"\\.")[[1]][1]
	#origfilename = paste(synbb.outfile,"-",min(al$x),"-",max(al$x),"-",min(al$y),"-",max(al$y), sep="")
	origfilename = paste(synbb.outfile,"-",sep="")
	#Track = as.matrix(read.table(file= locationsfile,head=T,sep=''))
	Track = as.matrix(al)
	#separate text refering to availability file from track data
	#locAvailFile <<- Track[,ncol(Track) - 1]
	locAvailFile <<- al$ExtentFile
	Track = apply(Track[,1:(ncol(Track)-1)],2,as.numeric)
	#(4) specify which variables will be used: use '1' to specify a variable that will be
	# used and '0' for variables that will not be used

	modfname = paste(workdir,"/models.r",sep="")
	source(modfname)

	#============================================================================
	#For temporally changing habitat values, create a list of availability grids
	AvailFileNames = table(locAvailFile)
	AvailList=list()
	for (i in 1:length(AvailFileNames)){
		#filename = paste(workdir,"/",names(AvailFileNames)[i],sep="")
		#habmat = as.matrix(read.table(file=filename,head=T,sep=''))
		habmat = as.matrix(ma)
		habmat = apply(habmat[,1:(ncol(habmat))],2,as.numeric)
		AvailList[i] = list(habmat)
	}
	names(AvailList)=names(AvailFileNames)
	#====================================================================================	
	#Standardize Covariates to range 0 to 1; this helps with likelihood calculations
	mins = matrix(NA,length(AvailList),ncol(AvailList[[1]]))
	maxs=mins
	#Enter absolute minimum and maximum values for each covariate from all of the input availability grids
	#The first 2 columns are for x and y coordinates and mins and maxs can be set to 0, these values will 
	#be ingored
	#Auto calculating min max added 07/01/10
	minmaxfilename = paste(workdir,"/",cluster,"_master_avail.txt",sep="")
	minmax = read.table(minmaxfilename, sep='\t',header=TRUE, as.is=TRUE) 
	#minmax <- ma
	#change to calculate off of every ma submitted, not just current one
	colmin <- c()
	colmax <- c()
	for(i in 1:ncol(minmax)) {
		colmin <- c(colmin, min(minmax[,i]))
		colmax <- c(colmax, max(minmax[,i]))
	}
	#end auto calc min max
	names(colmax) = c("mu.x", "mu.y", colnames(AvailList[[1]])[-c(1:2)])
	CoVarMinMax = matrix(0,2,ncol(AvailList[[1]]))
	CoVarMinMax[1,]=colmin
	CoVarMinMax[2,]=colmax
	colnames(CoVarMinMax)=names(colmax)
	rownames(CoVarMinMax)=c("min", "max")
	
	#Write output file of minimum and maximum values for each input variable
	outputfile = paste(synbb.outfile,"-",origfilename,"_CoVar_MinMax.txt", sep = "")
	data.frame(CoVarMinMax)
	write.table (CoVarMinMax, file = outputfile)
		
	#colmin = c(0,0,0,0,300)
	#colmax = c(0,0,4000,89,3000)
	#Standardize Covariates in Availability grids
	for (i in 1:length(AvailList)){
		for (col in 3:ncol(AvailList[[1]])){
			AvailList[[i]][,col]=(AvailList[[i]][,col]-colmin[col])/(colmax[col]-colmin[col])
		} # end column loop
	} # end list loop

	#Standardize Covariates in location data
	for (col in 3:ncol(AvailList[[1]])){
	Track[,col+2]=(Track[,col+2]-colmin[col])/(colmax[col]-colmin[col])
 	} # end column loop

	#===================================================================================
	#===================================================================================
	# Loop through candidate models;

	PrevSYNBBParamEsts = array(0,(ncol(AvailList[[1]])-1))
	names(PrevSYNBBParamEsts) = c("bbsd",colnames(AvailList[[1]])[-c(1:2)])
	#for synbb mode:
	#names(PrevBVNParamEsts) = c("bb.var",colnames(AvailList[[1]])[-c(1:2)])
	for(k in 1:length(ModelsList)) {
		#print("ts91")
		#print(k)
		lapply(1:length(cellgrid), function(i) {
		print(paste("TS94:",i))
		#lapply(1, function(i) {
			cellgrid[i][[1]][[1]] <<- cellgrid[i][[1]][[2]][1:2]
			z <- 1
		 	#because we want to ignore x and y cols
			for(j in ModelsList[[k]]) {
				if(ModelsList[[k]][z] == 1) {
					cellgrid[i][[1]][[1]] <<- merge(cellgrid[i][[1]][[1]], cellgrid[i][[1]][[2]][z+2], all.x=TRUE)
				}
				z <- z + 1
			}
		} ) #end insert filtered ma in each cellgrid
	}
	print(cellgrid[1][[1]])
	q()

	for (k in 1:length(ModelsList)){
    		#delete columns (i.e., variables) in Availability grids, and Track not used 
      		Nvariables = sum(ModelsList[[k]]) #number of variables in current model
		CurrentAList=list()
      		CurrentTrack=Track[,1:4] 					#keep x, y, time, sd
		for (i in 1:length(AvailList)){
	 	CurrentAList[[i]]=AvailList[[i]][,1:2]			#keep x and y
		} #end availability list loop
      		names(CurrentAList)=names(AvailFileNames)     	
		cc=3
		for (col in 3:ncol(AvailList[[1]])){
	  		if (ModelsList[[k]][col-2]==1){
          			CurrentTrack=cbind(CurrentTrack, Track[,col+2])
				colnames(CurrentTrack)[cc+1]=colnames(AvailList[[i]])[col]
	    			for (i in 1:length(AvailList)){
					CurrentAList[[i]]=cbind(CurrentAList[[i]],AvailList[[i]][,col])
             				colnames(CurrentAList[[i]])[cc]=colnames(AvailList[[1]])[col]	
	    			} #end availability list loop  
				cc = cc+1                  
	  		}# end if statement
		} #end variable column loopi
		print("TS112")
		#Kill each cellgrid ma col that we don't need
  		
		#-----------------------------------------------------------------------------------
		### Synoptic with brownian bridge
		# Get initial parameter values
		ThetaW = c(rep(0, ncol(CurrentTrack)-4))	#Initial RSF coeff. set to 0; no selection 
		if (k==1){
			#get initial bb standard
			bbsd = sqrt(bb.var)
			print("TSxxx")
		} else { #==> use estimated parameters of previous models (if they exist) for initial values
		print("TSxxx2")
			#get new estimate of bb.var
			mu = SBVN.fit$parTable[1:2,1]
			sdx = exp(SBVN.fit$parTable[3,1])
			bbsd = exp(SYNBB.fit$partable[1,1])
			sdy = exp(SBVN.fit$parTable[4,1])
			corrXY = (SBVN.fit$parTable[5,1])
			names(ThetaW) = colnames(CurrentAList[[1]])[3:ncol(CurrentAList[[1]])]
			for (i in 1:length(PrevBVNParamEsts)){
				for (j in 1:length(ThetaW)){
					if(names(PrevBVNParamEsts)[i]==names(ThetaW)[j]){
					ThetaW[j]=PrevBVNParamEsts[i]
	   				} #end if 
  				} #end ThetaW loop
 			} #end PrevBVNParamEsts loop
		} #end if

		lnbbsd = log(bbsd)
		paramSYNBB = c(lnbbsd, ThetaW)

		#print out each arg to compare with data
		print("TS145")
		print(CurrentAList)
		print("TS147")
		print(locAvailFile)

		#
		SYNBB.fit = synbbfit(CurrentTrack,CurrentAList,locAvailFile, start.val=paramSYNBB)
		#change sbvnle to synbble, names only

		#PrevBVNParamEsts[rownames(SBVN.fit$parTable)]=SBVN.fit$parTable[rownames(SBVN.fit$parTable),1]
		print("TS165")
		PrevSYNBBParamEsts[rownames(SYNBB.fit$parTable)]=SYNBB.fit$parTable[rownames(SYNBB.fit$parTable),1]

		#Transform back parameter estimates for sdx and sdy
		print("TS169")
		UnTransSYNBB.fit = SYNBB.fit
		UnTransSYNBB.fit$parTable[3:4,1:4]=exp(SYNBB.fit$parTable[3:4,1:4])
		rownames(UnTransSYNBB.fit$parTable)[3]="sdx"
		rownames(UnTransSYNBB.fit$parTable)[4]="sxy"
#print(ma)
		#Write output probability file to working directory
		for (i in 1:length(CurrentAList)){
			extentfile = names(AvailFileNames[i])
			outputfile = paste(synbb.outfile,"-Prob_", origfilename,"_BVN_Model",k,"_", extentfile,sep = 	"")
			write.table (UnTransSYNBB.fit[[7]][[i]], file = outputfile, col.names = TRUE, 	row.names = F, sep = ",")
		}

		#Write output table to working directory
		outputfile = paste(synbb.outfile,"-BVN_Model",k,"_Out.txt",sep = "")
		tmp.wid = getOption("width")  # save current width
		options(width=10000)
		sink(outputfile)              # redirect output to file
		print(UnTransSYNBB.fit[1:6])    # print the object
		sink()                        # cancel redirection
		options(width=tmp.wid)        # restore linewidth

	} #end loop through candidate models list

	#===================================================================================
	#===================================================================================
	# Loop through candidate models; Exponential power null model
	PrevExpParamEsts = array(0,(ncol(AvailList[[1]])+2))
	names(PrevExpParamEsts) = c("mu.x", "mu.y", "ln.scale", "ln.shape", colnames(AvailList[[1]])[-c(1:2)])

	for (k in 1:length(ModelsList)){
    	 #delete columns (i.e., variables) in Availability grids, and Track not used 
      	 Nvariables = sum(ModelsList[[k]]) #number of variables in current model
	 CurrentAList=list()
      	 CurrentTrack=Track[,1:3] #keep x, y and time
	 for (i in 1:length(AvailList)){
	  CurrentAList[[i]]=AvailList[[i]][,1:2] #keep x and y
	 } #end availability list loop
      	 names(CurrentAList)=names(AvailFileNames)     	
	 cc=3
	 for (col in 3:ncol(AvailList[[1]])){
	  if (ModelsList[[k]][col-2]==1){
          CurrentTrack=cbind(CurrentTrack, Track[,col+1])
          colnames(CurrentTrack)[cc+1]=colnames(AvailList[[1]])[col]
	    for (i in 1:length(AvailList)){
	       CurrentAList[[i]]=cbind(CurrentAList[[i]],AvailList[[i]][,col])
             colnames(CurrentAList[[i]])[cc]=colnames(AvailList[[1]])[col]	
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
 names(ThetaW) = colnames(CurrentAList[[1]])[3:ncol(CurrentAList[[1]])]
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

SEP.fit = seple(CurrentTrack,CurrentAList,locAvailFile,start.val=paramSEP)

PrevExpParamEsts[rownames(SEP.fit$parTable)]=SEP.fit$parTable[rownames(SEP.fit$parTable),1]

#Transform back parameter estimates for shape and scale parameters
 UnTransSEP.fit = SEP.fit
 UnTransSEP.fit$parTable[3:4,1:4]=exp(SEP.fit$parTable[3:4,1:4])
 rownames(UnTransSEP.fit$parTable)[3]="scale"
 rownames(UnTransSEP.fit$parTable)[4]="shape"


#Write output probability file to working directory
for (i in 1:length(CurrentAList)){
	extentfile = names(AvailFileNames[i])
	outputfile = paste(synbb.outfile,"-Prob_", origfilename,"_ExpPower_Model",k,"_", extentfile,sep = 	"")
	write.table (UnTransSEP.fit[[7]][[i]], file = outputfile, col.names = TRUE, 	row.names = F, sep = ",")
}

#Write output table to working directory
outputfile = paste(origfilename,"_ExpPower_Model",k,"_Out.txt",sep = "")
tmp.wid = getOption("width")  # save current width
options(width=10000)
sink(outputfile)              # redirect output to file
print(UnTransSEP.fit[1:6])    # print the object
sink()                        # cancel redirection
options(width=tmp.wid)        # restore linewidth

} #end loop through candidate models list

}#end syn function
