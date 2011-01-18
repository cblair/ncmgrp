# Synoptic Model for Analyzing Animal Space Use
# By: Jon Horne; jhorne@uidaho.edu

syn <- function(al) {
	al$ExtentFile <- NULL #blow away ExtentFile col, not needed anymore
	Track = as.matrix(al)
	#===================================================================================
	#===================================================================================
	# Loop through candidate models;

	#Previous SYNBB Paramater Estimates has colnames of the Covariants 
	# (non-required columns)
	RequiredVars <- c("x","y","time","sd")
	CovarColNames <- c(colnames(al[!colnames(al) %in% RequiredVars]))
	PrevSYNBBParamEsts = array(0,(length(CovarColNames))) 
	names(PrevSYNBBParamEsts) = CovarColNames
	
	#for synbb mode:
	for (k in 1:length(ModelsList)){
		print(ModelsList[[k]])
    		#delete columns (i.e., variables) in Track not used 
      		CurrentTrack=Track[,unique(c(RequiredVars,ModelsList[[k]]))] #keep x, y, time, sd
		
		print(CurrentTrack)
  		
	print("TS")
		#-----------------------------------------------------------------------------------
		### Synoptic with brownian bridge
		# Get initial parameter values
		print(ncol(CurrentTrack))
		ThetaW = c(rep(0, (ncol(CurrentTrack)-4)))	#Initial RSF coeff. set to 0; no selection 
		if (k==1){
			#get initial bb standard
			bbsd = sqrt(bb.var)
			print("TSxxx")
		} else { #==> use estimated parameters of previous models (if they exist) for initial values
			#get new estimate of bb.var
			bbsd = exp(SYNBB.fit$partable[1,1])
			names(ThetaW) = colnames(Track)
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

		SYNBB.fit = synbbfit(CurrentTrack,start.val=paramSYNBB, k=k)
		#change sbvnle to synbble, names only

		#PrevBVNParamEsts[rownames(SBVN.fit$parTable)]=SBVN.fit$parTable[rownames(SBVN.fit$parTable),1]
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
