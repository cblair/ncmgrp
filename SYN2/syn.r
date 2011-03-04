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
	SYNBB.fit <- data.frame()
	#for synbb mode:
	for (k in 1:length(ModelsList)){
		#print(ModelsList[[k]])
    		#delete columns (i.e., variables) in Track not used 
      		CurrentTrack=Track[,unique(c(RequiredVars,ModelsList[[k]]))] #keep x, y, time, sd
		
		#print(CurrentTrack)
  		
		#-----------------------------------------------------------------------------------
		### Synoptic with brownian bridge
		# Get initial parameter values
		#print(ncol(CurrentTrack))
		ThetaW = c(rep(0, (ncol(CurrentTrack)-4)))	#Initial RSF coeff. set to 0; no selection 
		print("TS31")
		print(ThetaW)
		if (k==1){
			#get initial bb standard
			bbsd = sqrt(bb.var)
			#print("TSxxx")
		} else { #==> use estimated parameters of previous models (if they exist) for initial values
			#get new estimate of bb.var
			print("TS42")
			#ThetaW <- array(0,(length(colnames(CurrentTrack))))
			ThetaW <- c()
			for(name in colnames(CurrentTrack)) {
				print(paste("TS44:",PrevSYNBBParamEsts[name]))
				if(!is.na(PrevSYNBBParamEsts[name])) {
					ThetaW <- c(ThetaW, PrevSYNBBParamEsts[name][[1]])
				}
			}
			#ThetaW <- c(ThetaW[1:length(ThetaW)])
			print("TS37:")
			print(head(CurrentTrack))
			print(ThetaW)
			bbsd = exp(SYNBB.fit$partable[1,1])
			#names(ThetaW) = colnames(Track)
			#names(ThetaW) = colnames(Track)
			#names(ThetaW) = colnames(CurrentTrack[,1:2])			
		} #end if

		print("TS58")
		print(ThetaW)

		lnbbsd = log(bbsd)
		paramSYNBB = c(lnbbsd, ThetaW)
		
		print("TS61: start.val")
		print(paramSYNBB)
		SYNBB.fit <- synbbfit(CurrentTrack,start.val=paramSYNBB, k=k)
		#change sbvnle to synbble, names only

		print("TS57: SYNBB.fit")
		print(SYNBB.fit[1])

		#PrevBVNParamEsts[rownames(SBVN.fit$parTable)]=SBVN.fit$parTable[rownames(SBVN.fit$parTable),1]
		PrevSYNBBParamEsts[rownames(SYNBB.fit$parTable)]=SYNBB.fit$parTable[rownames(SYNBB.fit$parTable),1]

		#Transform back parameter estimates for sdbb
		print("TS169")
		#UnTransSYNBB.fit = SYNBB.fit
		#UnTransSYNBB.fit$parTable[1,1:4]=exp(SYNBB.fit$parTable[1,1:4])
		#UnTransSYNBB.fit$parTable[1,2:5]=exp(SYNBB.fit$parTable[1,2:5])
		#rownames(UnTransSYNBB.fit$parTable)[1]="sdbb"
		print("TS169")
		print(SYNBB.fit$partable)
		#UnTransSYNBB.fit = SYNBB.fit$parTable[1,1:4]
		print("TSL67")
		#print(UnTransSYNBB.fit)
		print("wtf")
		UnTransSYNBB.fit=exp(SYNBB.fit$partable)
		print("TSL71")
		print(UnTransSYNBB.fit)
		print("TSL73")
		#rowname(UnTransSYNBB.fit[1])="sdbb"
		print(UnTransSYNBB.fit)
		#Write output probability file to working directory
		'
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
		'
	} #end loop through candidate models list

}#end syn function
