source("bb.r")
source("syn.r")
source("gen.r")
require(survival)
require(maptools)
require(sp)

################################################
#Handle Options
i <- 0
cluster <- ""
bb.only <- FALSE
syn.only <- FALSE
parallel <- FALSE
do.to.js <- FALSE #date offset to julian seconds
nodes <- 1
options <- commandArgs()
for(option in options) {
	i <- i + 1
	print(options[i])
    	if(option == "-c"){
       		cluster <- options[i + 1]
     	}
	if(option == "--bb-only") {
		bb.only <- TRUE
	}
	if(option == "--syn-only") {
		syn.only <- TRUE
	}
	if(option == "--do-to-js") {
		do.to.js <- TRUE
	}
	if(option == "-p") {
		parallel <- TRUE
		nodes <- options[i + 1]
	}
}

###############################################
#Parallel Setup
if(parallel == TRUE) {
	require(Rmpi)
	require(snow)
	c1 <- makeCluster(nodes, type="MPI")
}
source("step.r")

#################################################
#Set working directories
workdir = paste(getwd(),"/data/",cluster,sep="")
basedir = getwd()
print(paste("Using",workdir,"as workdir"))
setwd(workdir)

#################################################
#Profiling for performance
profile.fname = paste(workdir,"/profile.dat", sep="")
Rprof(profile.fname)

#################################################
#Setup data
global.Map.g.a <- matrix()
alfname = paste(workdir,"/",cluster,"_all_locations.txt",sep="")
#al = all locations
al = read.table(alfname, sep='\t',header=TRUE, as.is=TRUE)
#ma setup
ma.flist = unique(al$ExtentFile)

#ma.flist = list.files(workdir,pattern=paste(cluster,"_master_avail*",sep=""))
ma.list <- lapply(1:length(ma.flist), function(i) {
	masterfname = paste(workdir,"/",ma.flist[i],sep="")
	retval <- read.table(masterfname, sep='\t',header=TRUE, as.is=TRUE)
	retval$ExtentFile = ma.flist[i]
	return(retval)
} )

#################################################
#Process for each ma in the ma list
 
#make one big ma for use in normalization
ma.combined <- ma.list[[1]]
if(length(ma.list) > 1) {
	for(i in 2:length(ma.list)) {ma.combined <- rbind(ma.combined, ma.list[[i]])}
}


#al normalize
for(coln in names(al)) {
	if((coln != "x") && (coln != "y") && (coln != "time") && (coln != "sd") && (coln != "ExtentFile")) {
		#al[,i] <- (al[,i] - min(ma[,i-2])) /  (max(ma[,i-2]) - min(ma[,i-2]))
		if(coln %in% names(ma.combined)) {
			print(paste("Normalizing location data column",coln))
			al[,coln] <- (al[,coln] - min(ma.combined[,coln])) /  (max(ma.combined[,coln]) - min(ma.combined[,coln]))
		} else {
			print(paste("Warning: column",coln,"in location data, but not in habitat data."))
		}
	}
}
#end al normalize
#to add - tell them mins and max's and ask them if they want to input their own
#ma normalize
for(i in 1:length(ma.list)) {
	for(coln in names(ma.combined)) {
		if((coln != "x") && (coln != "y") && (coln != "time") && (coln != "sd") && (coln != "ExtentFile")) {
               		#ma[,i] <-  (ma[,i] - min(ma[,i])) /  (max(ma[,i]) - min(ma[,i])) 
			print(paste("Normalizing habitat data column",coln))
			ma.list[[i]][,coln] <-  (ma.list[[i]][,coln] - min(ma.combined[,coln])) /  (max(ma.combined[,coln]) - min(ma.combined[,coln])) 
		}
	}
}
#end ma normalize

#record what mins and maxes we used for normlization
#x <- data.frame(min=(min(ma.combined)),max=(max(ma.combined)))
#write.table(x,file="norm_data.r")

#import habitat (ma) models
source("models.r")
#add required values to models
#ModelsList <- lapply(1:length(ModelsList), function(z) {
#		return(c("x","y",ModelsList[z],recursive=TRUE))
#	} )

#build cellgrid
# calculate auxillary variables
bb.var <- get.bb.var(al) #bb variance
cellgrid <- list()
print("Starting cellgrid construction...")
starttime <- proc.time()[3]
if(parallel) {	
	clusterExport(c1, "al")
	clusterExport(c1, "ma.list")
	#cellgrid <- parLapply(c1, 1:length(al$x), get.cellgrid)
	print("Starting par cellgrid construction")
	cellgrid <- parLapply(c1, 1:length(al$x), get.cellgrid.with.mas)
} else {
	#cellgrid <- lapply(1:length(al$x), get.cellgrid)
	cellgrid <- lapply(1:length(al$x), get.cellgrid.with.mas)
}
cellgrid <- cellgrid[!is.na(cellgrid)]
endtime <- proc.time()[3]
print(paste("Cellgrid construction time was",endtime - starttime))


ma.gridsize <- get.ma.gridsize(ma.combined) 


#################################################
#Process options
if(do.to.js == TRUE) {
	al$Julian <- al$time * 24 * 60 * 60
	al = al[order(al$Julian),]
}

#else, we are doing bb and syn by triplicates
synbb.outfile = "step"

        
#################################################3
#Start for processing

#apply operations foreach element in the cellgrid
foreach.cellgrid <- function(i) {
	if(!is.na(cellgrid[i][[1]][[2]])) {
		synbb.outfile <<- paste("trip-",i,sep="")
		#print(paste("Running bb for cellgrid",i,"of",length(cellgrid)))
		#cellgrid[i][[1]][[3]]$Julian <- cellgrid[3][[1]][[3]]$time * 24 * 60
		#bb(cellgrid[i][[1]][[3]], bb.var)
	}
	#else {
		
	#}
}

#if(parallel) {	
#	clusterExport(c1, "cellgrid")
#	clusterExport(c1, "bb")
#	parLapply(c1, 1:length(cellgrid), foreach.cellgrid)
#} else {
	lapply(1:length(cellgrid), foreach.cellgrid)
#}


print("Processing syn using cellgrids...")
synbb.outfile = "syn"
syn(al)


#close parallel ops
if(parallel == TRUE) {
	stopCluster(c1)
}

#stop and summarize profiling
summaryRprof(profile.fname)
