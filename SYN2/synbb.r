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
alfname = paste(workdir,"/",cluster,"_all_locations.txt",sep="")
#al = all locations
al = read.table(alfname, sep='\t',header=TRUE, as.is=TRUE)
masterfname = paste(workdir,"/",cluster,"_master_avail.txt",sep="")
ma = read.table(masterfname, sep='\t',header=TRUE, as.is=TRUE)
#normalize
for(i in 3:ncol(ma)) { #start at 3, ignore x and y cols
                ma[,i] <- ma[,i] / max(ma[,i])
        }
#end normalize

#################################################
#Process
if(do.to.js == TRUE) {
	al$Julian <- al$time * 24 * 60 * 60
	al = al[order(al$Julian),]
}
if(bb.only == TRUE) {
	synbb.outfile <- cluster
	print("Processing bb without cellgrids...")
	bb(al)
	q()
}
if(syn.only == TRUE) {
	synbb.outfile <- cluster
	print("Processing syn without cellgrids...")
	syn(al,ma)
	q()
}

#else, we are doing bb and syn by triplicates
synbb.outfile = "step"

        #################################################3
        #Start for processing
        if(parallel) {
                clusterExport(c1, "al")
                cellgrid <- parLapply(c1, 1:length(al$x), get.cellgrid)
                #cellgrid <- lapply(1:length(al$x), get.cellgrid)

        } else {
                cellgrid <- lapply(1:length(al$x), get.cellgrid)
        }

cellgrid = step()

bb.var <- get.bb.var()

foreach.cellgrid <- function(i) {
	if(!is.na(cellgrid[i][[1]][[2]])) {
		synbb.outfile <<- paste("trip-",i,sep="")
		print(paste("Running bb for cellgrid",i,"of",length(cellgrid)))
		cellgrid[i][[1]][[3]]$Julian <- cellgrid[3][[1]][[3]]$time * 24 * 60
		bb(cellgrid[i][[1]][[3]], bb.var)
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
syn(al, ma)

#close parallel ops
if(parallel == TRUE) {
	stopCluster(c1)
}

#stop and summarize profiling
summaryRprof(profile.fname)
