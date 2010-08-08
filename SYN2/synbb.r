source("step.r")
source("bb.r")
source("syn.r")
require(survival)
require(maptools)
require(sp)
require(Rmpi)

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
	require(snow)
	c1 <- makeCluster(nodes, type="MPI")
}

#################################################
#Set working directories
workdir = paste(getwd(),"/data/",cluster,sep="")
basedir = getwd()
setwd(workdir)

#################################################
#Setup for processing
alfname = paste(workdir,"/",cluster,"_all_locations.txt",sep="")
#al = all locations
al = read.table(alfname, sep='\t',header=TRUE, as.is=TRUE)
masterfname = paste(workdir,"/",cluster,"_master_avail.txt",sep="")
ma = read.table(masterfname, sep='\t',header=TRUE, as.is=TRUE)

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
cellgrid = step()

lapply(1:length(cellgrid), function(i) {
	if(!is.na(cellgrid[i][[1]][[2]])) {
		synbb.outfile <<- paste("trip-",i,sep="")
		print(paste("Running bb for cellgrid",i,"of",length(cellgrid)))
		cellgrid[i][[1]][[3]]$Julian <- cellgrid[3][[1]][[3]]$time * 24 * 60
		bb(cellgrid[i][[1]][[3]])
	}
	else {
		
	}
 }
)

print("Processing syn using cellgrids...")
synbb.outfile = "syn"
syn(al, ma)

if(parallel == TRUE) {
	stopCluster(c1)
}
