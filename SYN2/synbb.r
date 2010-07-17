source("step.r")
source("bb.r")

################################################
#Handle Options
i <- 0
cluster <- ""
bb.only <- FALSE
jh.to.js <- FALSE
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
	if(option == "--jh-to-js") {
		jh.to.js <- TRUE
	}
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
if(jh.to.js == TRUE) {
	al$Julian <- al$Time * 60 * 60
	al = al[order(al$Julian),]
	al
}
if(bb.only == TRUE) {
	print("Need to test with original bb data.")
	bb(al)
	q()
}


cellgrid = step()

cellgrid[3][[1]][[3]]
bb(cellgrid[3][[1]][[3]])
