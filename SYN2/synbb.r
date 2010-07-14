source("step.r")
source("bb.r")

i <- 0
cluster <- ""
options <- commandArgs()
for(option in options) {
           i <- i + 1
           if(option == "-i"){
           cluster <- options[i + 1]
           }
}


#################################################3
#Setup for processing

workdir = paste(getwd(),"/data/",cluster,sep="")
basedir = getwd()
setwd(workdir)

cellgrid = step()

locs <- c()
cellgrid[3]
for(i in 1:ncol(cellgrid[3][[1]])) {
	print(i)
}
#"rbind:"
#rbind(cellgrid[3]$trip1, cellgrid[3]$trip2, cellgrid[3]$trip3)

#bb(

