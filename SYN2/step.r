library(plyr)

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
setwd(workdir)

alfname = paste(workdir,"/",cluster,"_all_locations.txt",sep="")
#al = all locations
al = read.table(alfname, sep='\t',header=TRUE, as.is=TRUE)
masterfname = paste(workdir,"/",cluster,"_master_avail.txt",sep="")
ma = read.table(masterfname, sep='\t',header=TRUE, as.is=TRUE)

##get the potential column size for future ma matrix based off of the potential y column
ma.getpotcolsize <- function() {
	i <- 1
	while(ma$y[i] == ma$y[1]) {
		i <- i + 1
	}
	return(i)
}

#################################################3
#Start for processing

cellgrid <- lapply(1:length(al$x), function(i) {
	if(i %% 2 == 0) {
		#print("Processing...")
		#print(i)
		xsection <- c(al$x[i-1], al$x[i], al$x[i+1])
		ysection <- c(al$y[i-1], al$y[i], al$y[i+1])
		xmin <- min(xsection)
		xmax <- max(xsection)
		ymin <- min(ysection)
		ymax <- max(ysection)
		newxmin <- xmin - (xmax - xmin)
		newxmax <- xmax + (xmax - xmin)
		newymin <- ymin - (ymax - ymin)
		newymax <- ymax + (ymax - ymin)
		#return(c(newxmin, newxmax, newymin, newymax, matrix(1, 1)))
		return(c(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0))
	}
	return(NA)
 }#end function
)#end lapply
#clean out NA values
cellgrid <- cellgrid[!is.na(cellgrid)]

#reset the appropriate cellgrid vector element to an ma variable
for(i in 1:length(cellgrid)) {
	#print(cellgrid[i][[1]][[5]])
	cellgrid[i][[1]][[5]] <- ma[1,]
	for(j in 1:ncol(cellgrid[i][[1]][[5]])) {
		cellgrid[i][[1]][[5]][,j] <- NA
	}
	#print(cellgrid[i][[1]][[5]])
}

#test <- ma[1,]
#test[1][[1]] <- NA
#test[2][[1]] <- NA
#test[3][[1]] <- NA
#test[4][[1]] <- NA
#test[5][[1]] <- NA
#testi <- function () {
#test2 <- data.frame(x=1, y=2, d2et=3, slope=4, elevm=5)
#test3 <- data.frame(x=6, y=7, d2et=8, slope=9, elevm=10)
#cellgrid[1][[1]][[5]] <- rbind(cellgrid[1][[1]][[5]], test2)
#cellgrid[1][[1]][[5]] <- rbind(cellgrid[1][[1]][[5]], test3)
#}
#testi()
#test <- rbind(test, test3)
#test
#cellgrid[1][[1]][[5]]
#q()

al.addToLocCells <- function(row) {
	for(j in 1:length(cellgrid)) {
		#print(cellgrid[j][[1]])
		#print(row$x)
		#print((row$x <= cellgrid[j][[1]]$xmax))
		if(((row$x >= cellgrid[j][[1]]$xmin) && (row$x <= cellgrid[j][[1]]$xmax) && (row$y >= cellgrid[j][[1]]$ymin) && (row$y <= cellgrid[j][[1]]$ymax))) {
			#print("inserting into cell:")
			#print(j)
			#print(cellgrid[j][[1]][[5]])
			cellgrid[j][[1]][[5]] <<- rbind(cellgrid[j][[1]][[5]], row)
		}
	}
}

#add each row of ma to the appropriate cell in the cell grid
z <- 1
by(ma, 1:nrow(ma), function(row) {
	starttime <- proc.time()[2]
	al.addToLocCells(row)
	runtime <- proc.time()[2] - starttime
	#print("cell grip proc time:")
	#print(runtime)
	z <<- z + 1
	print("loc done")
	print(z)
 }
)

cellgrid