step <- function() {
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
			return(c(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
			#return(list(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
		}
		return(NA)
	 }#end function
	)#end lapply
	#clean out NA values
	#print(cellgrid)
	#q()
	cellgrid <- cellgrid[!is.na(cellgrid)]

	#reset the appropriate cellgrid vector element to an ma variable
	for(i in 1:length(cellgrid)) {
		cellgrid[i][[1]][[5]] <- ma[1,]
		for(j in 1:ncol(cellgrid[i][[1]][[5]])) {
			cellgrid[i][[1]][[5]][,j] <- NA
			#cellgrid[i][[1]][[5]][,j] <- character(0)
		}
	}

	al.addToLocCells <- function(row) {
		for(j in 1:length(cellgrid)) {
			if(((row$x >= cellgrid[j][[1]]$xmin) && (row$x <= cellgrid[j][[1]]$xmax) && (row$y >= cellgrid[j][[1]]$ymin) && (row$y <= cellgrid[j][[1]]$ymax))) {
				if(is.na(cellgrid[j][[1]][[5]])) {
					cellgrid[j][[1]][[5]] <<- row
				}
				else {
					cellgrid[j][[1]][[5]] <<- rbind(cellgrid[j][[1]][[5]], row)
				}
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

	#write each ma to file
	for(row in cellgrid) {
		if(!is.na(row[[5]])) {
			write.table(row[[5]], file=paste(workdir,"/",row$xmin,"-",row$xmax,"-",row$ymin,"-",row$ymax,"_ma.txt",sep=""), append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
		}
	}

return(cellgrid)
}#end step function
