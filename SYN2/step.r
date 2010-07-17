step <- function() {
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
			#return(c(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
			return(list(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
		}
		return(NA)
	 }#end function
	)#end lapply
	#clean out NA values
	cellgrid <- cellgrid[!is.na(cellgrid)]
	#reset the appropriate cellgrid vector element to an ma variable
	for(i in 1:length(cellgrid)) {
		cellgrid[i][[1]][[2]] <- ma[1,]
		for(j in 1:ncol(cellgrid[i][[1]][[2]])) {
			cellgrid[i][[1]][[2]][,j] <- NA
		}
	}

	al.addToLocCells <- function(row) {
		for(j in 1:length(cellgrid)) {
			if(((row$x >= cellgrid[j][[1]][[1]]$xmin) && (row$x <= cellgrid[j][[1]][[1]]$xmax) && (row$y >= cellgrid[j][[1]][[1]]$ymin) && (row$y <= cellgrid[j][[1]][[1]]$ymax))) {
				if(is.na(cellgrid[j][[1]][[2]])) {
					cellgrid[j][[1]][[2]] <<- row
				}
				else {
					cellgrid[j][[1]][[2]] <<- rbind(cellgrid[j][[1]][[2]], row)
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
		if(!is.na(row[[2]])) {
		print(row[[1]])
			write.table(row[[2]], file=paste(workdir,"/",row[[1]]$xmin,"-",row[[1]]$xmax,"-",row[[1]]$ymin,"-",row[[1]]$ymax,"_ma.txt",sep=""), append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
		}
	}

return(cellgrid)
}#end step function
