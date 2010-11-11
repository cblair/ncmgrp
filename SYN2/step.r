get.cellgrid <- function(i) {
		if(i %% 2 == 0) {
			xsection <- c(al$x[i-1], al$x[i], al$x[i+1])
			ysection <- c(al$y[i-1], al$y[i], al$y[i+1])
			#if one of the trips does not exist, return
			if(is.na(al[i+1,]) || is.na(al[i+1,]) || is.na(al[i+1,])) {
				return(NA)
			}
			xmin <- min(xsection)
			xmax <- max(xsection)
			ymin <- min(ysection)
			ymax <- max(ysection)
			newxmin <- xmin - (xmax - xmin)
			newxmax <- xmax + (xmax - xmin)
			newymin <- ymin - (ymax - ymin)
			newymax <- ymax + (ymax - ymin)
			return(list(data.frame(xmin=newxmin, xmax=newxmax, ymin=newymin, ymax=newymax), 0, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
		}
		return(NA)
}#end function

#experiment for Jon
get.cellgrid.with.mas <- function(i) {
	if(i %% 2 == 0) {
		xsection <- c(al$x[i-1], al$x[i], al$x[i+1])
		ysection <- c(al$y[i-1], al$y[i], al$y[i+1])
		#if one of the trips does not exist, return
		if(is.na(al[i+1,]) || is.na(al[i+1,]) || is.na(al[i+1,])) {
			print(paste("cellgrid",i,"skipped"))
			return(NA)
		}
		xmin <- min(xsection)
		xmax <- max(xsection)
		ymin <- min(ysection)
		ymax <- max(ysection)
		newxmin <- xmin - (xmax - xmin)
		newxmax <- xmax + (xmax - xmin)
		newymin <- ymin - (ymax - ymin)
		newymax <- ymax + (ymax - ymin)

		#get ma values within x and y values of this cellgrid cell
		ma.temp <- NA
		ma.temp <- ma[(ma$x >= newxmin),]
		ma.temp <- ma.temp[(ma.temp$x <= newxmax),]
		ma.temp <- ma.temp[(ma.temp$y >= newymin),]
		ma.temp <- ma.temp[(ma.temp$y <= newymax),]
		
		print(paste("cellgrid",i,"done bin'ing"))
		#cellgrid structure:
		#	cellgrid[x][[1]] = filterred ma, used later for different ma models, reduction
		#			of columns on cellgrid[[2]], i.e different models
		#	cellgrid[x][[2]] = this cellgrid cell's ma values 
		#	cellgrid[x][[3]] = triplicate info dataframe from the location dataframe
		return(list(data.frame(), ma.temp, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
	}
	return(NA)
}#end function


step <- function() {
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
		print(paste("loc",z,"done"))
 	 }
	)

	#write each ma to file
	#for(row in cellgrid) {
	#	if(!is.na(row[[2]])) {
	#		write.table(row[[2]], file=paste(workdir,"/",synbb.outfile,"-ma-",row[[1]]$xmin,"-",row[[1]]$xmax,"-",row[[1]]$ymin,"-",row[[1]]$ymax,"_ma.txt",sep=""), append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
	#	}
	#}

return(cellgrid)
}#end step function
