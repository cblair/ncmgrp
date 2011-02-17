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
		starttime <- proc.time()[3]
		local.outlier <- FALSE
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


		#!
		bbsd = sqrt(bb.var) 
		StartX1 <- al$x[i - 1]
                StartY1 <- al$y[i-1]
                EndX3 <- al$x[i + 1]
                EndY3 <- al$y[i + 1] 
                StartSTD1 <- al$sd[i - 1] 
                EndSTD3 <- al$sd[i + 1]
                StartTime1 <- al$time[i - 1]
                EndTime3 <- al$time[i + 1] 
                TotalTime13 <- EndTime3 - StartTime1
		Time2 <- al$time[i] - al$time[i-1]

		Alpha <- Time2 / TotalTime13
		MeanXTime2 <- StartX1 + ((Alpha) * (EndX3 - StartX1))
                MeanYTime2 <- StartY1 + ((Alpha) * (EndY3 - StartY1))
		VarTime2 <- TotalTime13 * Alpha * (1 - Alpha) * (bbsd ^ 2) + ((1 - Alpha) ^ 2) * (StartSTD1 ^ 2) + (Alpha ^ 2) * (EndSTD3 ^ 2)
		SDTime2 <- sqrt(VarTime2)
		CI99 <- 3 * SDTime2
		#^here, in the future, may be an option of changing 3
		

		newxmin <- MeanXTime2 - CI99
		newxmax <- MeanXTime2 + CI99
		newymin <- MeanYTime2 - CI99
		newymax <- MeanYTime2 + CI99

		#check if pos 2 is out of extent
		if((al$x[i] > newxmax) || (al$x[i] < newxmin) 
		|| (al$y[i] > newymax) || (al$y[i] < newymin)
		) {
			print(paste("Middle location",i," is outside of extent"))
			local.outlier <- TRUE
		}


		#get ma values from the corresponding ms.list element
		# *build the ma based on the middle triplicate's ExtentFile
		ma.name <- ""
		ma.element <- data.frame()
		for(ma in ma.list) {
			if(ma$ExtentFile[1] == al$ExtentFile[i]) {
				ma.element <- ma
				ma.name <- ma$ExtentFile[1]
			}
		}

		#get ma values within x and y values of this cellgrid cell
		ma.clip <- NA
		ma.clip <- ma.element[(ma.element$x >= newxmin),]
		ma.clip <- ma.clip[(ma.clip$x <= newxmax),]
		ma.clip <- ma.clip[(ma.clip$y >= newymin),]
		ma.clip <- ma.clip[(ma.clip$y <= newymax),]

	
		#get ma cellsize
		ma.cellsize = max(abs(diff(ma.clip$x))) * max(abs(diff(ma.clip$y)))
		#^fix this!

		#create ma of Models
		ma.melement <- data.frame()
		#add x and y cordinates here
		ModelsListTemp <- lapply(1:length(ModelsList), function(z) {
               		return(c("x","y",ModelsList[z],recursive=TRUE))
		       } )
		ma.mlist <- lapply(1:length(ModelsListTemp), function(z) {
        		#this doesn't work like you think it does...
			for(name in ModelsListTemp[z]) {
				ma.melement <<- ma.clip[name]
        		}
        		return(ma.melement)
		} ) 

		endtime <- proc.time()[3]
		runtime <- endtime - starttime
		#print(paste("cellgrid",i,"done bin'ing in ",runtime,"seconds"))
		#cellgrid structure:
		#	cellgrid[x][[1]] = triplicate properties
		#			- $ma.name
		#			- $ma.cellsize, avg area of ma cells
		#			- $outlier - if 2nd location was out of extent
		#	cellgrid[x][[2]] = this cellgrid cell's ma values 
		#	cellgrid[x][[3]] = triplicate info dataframe from the location dataframe
		return(list(data.frame(ma.name=ma.name, ma.cellsize=ma.cellsize,outlier=local.outlier), ma.mlist, as.data.frame(rbind(al[i - 1,], al[i,], al[i + 1,]))))
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
