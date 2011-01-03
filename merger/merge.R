gps.act.merger <- data.frame()

get.gps.act.merger.list <- function(x,gps.act.merger) {
			col <- colnames(gps.act.merger)[x]
                        #print(paste("Processing column",col))
                        #setup local data.frame as the current column
                        local.gps.act.merger <- data.frame(gps.act.merger[col])

                        stime <- proc.time()[3]
                        for(i in 1:length(local.gps.act.merger[,col])) {
				if(is.na(local.gps.act.merger[i,col]) && i > 1) {
					print("TS")
                                        local.gps.act.merger[i,col] <- local.gps.act.merger[i - 1,col]
                                }
                        }
                        etime <- proc.time()[3]
			#print(paste("Replace Na's in row",i,"took",etime - stime))
return(local.gps.act.merger)
}


patt <- "Collar(.*)_merger.txt"
flist <- list.files(path="gps_act_merge/",pattern=patt,)

i <- 0
parallel <- FALSE
nodes <- 0
options <- commandArgs()
for(option in options) {
        i <- i + 1
        print(options[i])
        if(option == "-i"){
                flist <- c(options[i + 1])
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


print(paste("Processing",length(flist),"collar files"))
lapply(flist[1:length(flist)], function(file) {
	collarnum = unlist(strsplit(gsub(patt,"\\1",file,perl=TRUE),","))
	collar <- paste("Collar",collarnum,sep="")
	print(paste("Processing collar '",collar,"'",sep=""))

	print("Opening merger file...")
	merger <- read.table(file=paste("gps_act_merge/",collar,"_merger.txt",sep=""), header=T, row.names=NULL)
	print("Opening activity file...")
	act <- read.table(file=paste("gps_act_merge/",collar,"_Activity_data.txt",sep=""), header=T, row.names=NULL)
	print("Opening gps file...")
	gps <- read.table(file=paste("gps_act_merge/",collar,"_GPS_data.txt",sep=""), header=T, row.names=NULL)

	print("Merging gps and act...")
	gps.act <- merge(act,gps,by.y="LMT",all.x=T,incomparables=NA)
	print("Merging gps act and merger...")
	gps.act.merger <- merge(gps.act,merger,by.y="FixNum",all.x=T,incomparables=NA)

	if(dim(gps.act.merger) != dim(act)) {
		print(collar)
		print("Dims off...")
	}

	#gps.act.merger <- as.data.frame(gps.act.merger)
	#set nas to previous
	print("Changing NA's to previous values...") 
	print(paste("est. run time:",0.00829662146004851 * length(gps.act.merger[,1])))
	stime <- proc.time()[3]
	#map
	if(parallel) {
		print("Changing NA's in parallel")
		gps.act.merger.list <- parLapply(c1, 1:ncol(gps.act.merger), get.gps.act.merger.list, gps.act.merger=gps.act.merger)
	} else {
		gps.act.merger.list <- lapply(1:ncol(gps.act.merger), get.gps.act.merger.list, gps.act.merger=gps.act.merger)
	}
	etime <- proc.time()[3]
	ttime <- etime - stime
	pertime <- ttime / length(gps.act.merger[,1])
	print(paste("done replacing NA's in",ttime,",",pertime,'per gps.act.merger entry, entries:',length(gps.act.merger[,1])))
		
	#reduce
	gps.act.merger <- data.frame(gps.act.merger.list)
	print(paste("gps_act_merge-merged/",collar,".txt"))
	write.table(gps.act.merger, file=paste("gps_act_merge-merged/",collar,".txt",sep=""),sep="\t")
})


#close parallel ops
if(parallel == TRUE) {
	stopCluster(c1)
}
