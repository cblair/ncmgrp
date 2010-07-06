require(survival)
require(maptools)
require(sp)
require(snow)

##get input file
ifile <- ""
i <- 0
options <- commandArgs()
for(option in options) {
	   i <- i + 1
	   if(option == "-i"){
	   ifile <- options[i + 1]
	   }
}
##end get input file

print("Getting from")
print(ifile)

# Read in GPS locations of migraion
#locs = read.table("data/gps53_spring06.txt", sep='\t', 
#        header=TRUE, as.is=TRUE)

locs = read.table(ifile, sep='\t', 
        header=TRUE, as.is=TRUE)

# Convert time to Julian 
locs$Date =  as.Date(locs$Date, "%m/%d/%Y")
tmp = unlist(strsplit(locs$Hour, ":"))
tmp2 = rep(c(1:2), nrow(locs)/2)
Hour = as.numeric(tmp[tmp2 == 1])
Minute = as.numeric(tmp[tmp2 == 2])
locs$Julian = as.numeric(round((julian(locs$Date) - 1)*24*60 + 
                Hour*60 + Minute, 0))
locs = locs[order(locs$Julian),]
    

#-------------------------------------------------------------------
#   Brownian Bridge Function
#
#	INPUT:
#
#	X = vector of x coordinates (meters)
#	Y = vector of y coordinates (meters)
#   Time = vector of location times (julian minutes)
#   LocationError = mean location error
#
#	NOTE: X, Y, Time, and LocationError must be the same 
#       length, and all sorted according to Time (ascending).
#
#	cell.size = cell size for density grid.
#
#-------------------------------------------------------------------

BrownianBridge <- function(X, Y, Time, LocationError, cell.size=50){

	options(digits=20)

    max.lag = abs(min(diff(Time))+1)

    LocationError = LocationError/sqrt(pi/2)

	# Create grid
    system.time(xmin <- (min(X) - 0.1*(max(X)-min(X)))) 
    system.time(xmax <- (max(X) + 0.1*(max(X)-min(X))))
    system.time(ymin <- (min(Y) - 0.1*(max(Y)-min(Y))))
    system.time(ymax <- (max(Y) + 0.1*(max(Y)-min(Y))))
    x = seq(xmin, xmax, by=cell.size)
	y = seq(ymin, ymax, by=cell.size)
	x = rep(x, length(y))
	y = sort(rep(y, length(x)/length(y)))
	grid = data.frame(x, y)
	cat("Size of grid =", length(unique(grid$y)), "X", 
            length(unique(grid$x)), fill=TRUE)

	n.locs = length(X)
    cat("Number of animal locations =", n.locs, fill=TRUE)

    #--------------------------------------------------------------
	#Calculate Brownian motion variance "BMvar"
    LocationError = rep(LocationError, length(X))
    T.jump = alpha = ztz = ob = loc.error.1 = loc.error.2 = NULL
    i = 2
    while(i < n.locs){
    	if((Time[i + 1] - Time[i - 1]) / 2 > max.lag) {i = i + 1}
    	else {
    		ob = c(ob, i)
    		t = Time[i+1]-Time[i-1]
    		T.jump = c(T.jump, t)
    		a = (Time[i]-Time[i-1]) / t
    		alpha = c(alpha, a)
    		u = c(X[i-1], Y[i-1]) + a*(c(X[i+1], Y[i+1]) - c(X[i-1], Y[i-1]))
    		ztz = c(ztz, (c(X[i], Y[i]) - u)%*%(c(X[i], Y[i]) - u))
    		loc.error.1 = c(loc.error.1, LocationError[i-1])
    		loc.error.2 = c(loc.error.2, LocationError[i+1])
    	}
    	i = i + 2
    }

    #Likelihood (Eq.7 in Horne et al. 2007)
    likelihood <- function(var){	
    	v = T.jump*alpha*(1-alpha)*var + 
            ((1-alpha)^2)*(loc.error.1^2) + 
    		(alpha^2)*(loc.error.2^2)
    	#!!this is where ztz Null will kill the prog
    	l = (1/(2*pi*v))*exp(-ztz / (2*v))
    	-sum(log(l), na.rm=T)
    }
    BMvar = nlminb(start=10000, likelihood, lower=10)$par
    cat("Brownian Motion Variance (meters^2) =", round(BMvar), 
        fill=TRUE)
    BMvar = rep(BMvar, length(X))

	#Estimate Brownian bridge
	#Note: 5 minute time step (dt interval in Eq.5 Horne et al. 2007) is used.
	Time.Diff = c(diff(Time), NA)
	T.Total = Time[n.locs] - Time[1]
	probability = NULL
	int = 0
	for(i in 1:(n.locs-1)){
                proctime1 = proc.time()
		if(Time.Diff[i] <= max.lag){
			theta = NULL
			tm = 0
			while(tm <= Time.Diff[i]){
				alpha = tm/Time.Diff[i]
				mu.x = X[i] + alpha*(X[i+1] - X[i])
				mu.y = Y[i] + alpha*(Y[i+1] - Y[i])
				sigma.2 = Time.Diff[i]*alpha*(1-alpha)*BMvar[i] + 
					((1-alpha)^2)*(LocationError[i]^2) + 
					(alpha^2)*(LocationError[i+1]^2)
				ZTZ = (grid$x - mu.x)^2 + (grid$y - mu.y)^2
				theta = (1/sqrt(2*pi*sigma.2))*exp(-ZTZ/(2*sigma.2)) 
				int = int + theta
				tm = tm + 5
			}
                proctime2 <- proc.time()
                runtime <- proctime2[2] - proctime1[2]
                usertime <- proctime2[1] - proctime1[1]
                #cat("Runtime: ")
                #print(runtime)
                #cat("User Time: ")
                #print(usertime)
		}
	}

	# Scale probabilities to sum = 1.0
	probability = int/T.Total
	probability = probability/sum(probability)

	BB = data.frame(x = grid$x, y = grid$y, probability=probability)
	BB = BB[order(BB$x, BB$y),]

	# Create contour map
    prob = sort(BB$probability)
	cumprob = round(cumsum(prob), 4)
	contour99 = max(prob[cumprob <= (1-0.99)])
    #windows(height=6, width=8)
	z = matrix(BB$probability, nrow=length(unique(BB$x)), 
        ncol=length(unique(BB$y)), byrow=T)
    contour(x=unique(BB$x), y=unique(BB$y), z, levels=contour99, drawlabels=F, xlab="X", ylab="Y")
    points(X, Y, pch=19, cex=0.2, col="steelblue")
    lines(X, Y, lwd=0.5, col="steelblue")

  return(BB)

}


# Create output ASCII file with probabilities (scaled so max=1000)
print("Brown Bridge Runtime...")
system.time(density <- BrownianBridge(locs$X, locs$Y, Time=locs$Julian, 
            LocationError=20, cell.size=30)
	    )
bb = density

print("round func runtime...")
system.time(bb$probability <- round(bb$probability*1000/max(bb$probability), 0))

print("data.fram runtime...")
system.time(m <- data.frame(bb))

print("SpatialPixelDataFrame runtime...")
system.time(m <- SpatialPixelsDataFrame(points = m[c("x", "y")], data=m))

print("SpatialGridDataFrame runtime...")
system.time(m <- as(m, "SpatialGridDataFrame"))

print("writeAsciiGrid runtime...")
system.time(writeAsciiGrid(m, paste(ifile, "example_grid.asc", sep=""), attr=3))