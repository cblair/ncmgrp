require(survival)
require(maptools)
require(sp)
require(snow)

#create cluster for processing
cluster.nodes <- c("localhost", "adam@mot", "adam@dataserver")
#cluster.nodes <- c("localhost", "adam@dataserver")
c1 <- makeCluster(cluster.nodes, type="MPI")



stopCluster(c1)