###################################################
## 
## This is the interface to the mcmc district plan
## sampling routine. Users should only interact
## with this function. 
##
## Santiago Olivella, 
## Washington University in St. Louis, 2011
##
###################################################


sample.plans <- function(base.polygons,
						base.data,
						n.districts,
						user.functions,
						fun.weights,
						burn.in=1000,
						thin.interval=10,
						mcmc.samples=50000){

	#library(sp)
	#library(spdep)
	#library(igraph)
	#library(MCMCpack)
	#library(gpclib)
	#library(maptools)

	gpclibPermit()

        cat("Preparing Spatial Data...\n")
	base.DataFrame <- SpatialPolygonsDataFrame(base.polygons,base.data)


	base.nb <- poly2nb(base.DataFrame)
	base.nb.dist <- nbdists(base.nb,coordinates(	base.DataFrame),longlat=TRUE)
	base.idw <-  lapply(base.nb.dist,function(x) x)
	base.nb.w <- nb2listw(base.nb,glist=base.idw,style="B")
	base.mat <- listw2mat(base.nb.w)

	base.graph <- graph.adjacency(base.mat,weighted=TRUE,mode="undirected")

	edge.n <- ecount(base.graph)
	vertex.n <- vcount(base.graph)

	edge.matrix <- cbind(base.graph[[3]],base.graph[[4]])

        cat("Sampling characteristic points within basic spatial units...\n")
	eigenPoints <- matrix(as.numeric(data.frame(lapply(base.DataFrame@polygons, 
									spsample,
									n=1,
									type="random"))),
						ncol=2,
						byrow=TRUE)

        cat("Beginning MCMC sampling routine.\n")
	post.samples <- MCMCmetropSO(districting,
								theta.init=sample(0:1,edge.n,replace=TRUE),
								cand.fun=gen.cand,
								fullGraph=base.graph,
								basePoly=base.DataFrame,
								nDist=n.districts,
								basePoints=eigenPoints,
								fullEdgeMat=edge.matrix,
								objFunctions=user.functions,
								weights=fun.weights,
								burnin=burn.in,
								thin=thin.interval,
								mcmc=mcmc.samples,
								verbose=1)

	mcmc.samples <- mcmc(post.samples[,1:edge.n],
						start=burn.in+1,
						end=mcmc.samples+burn.in,
						thin=thin.interval)

	sample.values <- post.samples[,edge.n+1]

	return(list(mcmc.samples=mcmc.samples,
				sample.values=sample.vales))
}





