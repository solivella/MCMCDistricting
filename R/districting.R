districting <- function(cutVector,
                        nDist,
                        basePoly,
                        basePoints = basePoints,
                        fullGraph,
                        fullEdgeMat,
                        objFunctions,
                        weights
                        ){
                        
  cutVector <- as.logical(cutVector)
  if(length(weights)>length(objFunctions)) stop("Number of weights cannot be greater than number of objective functions.\n",call.=FALSE)
    
  new.graph <- delete.edges(fullGraph,
                            E(fullGraph)[fullEdgeMat[cutVector,1]%--%
                                         fullEdgeMat[cutVector,2]])
  
  dist.info <- clusters(new.graph)

  if(dist.info$no!=nDist){
    return(-9999999999999)
  }

  district.index <- dist.info$membership
  
  district.plan <- unionSpatialPolygons(basePoly,district.index)
  
  points.data <- SpatialPointsDataFrame(basePoints,
  										as.data.frame(basePoly))
  
  new.data <- overlay(points.data,
  					  district.plan,
  					  fn=colSums)
  					  
  rownames(new.data) <- as.numeric(sub("[a-z]*",
  							"",
  							rownames(new.data),
  							ignore.case=TRUE))-1
  
  
  d.plan.data <- SpatialPolygonsDataFrame(district.plan,
  										  new.data)
  
  obj.eval <- sapply(objFunctions,
  					do.call,
  					args=list(plan=d.plan.data))
  
  result <- as.double(sum(obj.eval*weights))
  
  return(result)
}