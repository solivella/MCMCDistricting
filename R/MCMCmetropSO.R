##########################################################################
## The function takes samples of an objective distribution defined 
## over a space of feasible districting plans using the metropolis
## algorithm. 
##
## It is based on the function MCMCMetrop1R in package 'MCMCPack'.
## 
## Santiago Olivella. Washington University in St. Louis, 2011.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMCmetropSO" <- function(fun, theta.init,
                           cand.fun,
                           burnin, mcmc, thin,
                           verbose, seed=NA, ...){

  
  ## error checking here
  MCMCpack:::check.mcmc.parameters(burnin, mcmc, thin)

 
  
  ## form seed
  seeds <- MCMCpack:::form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  
  ## setup the environment so that fun can see the things passed as ...
  userfun <- function(ttt) fun(ttt, ...)
  my.env <- environment(fun = userfun)

  #candfun <- function(fff) cand.fun(fff, ...)
  my.env2 <- environment(fun = cand.fun)

  
  ## Call the C++ function to do the MCMC sampling 
  mcmc.sample <- .Call("MCMCmetropSO_cc",
                       userfun,
                       as.double(theta.init),
                       my.env,
                       cand.fun,
                       my.env2,
                       as.integer(burnin),
                       as.integer(mcmc),
                       as.integer(thin),
                       as.integer(verbose),
                       lecuyer=as.integer(lecuyer), 
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       PACKAGE="EITMdistricting")
  
  return(mcmc.sample)
}
 
