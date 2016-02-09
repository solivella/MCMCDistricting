.onAttach <- function(...) {
 
   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   cat("##\n## Empirical Implications of Theoretical Models 
   \n## Districting Package (EITMdistricting).
   \n##
   \n## THIS SOFTWARE USES THE gpclib PACKAGE, WHICH 
   \n## MAY ONLY BE USED FREELY FOR NON-COMMERICAL PURPOSES. 
   \n## ALSO, AND IN ACCORDANCE WITH THE GPL-2, NOTE THAT
   \n## THIS PACKAGE CONTAINS A FUNCTION CALLED `MCMCmetropSO',
   \n## WHICH BUILDS OFF OF `MCMCmetrop1R' IN PACKAGE MCMCpack.
   \n##
   \n## (Santiago Olivella, WashU, 2011)
   \n##\n")
   
   require(coda, quietly=TRUE)
   require(sp, quietly=TRUE)
   require(spdep, quietly=TRUE)
   require(gpclib, quietly=TRUE)
   require(MCMCpack, quietly=TRUE)
   require(maptools, quietly=TRUE)
   require(igraph, quietly=TRUE)      
}

.onUnload <- function(libpath) {
    library.dynam.unload("EITMdistricting", libpath)
}

