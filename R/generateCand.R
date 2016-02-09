gen.cand <- function(curr.cuts){
  cand.cuts <- curr.cuts
  cand.cuts[sample(length(curr.cuts),1)] <- sample(0:1,1)
  return(as.double(cand.cuts))
}