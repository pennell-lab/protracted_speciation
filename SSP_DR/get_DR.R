

# auxiliary functions to estimate DR
# from the SM of Jetz etal. 2012:
#       The ES measure for a focal tip on a rooted bifurcating tree is the sum of the edge lengths from the species i
# to the root, with each consecutive edge discounted by a factor of 1/2:
#   ESi = sum (from j = 1 to Ni) [ lj * ( 1/2^(j-1) ) ] (Supplementary equation 1)
# Ni is the number of edges on the path from species i and the root, and lj is the length of the edge j, with j=1
# being the pendant edge leading to the species and j=Ni being the edge nearest the root on that path. For a more
# general equation for non-bifurcating trees, see Redding et al. (2008).
#        The inverse of this measure can be seen as a measure of the splitting rate of the path to a tip: species in 
# rapidly-diversifying clades will have short edge lengths shared among many species and low ES values, while isolated 
# species on a tree have no evidence of recent diversification and large ES values. We term the 1/ES metric for a
# species its species-level lineage diversification rate, or DR. 
require(ape)
get_DR = function(phy){
  nt = Ntip(phy)
  path = get_path(1:nt, phy$edge, nt+1)
  ES = estimateES(path, edge.length = phy$edge.length)
  DR = 1/ES
  names(DR) = phy$tip.label
  
  return(DR)
}
Get_Path = function(tip, edge, root){
  out = c()
  while(tip != root){
    out = c(which(edge[ , 2] == tip), out)
    tip = edge[out[1], 1]
  }
  
  return(rev(out))
}
get_path = Vectorize(FUN = Get_Path, vectorize.args = "tip")
EstimateES = function(ind, edge.length){
  lj = edge.length[ind]
  weight = 2^( (1:length(ind)) - 1 )
  out = sum(lj * ( 1/weight ))
  
  return(out)
}
estimateES = Vectorize(FUN = EstimateES, vectorize.args = "ind")