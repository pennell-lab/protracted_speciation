


# from Phillimore 2010 "SUBSPECIES ORIGINATION AND EXTINCTION IN BIRDS"
# Given a subspeciation rate (lambda) and subspecies extinction rate (mu), I estimated the likelihood (l_i = P(n|t)) of observing n subspecies in species i of age t following Equation 3 and Equation 4 (Bokma 2003; Ricklefs 2007, 2009), where the relative extinction rate epsilon = mu / lambda.
#   E(n) = e^(lambda * (1 - epsilon) * t) = epsilon^((lambda - mu) * t)                Eq.3
#   P(n|t) = (1 - epsilon) * ( [E(n) - 1]^(n-1) / [E(n) - epsilon]^n )                  Eq.4
phillimore_model = function(n, t, rate, epsilon){
  LL = function(rate, epsilon){
    En = exp(rate*t)
    Pnt = (1 - epsilon) * ((En - 1)^(n - 1)) / ((En - epsilon)^n)
    return(-sum(log(Pnt)))
  }
  
  return(LL)
}





# https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
# Posted on August 9, 2012 by susanejohnston
ggplotRegression <- function (fit) {
  require(ggplot2)
  if(length(class(fit)) == 1){ # it's a regular "lm"
    R2 = summary(fit)$adj.r.squared
  } else{
    R2 = cor(fit$y,predict(fit))^2
  }
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(R2, 5),
                       "Intercept =",signif(coef(fit)[1],5 ),
                       " Slope =",signif(coef(fit)[2], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) 
}





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

######### Gets the ssp given a full matrix "L"
get_ssp = function(mat, taxa = NULL, tree = NULL, patt = NULL){
  if(!is.null(taxa)){
    out = sapply(taxa, FUN = function(s) length(which(mat[ , 6] == s)) - 1)
    nn = taxa
  } else if(!is.null(tree)){
    if(is.null(patt)){
      patt = "sp."
    }
    taxa = gsub(pattern = patt, replacement = "", tree$tip.label)
    out = sapply(taxa, FUN = function(s) length(which(mat[ , 6] == s)) - 1)
    nn = paste0(patt, taxa)
  }
  if(ncol(mat) == 7){
    out2 = sapply(taxa,  FUN = function(s){
      aa = mat[which(mat[ , 6] == s), 7]
      aa[1]
    })
    out = data.frame(ssp = out, dynamics = out2, check.rows = TRUE)
    row.names(out) = nn
  } else{
    names(out) = nn
  }
  
  return(out)
}

######### Gets the SP age (i.e., most recent branch length)
get_sp_age = function(phy){
  n = Ntip(phy)
  id = which(phy$edge[ , 2] %in% 1:n)
  out = phy$edge.length[id]
  names(out) = phy$tip.label
  
  return(out)
}