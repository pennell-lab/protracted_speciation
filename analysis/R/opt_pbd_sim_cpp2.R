
source("get_phylo.R")
try = lapply(1:10, function(i) opt_pbd_sim_cpp2(1, 0, 1, 10, 0, 0.5, taxa = 100))
Ls = lapply(try, "[[", 1)
step1 = lapply(Ls, twoL)
step2 = lapply(step1, function(x) lapply(x, get_phylo))
phy = mapply(combineL, phy1 = lapply(step2, "[[", 1), phy2 = lapply(step2, "[[", 2), SIMPLIFY = FALSE)
plot(phy[[1]]);axisPhylo()


#phy1=aa[[1]];phy2=aa[[2]]
combineL = function(phy1, phy2){
  Nnode = phy1$Nnode + phy2$Nnode + 1
  tips = sapply(list(phy1, phy2), Ntip)
  # use the 1st phylo as reference
  edge = phy1$edge
  
  # modify "edge" to include the nodes and tips from the other phylo
  edge[ , 1] = edge[ , 1] + tips[2] + 1
  edge[ , 2] = sapply(edge[ , 2], function(x) ifelse(x > tips[1], x + tips[2] + 1, x))
  edge = rbind(c(sum(tips) + 1, sum(tips) + 2), edge) # add node that conects the 2 phylos
  # modify the edge from phylo2 to include the nodes and tips from the other phylo
  ee = phy2$edge
  ee[ , 1] = ee[ , 1] + max(edge) - tips[2]
  ee[ , 2] = sapply(ee[ , 2], function(x) ifelse(x > tips[2], x + max(edge) - tips[2], x + tips[1]))
  edge = rbind(edge, c(sum(tips) + 1, max(edge) + 1), ee) # paste the 2 edges, adding the connection to the second phylo
  
  # creates the edge.length
  edge.length = c(phy1$root.edge, phy1$edge.length, phy2$root.edge, phy2$edge.length)
  
  out = list("edge" = edge, "edge.length" = edge.length, "Nnode" = Nnode, "tip.label" = paste0("sp.", 1:sum(tips)))
  class(out) = "phylo"
  
  return(out)
}
#aa=twoL(try[[1]][[1]])
twoL = function(L){
  out = list()
  for(i in 1:2){
    ind = sort(open_tree(L, i))
    out[[i]] = L[c(i, ind), ]
    sp = length(unique(out[[i]][ , 6]))
    or = order(out[[i]][ , 6])
    mul = table(out[[i]][ , 6])
    out[[i]][or, 6] = rep(1:sp, mul)
  }
  out[[2]][1, 1] = 1
  out[[2]][ , 2] = sapply(out[[2]][ , 2], function(x) ifelse(x == 2, 1, x))
  
  names(out) = c("L1", "L2")
  
  return(out)
}
#mat=L;ind=1
open_tree = function(mat, ind){
  out = which(mat[ , 2] == mat[ind, 1])
  if(length(out) > 0){
    for(i in 1:length(out)){
      aux2 = open_tree(out[i], mat = mat)
      out = c(out, aux2)
    }
  }
  return(out)
}


opt_pbd_sim_cpp2 = function (b1, mu1, la1, la2, b2, mu2,
                             age = NULL, taxa = NULL, ntry = 1, soc = 2) 
  # ntry - controls the number of times the simulation will try to generate a valid phylogeny
  # soc	- sets whether the 'age' is the stem (1) or crown (2) age
{
  require(ape)
  require(PBD)
  
  require(Rcpp)
  require(RcppArmadillo)
  
  sourceCpp("opt_pbd_sim_taxa2tts.cpp")
  
  
  # la1 = pars[1]
  # la2 = pars[2]
  # la3 = pars[3]
  # mu1 = pars[4]
  # mu2 = pars[5]
  # need to be ordered accordingly 
  pars = c(b1, mu1, la1, la2, b2, mu2)
  
  L = pbdLoop_taxa2tts(pars, taxa, ntry)
  age = L[1, 4]
  L[1, 4] = 0
  
  
  if(sum(L[1, ]) == 0){
    return(NULL)
  }
  
  L0 = L
  # absL = L
  # absL[, 2] = abs(L[, 2])
  # absL[1, 2] = 0
  # ttree = try(PBD::detphy(absL, age), silent = TRUE)
  # if(class(ttree) == "try-error"){
  #   return(NULL)
  # }
  # rtree = try(ape::read.tree(text = ttree))
  # if(is.null(rtree) | class(rtree) == "try-error"){
  #   return(NULL)
  # }
  # tree = ape::as.phylo(rtree)
  
  # sL_random = PBD::sampletree(absL, age, samplemethod = "random")
  # stree_random = ape::as.phylo(ape::read.tree(text = PBD::detphy(sL_random, age)))
  # sL_oldest = PBD::sampletree(absL, age, samplemethod = "oldest")
  # stree_oldest = ape::as.phylo(ape::read.tree(text = PBD::detphy(sL_oldest, age)))
  # sL_youngest = PBD::sampletree(absL, age, samplemethod = "youngest")
  # stree_youngest = ape::as.phylo(ape::read.tree(text = PBD::detphy(sL_youngest, age)))
  # sL_random[, 3:5][which(sL_random[, 3:5] == -1)] = age + 1
  # sL_random[, 3:5] = age - sL_random[, 3:5]
  # sL_random = sL_random[order(sL_random[, 1]), ]
  # sL_oldest[, 3:5][which(sL_oldest[, 3:5] == -1)] = age + 1
  # sL_oldest[, 3:5] = age - sL_oldest[, 3:5]
  # sL_oldest = sL_oldest[order(sL_oldest[, 1]), ]
  # sL_youngest[, 3:5][which(sL_youngest[, 3:5] == -1)] = age + 1
  # sL_youngest[, 3:5] = age - sL_youngest[, 3:5]
  # sL_youngest = sL_youngest[order(sL_youngest[, 1]), ]
  # reconL = PBD::pbd_reconstruct(L0)
  # recontree = ape::as.phylo(ape::read.tree(text = PBD::detphy(reconL, age)))
  L[, 3:5][which(L[, 3:5] == -1)] = age + 1
  L[, 3:5] = age - L[, 3:5]
  L = L[order(L[, 1]), ]
  
  # reconL[, 3:5][which(reconL[, 3:5] == -1)] = age + 1
  # reconL[, 3:5] = age - reconL[, 3:5]
  # reconL = reconL[order(reconL[, 1]), ]
  
  Ltreeslist = list(# tree = tree,
    # stree_random = stree_random, 
    # stree_oldest = stree_oldest,
    # stree_youngest = stree_youngest, 
    # recontree = recontree, 
    # reconL = reconL,
    # sL_random = sL_random,
    # sL_oldest = sL_oldest, 
    # sL_youngest = sL_youngest,
    L = L,
    L0 = L0)
  return(Ltreeslist)
}









