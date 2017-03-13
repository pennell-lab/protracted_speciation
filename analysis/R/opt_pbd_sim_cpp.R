




opt_pbd_sim_cpp = function (pars, age = NULL, taxa = NULL, ntry = 1, soc = 2) 
  # ntry - controls the number of times the simulation will try to generate a valid phylogeny
  # soc	- sets whether the 'age' is the stem (1) or crown (2) age
{
  require(ape)
  require(PBD)
  
  require(Rcpp)
  require(RcppArmadillo)
  sourceCpp("opt_pbd_sim.cpp")
  sourceCpp("opt_pbd_sim_taxa.cpp")
  
  
  # la1 = pars[1]
  # la2 = pars[2]
  # la3 = pars[3]
  # mu1 = pars[4]
  # mu2 = pars[5]
  # need to be ordered accordingly 
  pars = pars[c(1, 4, 2, 3, 5)]
  
  if(!is.null(age)){
    L = pbdLoop(pars, age, ntry, soc)
    L[1, 4] = 0
  } else if(!is.null(taxa)){
    L = pbdLoop_taxa(pars, taxa, ntry)
    age = L[1, 4]
    L[1, 4] = 0
  } else{
    stop("Either 'age' or 'taxa' must be given!")
  }
  
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









