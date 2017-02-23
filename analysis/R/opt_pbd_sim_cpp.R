library(ape)
library(PBD)
source("analysis/R/opt_pbd_sim.cpp")

library(Rcpp)
library('inline')


pars = c(0.3, 2, 0, 0.1, 0.1)
age = 10
soc = 2
set.seed(666)
(mm = opt_pbd_sim_cpp(pars = c(0.3, 2, 0, 0.1, 0.1), age = 5))
(tt = opt_pbd_sim_cpp(pars = c(0.3, 2, 0, 0.1, 0.1), taxa = 10))
plot(tt[[1]])
set.seed(660)
pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 3, soc = 1, plotit = FALSE)

library(microbenchmark)
compare <- microbenchmark(opt_pbd_sim_cpp(pars = c(0.3, 2, 0, 0.1, 0.1), age = 5),
                          pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 5, soc = 1, plotit = FALSE),
                          times = 100)
cpp = list()
pbd = list()
pp = c(0.3, 2, 0, 0.1, 0.1)
for(i in 1:100){
  cpp[[i]] = opt_pbd_sim_cpp(pars = pp, age = 5)
  pbd[[i]] = pbd_sim(pars = pp, age = 5, soc = 1, plotit = FALSE)
}
length(cpp)
length(pbd)
sum(sapply(lapply(cpp, '[[', 1), is.null))/100
sum(sapply(lapply(pbd, '[[', 1), is.null))/100


sourceCpp("analysis/R/opt_pbd_sim.cpp")
sourceCpp("analysis/R/opt_pbd_sim_taxa.cpp")
opt_pbd_sim_cpp = function (pars, age = NULL, taxa = NULL) # soc = 2, ntry = 1) 
  # ntry - controls the number of times the simulation will try to generate a valid phylogeny
{
  # la1 = pars[1]
  # la2 = pars[2]
  # la3 = pars[3]
  # mu1 = pars[4]
  # mu2 = pars[5]
  # need to be ordered accordingly 
  pars = pars[c(1, 4, 2, 3, 5)]
  
  if(!is.null(age)){
    L = pbdLoop(pars, age)
  } else if(!is.null(taxa)){
    L = L = pbdLoop_taxa(pars, taxa)
    age = L[1, 4]
    L[1, 4] = 0
  } else{
    stop("Either 'age' or 'taxa' must be given!")
  }
  
  L0 = L
  absL = L
  absL[, 2] = abs(L[, 2])
  ttree = try(PBD::detphy(absL, age), silent = TRUE)
  if(class(ttree) == "try-error"){
    return(NULL)
  }
  rtree = try(ape::read.tree(text = ttree))
  if(is.null(rtree) | class(rtree) == "try-error"){
    return(NULL)
  }
  tree = as.phylo(rtree)
  L[, 3:5][which(L[, 3:5] == -1)] = age + 1
  L[, 3:5] = age - L[, 3:5]
  L = L[order(L[, 1]), ]
  Ltreeslist = list(tree = tree, 
                    L = L,
                    L0 = L0)
  return(Ltreeslist)
}
