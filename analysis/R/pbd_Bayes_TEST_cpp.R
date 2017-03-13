




library(PBD)
library(ape)
library(doParallel)

path = "data/"
#sim.files = list.files(path = path, pattern = ".RData")
load(paste0(path, "pbd_sim_cpp_full.RData"))
parameters = names(sim)

#source("analysis/R/get_phylo.R")
source("get_phylo.R")
bool = lapply(sim, sapply, function(x) sum(x[[1]])==0)
clean = mapply(function(x, b){
  if(sum(b) == 100){
    return(NULL)
  } else{
    return(x[!b])
  }
}, x = sim, b = bool, SIMPLIFY = FALSE)
as.numeric(sapply(clean, length))
clean = clean[as.numeric(sapply(clean, length)) > 0]
as.numeric(sapply(clean, length))
phylos = lapply(clean, lapply, function(x){
  age = x$L[1, 3]
  phy = get_phylo(x$L, age)
  return(phy)
})
branches = lapply(phylos, function(x) sapply(x, branching.times))

save(clean, phylos, file = "pbd_sim_cpp_phy.RData")



rep = 5e+4
chains = 10
files = sapply(1:length(phylos), function(x) paste0("Bayes_test_param_", parameters[x], "_phy_", 1:length(phylos[[x]])))
#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")
for(i in 1:length(branches)){
  sim.pars = as.numeric(strsplit(parameters[i], "_")[[1]])
  
  prior_b = function(b){
    dexp(b, rate = 1/sim.pars[1], log = TRUE)
  }
  prior_mu1 = function(mu1){
    dexp(mu1, rate = 1/max(0.01, sim.pars[4]), log = TRUE)
  }
  prior_la1 = function(la1){
    dlnorm(la1, meanlog = log(sim.pars[2]), sdlog = log(15), log = TRUE)
  }
  prior_mu2 = function(mu2){
    dexp(mu2, rate = 1/max(0.01, sim.pars[5]), log = TRUE)
  }
  
  for(p in 1:length(phylos[[i]])){
    mclapply(1:chains, FUN = function(xxx){
      pbd_Bayes(brts = branches[[i]],
                initparsopt = rep(k, 4),
                prior_b = prior_b, prior_mu1 = prior_mu1, 
                prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                step = 0.7, rep = rep, file = paste0(files[[i]][p], ".", xxx))
    }, mc.cores = chains)
  }
}
