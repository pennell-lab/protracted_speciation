
library(doParallel)
library(PBD)

#To assess the robustness (particularly bias) of the ML approach in yielding the correct parameter values, we simulated 1000 phylo- genetic trees under the protracted speciation model for various sets of parameters (b=0.5, λ=0.1,0.3,1, μ1 =μ2 =μ= 0, 0.1, 0.2), and a fixed crown age of 5, 10, or 15 My, conditional on the realized tree retaining the initial root (i.e., survival of both original crown lineages). We used the function pbd_sim_cpp to perform these simulations. This function uses the coalescent process to efficiently simulate branching times (see Lambert et al. 2014). We then estimated the parameters using ML

cores = 12
Nphy = 1000
parameters = expand.grid(data.frame(b1 = 0.5,
                                    mu1 = c(0,0.1,0.2),
                                    lambda = c(0.1,0.3,1)))
parameters$mu2 = parameters$mu1
simulations = mclapply(1:nrow(parameters), FUN = function(k){
  sim = list()
  for(i in 1:Nphy){
    while(1){
      aux = try(pbd_sim_cpp(pars = parameters[k , ], soc = 2, age = 15, plotltt = 0))
      if(class(aux) == "numeric") if(max(aux) == 15){
        sim[[i]] = aux
        break()
      }
    }
  }
  return(sim)
}, mc.cores = cores)
save(simulations, file = "etienne_simulations3par.RData")

ml.estimates = mclapply(simulations, FUN = function(b){
  do.call(rbind, lapply(b, FUN = pbd_ML,
                        initparsopt = c(0.2,0.1,1),
                        exteq = 1,
                        btorph = 0,
                        soc = 2))
}, mc.cores = cores)
save(ml.estimates, file = "etienne_ML3par.RData")
rm(simulations)


