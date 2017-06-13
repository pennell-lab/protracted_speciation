
library(doParallel)
library(PBD)

cores = 16
parameters = expand.grid(data.frame(b = 0.5, mu1 = c(0,0.1,0.2), lambda = c(0.1,0.3,1), mu2 = c(0,0.1,0.2)))
simulations = mclapply(1:nrow(parameters), FUN = function(k){
  sim = list()
  for(i in 1:1000){
    while(1){
      aux = try(pbd_sim_cpp(pars = parameters[k , ], soc = 2, age = 15, plotltt = 0))
      if(class(aux) == "numeric") if(max(aux) == 15){
        sim[[i]] = aux
        break()
      }
    }
  }
}, mc.cores = cores)
save(simulations, file = "etienne_simulations.RData")

ml.estimates = mclapply(simulations, FUN = function(b){
  do.call(rbind, lapply(b, FUN = pbd_ML, initparsopt = c(0.2,0.1,1,0.1), exteq = 0, btorph = 0, soc = 2))
}, mc.cores = cores)
save(ml.estimates, file = "etienne_ML.RData")
