

library(doParallel)
library(PBD)


#### pbd_sim()
#   -pars --> Vector of parameters: 
# pars[1] corresponds to b_1, the speciation-initiation rate of good species 
# pars[2] corresponds to la_1, the speciation-completion rate 
# pars[3] corresponds to b_2, the speciation-initiation rate of incipient species 
# pars[4] corresponds to mu_1, the extinction rate of good species 
# pars[5] corresponds to mu_2, the extinction rate of incipient species 
#   -age --> Sets the age for the simulation
#   -soc --> Sets whether this age is the stem (1) or crown (2) age
#   -plotit --> Sets whether the various trees produced by the function should be plotted or not

parameters = matrix(data = c(c(0.3,1,0,0.1,0.1),
                             c(0.3,1,0,0.2,0.1),
                             c(0.5,1,0,0.1,0.7),
                             c(0.3,0.5,0,0.1,0.1)),
                    ncol = 5, byrow = TRUE)


nphy = 10
sim = lapply(1:nrow(parameters), FUN = function(s){
  cat(paste0("parameter", s, "\t", Sys.Date(), "\n"))
  #out = mclapply(X = 2:nphy, FUN = function(x){
  out = lapply(X = 1:nphy, FUN = function(x){
    cat(paste0(x, "\t", Sys.Date()))
    pbd_sim(pars = parameters[s, ], age = 50,
            soc = 1, plotit = FALSE)
    cat("\n")
  })
  save(out, file = paste0("pbd_sim", s, ".RData"))
  cat(paste0("\n", "DONE", s, "\n\n\n"))
  #mc.silent = TRUE, mc.set.seed = 666, mc.cores = nphy/2)
  return(out)
})


phy = lapply(sim, `[[`, 1)


branches = lapply(phy, FUN = branching.times)



#save(branches, parameters, phy, sim, file = "~/Dropbox/PhD/gitHub/protracted_sp/analysis/data/pbd_sim_best_scenario.RData")


save(branches, parameters, phy, sim, file = "pbd_sim_best_scenario.RData")






