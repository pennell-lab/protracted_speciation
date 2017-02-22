







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

theta0 = list("b1" = c(0.3, 0.5, 1),
              "la1" = c(5, 2, 1, 0.5),
              "b2" = c(0),
              "mu1" = c(0.1, 0.2, 0.4, 1),
              "mu2" = c(0.1, 0.7, 1))

parameters = expand.grid(theta0)



nphy = 10
sim = lapply(1:nrow(parameters), FUN = function(s){
  out = mclapply(X = 1:nphy, FUN = function(x){
    pbd_sim(pars = parameters[s, ], age = 50,
            soc = 1, plotit = FALSE)
  }, mc.silent = TRUE, mc.cores = nphy)

  save(out, file = paste0("pbd_sim", s, ".RData"))
  return(out)
})


phy = lapply(sim, `[[`, 1)


branches = lapply(phy, FUN = branching.times)



save(branches, parameters, phy, sim, file = "pbd_sim_expanded.RData")













