

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

parameters = matrix(data = c(c(0.2,1,0.2,0.1,0.1),
                             c(0.5,1,0.1,0.2,0.1),
                             c(0.5,1,0.2,0.1,0.7)),
                    ncol = 5, byrow = TRUE)

set.seed(666)
sim = lapply(1:3, FUN = function(x){
  out = list()
  out[[1]] = pbd_sim(pars = parameters[x, ], age = 10,
                     soc = 1, plotit = FALSE)
  out[[2]]  = pbd_sim(pars = parameters[x, ], age = 10,
                      soc = 2, plotit = FALSE)
  out = lapply(out, `[[`, 1)
  cat(paste0("DONE: parameter=", x, "\n"))
  return(out)
})


phy = c(lapply(sim, `[[`, 1), lapply(sim, `[[`, 2))

library(ape)
mltt.plot(phy = phy)


branches = lapply(phy, FUN = branching.times)



save(branches, phy, sim, file = "~/Dropbox/PhD/gitHub/protracted_sp/analysis/data/pbd_sim.RData")
