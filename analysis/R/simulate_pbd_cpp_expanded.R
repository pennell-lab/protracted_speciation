


library(ape)
library(PBD)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)


source("opt_pbd_sim_cpp.R")
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

theta0 = list("b1" = c(0.3, 0.5, 1, 3),
              "la1" = c(5, 2, 1, 0.5, 0.3),
              "b2" = c(0),
              "mu1" = c(0, 0.1, 0.2, 0.4, 1),
              "mu2" = c(0, 0.1, 0.7, 1, 2))

parameters = expand.grid(theta0)


tries = 5
ntaxa = 100
nphy = 10
ncore = 5
sim = lapply(1:nrow(parameters), FUN = function(s){
  cat(paste0(Sys.time(), "\t parameters: ", s, "\n"))
  cat(paste0(parameters[s,]))
  out = mclapply(X = 1:nphy, FUN = function(x){
    opt_pbd_sim_cpp(pars = as.matrix(parameters[s, ]), taxa = ntaxa, ntry = tries)
  }, mc.silent = TRUE, mc.cores = ncore)
  
  cat("\n\n")
  save(out, file = paste0("~/rcpp/pbd_sim_cpp", paste0(parameters[s,], collapse = "_"), ".RData"))
  return(out)
})











