
library(mailR)
send = function(x){
  send.mail(from = "brejafimano@gmail.com",
            to = "maurotcs@gmail.com",
            subject = "First_try",
            body = x,
            smtp = list(host.name = "smtp.gmail.com", port = 465, 
                        user.name = "brejafimano", passwd = "carnivora1", 
                        ssl = TRUE),
            authenticate = TRUE,
            send = TRUE)
}

start = Sys.time()

library(PBD)
library(TreeSim)
library(parallel)
set.seed(666)
constant = sim.bd.age(age = 10, numbsim = 50, 
             lambda = 1.5, mu = 1, 
             complete = FALSE)
constant = constant[sapply(constant, class) == "phylo"]
const.branch = lapply(constant, FUN = function(x) x$edge.length)
init = c(2, 1, 1)
const.estimates = lapply(const.branch, FUN = pbd_ML, initparsopt = init)
const.df = do.call(rbind, const.estimates)
head(const.df)
save(constant, const.df, file = "analysis/data/sim.bd.age_age10_lamb1.5_mu1.RData")


loop1 = Sys.time()
send(x = paste("Simulation of phylogenies finished at", loop1, "with", length(constant), "phylogenies."))







rep = 1e+5
initial = list()
set.seed(666)
for(i in 1:nrow(const.df)){
  random = c(runif(2), runif(1, min = 0, max = 15), runif(1))
  initial[[i]] = const.df[i, 1:4] + random
}
save(const.branch, initial, file = "analysis/data/parameters_first_try.RData")


library(msm)
first_try = list()
source("analysis/R/pbd_Bayes.R")
for(k in 1:length(const.branch)){
  prior_b = function(b){
    dtnorm(b, mean = as.numeric(initial[[k]][1]), 
           sd = 15, log = TRUE, lower = 0)
  }
  prior_mu1 = function(mu1){
    dtnorm(mu1, mean = as.numeric(initial[[k]][2]), 
           sd = 15, log = TRUE, lower = 0)
  }
  prior_la1 = function(la1){
    dtnorm(la1, mean = as.numeric(initial[[k]][3]),
           sd = 15, log = TRUE, lower = 0)
  }
  prior_mu2 = function(mu1){
    dtnorm(mu1, mean = as.numeric(initial[[k]][4]),
           sd = 15, log = TRUE, lower = 0)
  }
  
  first_try[[k]] = pbd_Bayes(brts = const.branch[[k]],
                             initparsopt = initial[[k]],
                             prior_b = prior_b, prior_mu1 = prior_mu1, 
                             prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                             step = 0.5, rep = rep)
  
  save(first_try, file = "analysis/data/results_first_try.RData")
}


loop2 = Sys.time()
send(x = paste("MCMC has finished at", loop2, ".", length(first_try), "phylogenies were analyzed."))




sapply(first_try, FUN = function(x) sum(x$accepted)/rep)




