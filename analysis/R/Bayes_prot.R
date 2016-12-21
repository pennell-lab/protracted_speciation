




library(PBD)
library(TreeSim)
set.seed(666)
constant = sim.bd.age(age = 10, numbsim = 100, 
                      lambda = 1.5, mu = 1, 
                      complete = FALSE)
constant = constant[sapply(constant, class) == "phylo"]
const.branch = lapply(constant, FUN = function(x) x$edge.length)
init = c(2, 1, 1)
const.estimates = lapply(const.branch, FUN = pbd_ML, initparsopt = init)
const.df = do.call(rbind, const.estimates)



library(msm)
prior_b = function(b){
  dtnorm(b, mean = 2, sd = 15, log = TRUE, lower = 0)
}
curve(prior_b(x), -1, 100)
prior_mu1 = function(mu1){
  dtnorm(mu1, mean = 1, sd = 15, log = TRUE, lower = 0)
}
curve(prior_mu1(x), -1, 100)
prior_la1 = function(la1){
  dtnorm(la1, mean = 1, sd = 15, log = TRUE, lower = 0)
}
curve(prior_la1(x), -1, 100)


source("pbd_Bayes.R")
partiu = pbd_Bayes(brts = const.branch[[1]],
                   initparsopt = c(2, 1, 1),
                   prior_b, prior_mu1, prior_la1,
                   step = 0.5, rep = 1e+4)
res = partiu[partiu$converged, ]
head(res)
tail(res)
plot(partiu$posterior, type = "l")
plot(partiu$logLik, type = "l")
plot(partiu$prior, type = "l")
plot(partiu$b, type = "l")
plot(partiu$mu1, type = "l")
plot(partiu$la1, type = "l")
plot(partiu$mu2, type = "l")
plot(cumsum(partiu$accepted), type = "l")







