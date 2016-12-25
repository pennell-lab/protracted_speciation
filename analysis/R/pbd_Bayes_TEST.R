




library(PBD)
library(TreeSim)

load("analysis/data/pbd_sim.RData")
sim.pars = apply(parameters, MARGIN = 2, mean)


library(msm)
prior_b = function(b){
  dtnorm(b, mean = sim.pars[1], sd = 15, log = TRUE, lower = 0)
}
curve(prior_b(x), -1, 100)
curve(exp(prior_b(x)), 0, 10)
prior_mu1 = function(mu1){
  dtnorm(mu1, mean = sim.pars[4], sd = 15, log = TRUE, lower = 0)
}
curve(prior_mu1(x), -1, 100)
prior_la1 = function(la1){
  dtnorm(la1, mean = sim.pars[2], sd = 15, log = TRUE, lower = 0)
}
curve(prior_la1(x), -1, 100)
prior_mu2 = function(mu2){
  dtnorm(mu2, mean = sim.pars[5], sd = 15, log = TRUE, lower = 0)
}
curve(prior_mu2(x), -1, 100)


rep = 5e+3
initial = sim.pars[c(1,4,2,5)]
source("analysis/R/pbd_Bayes.R")
partiu = pbd_Bayes(brts = branches[[1]],
                   initparsopt = initial,
                   prior_b = prior_b, prior_mu1 = prior_mu1, 
                   prior_la1 = prior_la1, prior_mu2 = prior_mu2,
                   step = 0.5, rep = rep)
sum(partiu$accepted)/rep
head(partiu)
tail(partiu)
library(reshape2)
res = melt(partiu, id.vars = c("accepted", "proposed.par", "converged"))
res$sim = rep(1:(rep+1), 7)
head(res)
library(ggplot2)
ggplot(res, aes(sim, value)) + geom_line() +
  facet_grid(variable ~ ., scales = "free_y")


plot(partiu$posterior, type = "l")
plot(partiu$logLik, type = "l")
plot(partiu$prior, type = "l")
plot(partiu$b, type = "l")
plot(partiu$mu1, type = "l")
plot(partiu$la1, type = "l")
plot(partiu$mu2, type = "l")

hist(partiu$b, freq = FALSE, xlim = c(0, 10))
curve(exp(prior_b(x)), col = "red", add = TRUE)






