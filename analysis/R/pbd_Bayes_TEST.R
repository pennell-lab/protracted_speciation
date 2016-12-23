




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
head(const.df)




library(msm)
prior_b = function(b){
  dtnorm(b, mean = const.df[1, 1], sd = 15, log = TRUE, lower = 0)
}
curve(prior_b(x), -1, 100)
curve(exp(prior_b(x)), 0, 10)
prior_mu1 = function(mu1){
  dtnorm(mu1, mean = const.df[1, 2], sd = 15, log = TRUE, lower = 0)
}
curve(prior_mu1(x), -1, 100)
prior_la1 = function(la1){
  dtnorm(la1, mean = const.df[1, 3], sd = 15, log = TRUE, lower = 0)
}
curve(prior_la1(x), -1, 100)
prior_mu2 = function(mu1){
  dtnorm(mu1, mean = const.df[1, 4], sd = 15, log = TRUE, lower = 0)
}
curve(prior_mu2(x), -1, 100)


rep = 5e+3
initial = const.df[1, 1:4] + c(-0.7, +0.13, -90, +1.02)
source("analysis/R/pbd_Bayes.R")
partiu = pbd_Bayes(brts = const.branch[[1]],
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






