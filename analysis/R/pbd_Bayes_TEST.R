




library(PBD)
library(TreeSim)
library(doParallel)

#load("analysis/data/pbd_sim.RData")
load("pbd_sim.RData")
sim.pars = apply(parameters, MARGIN = 2, mean)


# library(msm) for "dtnorm" - truncated normal
prior_b = function(b){
  dexp(b, rate = 1/sim.pars[1], log = TRUE)
}
prior_mu1 = function(mu1){
  dexp(mu1, rate = 1/sim.pars[4], log = TRUE)
}
prior_la1 = function(la1){
  dlnorm(la1, meanlog = log(sim.pars[2]), sdlog = log(15), log = TRUE)
}
prior_mu2 = function(mu2){
  dexp(mu2, rate = 1/sim.pars[5], log = TRUE)
}


rep = 1e+4
chains = 3
initial = sim.pars[c(1,4,2,5)]
files = paste0("Bayes_test_phy", 1:6)
times = vector(mode = "numeric")
#source("analysis/R/pbd_Bayes.R")
source("pbd_Bayes.R")
for(i in 1:6){
  for(k in 1:chains){
    foo = Sys.time()
    cat(paste0(foo, "\n"))
    pbd_Bayes(brts = branches[[i]],
              initparsopt = initial,
              prior_b = prior_b, prior_mu1 = prior_mu1, 
              prior_la1 = prior_la1, prior_mu2 = prior_mu2,
              step = 0.7, rep = rep, file = paste0(files[i], ".", k))
    times[i] = Sys.time() - foo
    cat(paste0(times[i], "\n"))
  }
}

# sum(partiu$accepted)/rep
# head(partiu)
# tail(partiu)

# library(reshape2)
# res = melt(partiu, id.vars = c("accepted", "proposed.par", "converged"))
# res$sim = rep(1:(rep+1), 7)
# head(res)
# library(ggplot2)
# ggplot(res, aes(sim, value)) + geom_line() +
#   facet_grid(variable ~ ., scales = "free_y")
# 
# 
# plot(partiu$posterior, type = "l")
# plot(partiu$logLik, type = "l")
# plot(partiu$prior, type = "l")
# plot(partiu$b, type = "l")
# plot(partiu$mu1, type = "l")
# plot(partiu$la1, type = "l")
# plot(partiu$mu2, type = "l")
# 
# hist(partiu$b, freq = FALSE, xlim = c(0, 10))
# curve(exp(prior_b(x)), col = "red", add = TRUE)






