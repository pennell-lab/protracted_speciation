

library(ape)
library(coda)
library(reshape2)
library(ggplot2)

load("analysis/data/pbd_sim.RData")
sim.pars = apply(parameters, MARGIN = 2, mean)
source("analysis/R/betterPairs.R")

##### PRIORS
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
priors = list(prior_b, prior_mu1, prior_la1, prior_mu2)

plot.prior.post = function(coda, priors, par, ...){
  lcoda = length(coda)
  for(i in 1:4){
    plot(density(coda[[1]][, i]), main = colnames(coda[[1]])[i], ...)
    for(k in 2:lcoda){
      par(new = TRUE)
      plot(density(coda[[k]][, i]), lty = k,
           axes = FALSE, xlab = "", ylab = "",
           main = "", ...)
    }
    abline(v = parameters[i], col = "green")
    curve(exp(priors[[i]](x)), add = TRUE, col = "red")
    legend("topright", legend = c("prior", paste0("chain",1:lcoda)),
           lty = c(1, 1:lcoda), col = c("red", rep("black", lcoda)),
           cex = 0.7)
  }
}





for(i in 1:6){
  fuck1 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".1_PBD_Bayes.txt"), sep = "\t")
  acc1 = round(sum(fuck1$accepted, na.rm = 1)/nrow(fuck1), 3)
  fuck2 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".2_PBD_Bayes.txt"), sep = "\t")
  acc2 = round(sum(fuck2$accepted, na.rm = 1)/nrow(fuck2), 3)
  fuck3 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".3_PBD_Bayes.txt"), sep = "\t")
  acc3 = round(sum(fuck3$accepted, na.rm = 1)/nrow(fuck3), 3)
  
  
  
  fucks = list(fuck1[, 4:7], fuck2[, 4:7], fuck3[, 4:7])
  coda.fuck0 = lapply(X = fucks, FUN = mcmc, start = 3000, thin = 100)
  coda.fuck = mcmc.list(coda.fuck0)
  
  
  
  fuck.all = rbind(fuck1[, 1:7], fuck2[, 1:7], fuck3[, 1:7])
  #fuck.all = mcmc(fuck1[, 1:7], start = 3000, thin = 100)
  r = data.frame(fuck.all, "sim" = rep(1:length(fucks), each = nrow(fuck1)))
  r = melt(r, id = "sim")
  r$sim = as.factor(r$sim)
  r$steps = rep(1:nrow(fuck1), length(fucks))
  head(r)
  
  
  
  par.ind = ceiling(i/2)
  pdf(file = paste0("analysis/output/Bayes_test_phy_", i, ".pdf"))
  # the different chains combined
  plot(coda.fuck)
  
  # convergence/mix per variable - PLOT
  gelman.plot(coda.fuck0)
  # convergence/mix per variable - VALUES
  par(mfrow=c(1,1))
  plot(0, 0, xlim = c(0, 100), ylim = c(0, 100),
       type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "Convergence diagnostics - coda::gelman.diag")
  text(40+c(-20, 0, 20), 95, labels = c("","Point est.", "Upper C.I."),
       pos = 4)
  text(rep(40+c(0, 20), each = 4), 85+rep((-10)*0:3, 2),
       labels = as.character(round(gelman.diag(coda.fuck0)$psrf, 3)),
       pos = 4)
  text(20, 85+rep((-10)*0:3, 2), label = c("b", "mu1", "la1", "mu2")
       , pos = 4)
  text(20, 30, labels = paste0("Multivariate psrf = ",
                               round(gelman.diag(coda.fuck0)$mpsrf, 3)),
       pos = 4)
  text(20, 5, labels = paste("Acceptance ratio =", 
                             acc1, acc2, acc3, sep = "   "),
       pos = 4)
  
  # correlation between variables
  par(mfrow=c(1,1))
  betterPairs(fucks[[1]])
  betterPairs(fucks[[1]])
  betterPairs(fucks[[1]])
  
  # phylogeny
  plot(phy[[i]])
  axisPhylo()
  
  # comparison between priors and posteriors
  layout(matrix(1:4, nrow = 4))
  plot.prior.post(coda.fuck0, priors,
                  par = parameters[par.ind, c(1,4,2,5)])
  
  
  dev.off()
  
  
  
  # trace plots for each chain separately
  par(mfrow=c(1,1))
  ggplot(r, aes(steps, value, colour = sim)) + geom_line(linetype = 2) +
    facet_grid(variable ~ ., scales = "free_y")
  ggsave(file = paste0("analysis/output/Bayes_test_trace_plots_phy_", i, ".pdf"))
}




