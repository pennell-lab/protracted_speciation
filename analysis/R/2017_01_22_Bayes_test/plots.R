


library(coda)

load("analysis/data/pbd_sim.RData")

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

plot.prior.post = function(coda, priors, par){
  layout(matrix(1:4, nrow = 4))
  for(i in 1:4){
    plot(density(coda[, i]), main = colnames(coda)[i])
    curve(exp(priors[[i]](x)), add = TRUE, col = "red")
    abline(v = parameters[i], col = "green")
  }
}



for(i in 1:6){
  fuck1 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".1_PBD_Bayes.txt"), sep = "\t")
  sum(fuck1$accepted, na.rm = 1)/nrow(fuck1)
  fuck2 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".2_PBD_Bayes.txt"), sep = "\t")
  sum(fuck2$accepted, na.rm = 1)/nrow(fuck2)
  fuck3 = read.csv(paste0("analysis/R/2017_01_22_Bayes_test/Bayes_test/Bayes_test_phy", i, ".3_PBD_Bayes.txt"), sep = "\t")
  sum(fuck3$accepted, na.rm = 1)/nrow(fuck3)
  
  fucks = list(fuck1[, 4:7], fuck2[, 4:7], fuck3[, 4:7])
  
  
  coda.fuck0 = lapply(X = fucks, FUN = mcmc, start = 3000, thin = 100)
  coda.fuck = mcmc.list(coda.fuck0)
  
  # plot(coda.fuck0[[1]])
  # parameters[1, ]
  # 
  # pairs(fuck1[, 4:7])
  # betterPairs(fuck[, 4:7])
  
  par.ind = ceiling(i/2)
  pdf(file = paste0("analysis/output/Bayes_test_phy_", i, ".pdf"))
  plot(coda.fuck)
  plot.prior.post(coda.fuck0[[1]], priors, 
                  par = parameters[par.ind, c(1,4,2,5)])
  plot.prior.post(coda.fuck0[[2]], priors,
                  par = parameters[par.ind, c(1,4,2,5)])
  plot.prior.post(coda.fuck0[[3]], priors,
                  par = parameters[par.ind, c(1,4,2,5)])
  par(mfrow=c(1,1))
  plot(phy[[i]])
  axisPhylo()
  dev.off()
}


