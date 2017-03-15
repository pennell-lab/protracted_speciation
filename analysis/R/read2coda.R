


library(coda)


files = list.files(pattern = "*PBD_Bayes.txt")


chains = 10
my.mcmc = list()
nphy = length(files)/chains
pp = matrix(ncol = 5, nrow = nphy, dimnames = list(1:nphy, c("b1", "la1", "b2", "mu1", "mu2")))
for(k in 1:nphy){
  foo = gsub(pattern = "Bayes_test_param_", replacement = "", files[((k-1)*chains)+1])
  foo = gsub(pattern = "_PBD_Bayes.txt", replacement = "", foo)
  pp[k, ] = as.numeric(strsplit(foo, "_")[[1]][1:5])
  aux = list()
  for(i in 1:chains){
    aux[[i]] = read.csv(files[((k-1)*chains)+i], sep = "\t", row.names = NULL)
    aux[[i]] = as.mcmc(aux[[i]][ , 1:7])
  }
  my.mcmc[[k]] = as.mcmc.list(aux)
}
names(my.mcmc) = apply(pp, MARGIN = 1, paste, collapse = "_")

source("analysis/R/plot_functions.R")
plot.prior.post(aux[[1]])
