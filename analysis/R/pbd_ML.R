###########################################
#   Function to fit the protracted sp model
# using Maximum Likelihood
##########################################


pbd_ML
#   "idparsopt" The ids of the parameters that must be optimized, e.g. 1:4 for all parameters. The ids are defined as follows:
#     - id == 1 corresponds to b (speciation-initiation rate)
#     - id == 2 corresponds to mu_1 (extinction rate of good species)
#     - id == 3 corresponds to la_1 (speciation-completion rate)
#     - id == 4 corresponds to mu_2 (extinction rate of incipient species)
#   "idparsfix" The ids of the parameters that should not be optimized, e.g. c(2,4) if mu_1 and mu_2 should not be optimized, but only b and la_1. In that case idparsopt must be c(1,3).
#   "parsfix" The values of the parameters that should not be optimized
#   "exteq" Sets whether incipient species have the same (1) or different (0) extinction rate as good species. If exteq = 0, then idparsfix and idparsopt should together have all parameters 1:4
#   "parsfunc" Specifies functions how the rates depend on time, default functions are constant functions missnumspec The number of species that are in the clade but missing in the phylogeny
#   "cond" Conditioning:
#     - cond == 0 : conditioning on stem or crown age
#     - cond == 1 : conditioning on stem or crown age and non-extinction of the phylogeny
#     - cond == 2 : conditioning on stem or crown age and number of extant taxa
#   "btorph" Sets whether the likelihood is for the branching times (0) or the phylogeny (1)
#   "soc" Sets whether the first element of the branching times is the stem (1) or the crown age methode Sets which method should be used in the ode-solver. Default is ’lsoda’.
#   "n_low" Sets the lower bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)
#   "n_up" Sets the upper bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)
#   "tol" Sets the tolerances in the optimization. Consists of:
#     - reltolx = relative tolerance of parameter values in optimization
#     - reltolf = relative tolerance of function value in optimization
#     - abstolx = absolute tolerance of parameter values in optimization
#   "maxiter" Sets the maximum number of iterations in the optimization
#   "optimmethod" Method used in optimization of the likelihood. Current default is ’subplex’. Alternative is ’simplex’ (default of previous versions)

#### Value
#   "b" gives the maximum likelihood estimate of b
#   "mu_1" gives the maximum likelihood estimate of mu_1
#   "la_1" gives the maximum likelihood estimate of la_1
#   "mu_2" gives the maximum likelihood estimate of mu_2
#   "loglik" gives the maximum loglikelihood
#   "df" gives the number of estimated parameters, i.e. degrees of feedom
#   "conv" gives a message on convergence of optimization; conv = 0 means convergence





library(PBD)
example("pbd_ML")
# pbd_ML(1:10,initparsopt = c(0.2,0.01,0.3), exteq = 1)

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
library(geiger)
const.df$ntip = sapply(constant, Ntip)
library(reshape2)
const2plot = melt(const.df, id.vars = c("loglik", "df", "conv", "ntip"))
library(ggplot2)
ggplot(const2plot, aes(ntip, value)) + geom_point() +
  facet_grid(variable ~ ., scales = "free_y")




set.seed(666)
slowdown = sim.bd.age(age = 10, numbsim = 100, 
                      lambda = 1.5, mu = 1,
                      K = 300, # sp limit
                      complete = FALSE)
slowdown = slowdown[sapply(slowdown, class) == "phylo"]
slow.branch = lapply(slowdown, FUN = function(x) x$edge.length)

init = c(2, 1, 1)
slow.estimates = lapply(slow.branch, FUN = pbd_ML, initparsopt = init)
slow.df = do.call(rbind, slow.estimates)
slow.df$ntip = sapply(slowdown, Ntip)
slow2plot = melt(slow.df, id.vars = c("loglik", "df", "conv", "ntip"))
ggplot(slow2plot, aes(ntip, value)) + geom_point() +
  facet_grid(variable ~ ., scales = "free_y")





both = rbind(const2plot, slow2plot)
both$model = c(rep("constant", nrow(const2plot)), 
               rep("slowdown", nrow(slow2plot)))
ggplot(both, aes(ntip, value, colour = model)) + geom_point() +
  facet_grid(variable ~ ., scales = "free_y")




