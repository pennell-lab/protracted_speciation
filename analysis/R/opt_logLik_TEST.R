
# eg
# pbd_loglik(pars1 = c(0.2,0.1,1,0.1), pars2 = c(1,1,2,0,"lsoda"),brts = 1:10)

## DATA
input = matrix(data = c(0.2,0.1,1,0.1,
                        0.18,0.1,1,0.15,
                        0.18,0.07,1.2,0.1), 
               ncol = 4, byrow = TRUE)
input
mat14 = matrix(c(0, 1, 2, 0,
                 1, 0, 2, 0,
                 1, 1, 1, 0,
                 2, 1, 1, 0,
                 1, 1, 2, 0,
                 1, 1, 2, 0),
               ncol = 4, byrow = TRUE)
mat67 = matrix(0, nrow = nrow(mat14),
               ncol = 2, byrow = TRUE)
pars2 = data.frame(mat14,
                   "method" = c(rep("lsoda", 4), "lsodes", "vode"),
                   mat67,
                   stringsAsFactors = FALSE)
pars2 # the DEFAULT is c(1, 1, 2, 1, "lsoda", 0, 0)


load("analysis/data/sim.bd.age_age10_lamb1.5_mu1.RData")
library(geiger)
branch = lapply(constant, FUN = branching.times)
branch[[length(branch)+1]] = 1:10



# MY function
source("analysis/R/opt_logLik.R")
# PBD's function
library(PBD)
my.results = numeric()
pbd.results = numeric()
check = data.frame("pars2" = numeric(nrow(pars2)*length(branch)),
                   "branch" = numeric(nrow(pars2)*length(branch)),
                   "equal" = character(nrow(pars2)*length(branch)),
                   stringsAsFactors = FALSE)
for(i in 1:nrow(pars2)){
  for(k in 1:length(branch)){
    # MY function
    my_function = opt_loglik(brts = branch[[k]], pars2 = pars2[i, ])
    my = apply(input, MARGIN = 1, FUN = my_function)
    my.results = c(my.results, my)
    # PBD's function
    pbd = apply(input, MARGIN = 1, FUN = pbd_loglik,
                brts = branch[[k]],
                pars2 = as.vector(pars2[i, ], mode = "character"))
    pbd.results = c(pbd.results, pbd)
    
    ind = ((i-1)*length(branch)) + k
    check[ind, 1:2] = c(i, k)
    check[ind, 3] = all.equal(my, pbd)
  }
}

check
sum(check$equal != "TRUE")









