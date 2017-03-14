
context("opt_sim_cpp")

# la1 = pars[1]
# la2 = pars[2]
# la3 = pars[3]
# mu1 = pars[4]
# mu2 = pars[5]
input = data.frame("la1" = c(0.2, 0.1, 1, 0.1),
                   "la2" = c(0.2, 0.2, 0.2, 0.1),
                   "la3" = c(0, 0, 0, 0),
                   "mu1" = c(0, 0, 0.5, 0),
                   "mu2" = c(1, 1, 0.02, 0.2))
Ntaxa = 50
Ttime = 15
Ntry = 10

### simulate for taxa 
sim.taxa = lapply(1:nrow(input), function(x) opt_pbd_sim_cpp(pars = as.numeric(input[x, ]), taxa = Ntaxa, ntry = Ntry))
sim.taxa = lapply(sim.taxa, '[[', 1)
age.taxa = sapply(sim.taxa, function(x) x[1, 3])
phy.taxa = mapply(get_phylo, L = sim.taxa, age = age.taxa, SIMPLIFY = FALSE)

### simulate for time
# CROWN age
sim.time.crown = lapply(1:nrow(input), function(x)
  opt_pbd_sim_cpp(pars = as.numeric(input[x, ]), age = Ttime, ntry = Ntry, soc = 2))
sim.time.crown = lapply(sim.time.crown, '[[', 1)
# STEM age
sim.time.stem = lapply(1:nrow(input), function(x) 
  opt_pbd_sim_cpp(pars = as.numeric(input[x, ]), age = Ttime, ntry = Ntry, soc = 1))
sim.time.stem = lapply(sim.time.stem, '[[', 1)


library(ape)
library(Rcpp)
library(RcppArmadillo)

testthat("TAXA: L matrix first column (incipient species id) is sequential",{
  incip = lapply(sim.taxa, function(x) x[ , 1])
  sequential = lapply(sim.taxa, function(x) 1:nrow(x))
  expect_equal(incip, sequential)
})

testthat("TIME: L matrix first column (incipient species id) is sequential [crown]",{
  incip = lapply(sim.time.crown, function(x) x[ , 1])
  sequential = lapply(sim.time.crown, function(x) 1:nrow(x))
  expect_equal(incip, sequential)
})

testthat("TIME: L matrix first column (incipient species id) is sequential [stem]",{
  incip = lapply(sim.time.stem, function(x) x[ , 1])
  sequential = lapply(sim.time.stem, function(x) 1:nrow(x))
  expect_equal(incip, sequential)
})


testthat("TAXA: number of tips in phylo is the same as input 'taxa'",{
  tips = sapply(lapply(phy.taxa, drop.fossil), Ntip)
  expect_equal(incip, rep(Ntaxa, times = length(tips)))
})

testthat("TIME: age of simulation is the same as input 'age' [crown]",{
  age.crown = lapply(sim.time.crown, function(x){
    age = x$L[1, 3]
    get_phylo(x$L, age)
  })
  age.crown = sapply(age.crown, function(x){
    max(branching.times(x))
  })
  exp = rep(Ttime, length(age.crown))
  expect_equal(incip, sequential)
})

testthat("TIME: age of simulation is the same as input 'age' [stem]",{
  age.stem = lapply(sim.time.stem, function(x){
    age = x$L[1, 3]
    get_phylo(x$L, age)
  })
  age.stem = sapply(age.stem, function(x){
    max(branching.times(x)) + x$root.edge
  })
  exp = rep(Ttime, length(age.stem))
  expect_equal(incip, sequential)
})
