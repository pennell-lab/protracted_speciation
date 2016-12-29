###########################################
#   Fitting the protracted sp model to the
# bird tree of Jetz etal. 2014 - Hackett 
# backbone. Here are 1,000 trees of a total
# of 10k available.
#
# http://birdtree.org/downloads/
###########################################

library(ape)
jetz = read.tree("analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1.tre")

class(jetz)
summary(jetz)


jetz_branch = lapply(jetz, branching.times)
jetz_age = sapply(jetz_branch, max)/5


library(PBD)
jetz_ML = lapply(jetz_branch, pbd_ML)

ML2sim = function(vector, b2 = 0){
  vector = as.numeric(vector)
  out = c(vector[1], vector[3], b2, vector[2], vector[4])
  return(out)
}
pars2sim = lapply(jetz_ML, FUN = ML2sim)

sim = list()
for(i in 1:length(jetz_ML)){
  numbsim = 100
  repeat{
    aux = lapply(rep(2, numbsim), FUN = pbd_sim, pars = pars2sim[[i]],
                 age = jetz_age[i], plotit = FALSE)
    ind = sapply(aux, FUN = function(x){ 
      if(class(x[[1]]) == "phylo"){
        return(Ntip(x[[1]]) >= 3)
      } else return(FALSE)
      })
    simi = aux[ind]
    numbsim = numbsim - sum(ind)
    if(numbsim <= 0){
      sim[[i]] = simi
      break
    }
  }
}


save(jetz_ML, jetz_age, jetz_branch, sim, file = "analysis/data/jetz_simulated_data.RData")



library(parallel)
parallel::detectCores()

first = pbd_ML(jetz_branch)
sim0 = pbd_sim(pars = ML2sim(first), age = jetz_age[1], soc = 2, plotit = FALSE)




