

#setwd("~/Desktop/gitHub/protracted_sp/SSP_DR/")

library(doParallel)
library(PBD)
library(ape)
library(ggplot2)
#library(gridExtra)

#source("../analysis/R/get_phylo.R")
source("get_phylo.R")
#source("../analysis/R/opt_pbd_sim_cpp.R")
source("opt_pbd_sim_cpp.R")
#source("R/auxiliary_functions.R")
source("auxiliary_functions.R")

#### pbd_sim()
#   -pars --> Vector of parameters: 
# pars[1] corresponds to b_1, the speciation-initiation rate of good species 
# pars[2] corresponds to la_1, the speciation-completion rate 
# pars[3] corresponds to b_2, the speciation-initiation rate of incipient species 
# pars[4] corresponds to mu_1, the extinction rate of good species 
# pars[5] corresponds to mu_2, the extinction rate of incipient species 
#   -age --> Sets the age for the simulation
#   -soc --> Sets whether this age is the stem (1) or crown (2) age
#   -plotit --> Sets whether the various trees produced by the function should be plotted or not

theta0 = list("b1" = c(0.3, 1, 3), #0.5
              "la1" = c(10, 5, 1, 0.3), #2,0.5
              "la2" = c(10, 5, 1, 0.3), #2,0.5
              "b2" = c(0),
              "mu1" = c(0, 0.1, 0.4, 1), #0.2
              "mu2" = c(0, 0.1, 1, 2)) #0.7

parameters = expand.grid(theta0)
parameters = parameters[which(parameters$b1 > parameters$mu1), ] # r1 positive

transition = 0.1
nphy = 10
ntips = 100
ntries = 1000
cores = 3
# VARIABLE
try = lapply(1:nrow(parameters), FUN = function(xxx){
  ppp = parameters[xxx, ]
  aa = lapply(1:nphy, function(i) {
    out = opt_pbd_sim_cpp(b1 = ppp$b1, mu1 = ppp$mu1, la1 = ppp$la1,
                          la2 = ppp$la2, b2 = ppp$b2, mu2 = ppp$mu2, 
                          trans = transition, 
                          taxa = ntips, ntry = ntries)
    return(out)})
  
  ls = lapply(aa, "[[", 1)
  full = lapply(ls, get_phylo)
  phy = lapply(full, drop.fossil)
  dr = lapply(phy, get_DR)
  ssp = mapply(get_ssp, mat = ls, tree = phy, SIMPLIFY = FALSE)
  my.data = do.call(rbind, 
                    mapply(function(x, y) 
                      data.frame(dr=x, ssp=y[ , 1], dyn=y[ , 2]),
                      x = dr, y = ssp, SIMPLIFY = F))
  my.data$dyn = as.factor(my.data$dyn)
  save(my.data, aa, ppp, file = paste0("SSP_DR_sim_trans_01_", xxx, ".RData"))
  
  return(my.data)
})#, mc.cores = cores)
save(try, "SSP_DR_sim_trans_01.RData")


### PLOT
hist_top = ggplot(my.data, aes(dr))+
  geom_histogram(aes(y = ..density..))
empty = ggplot() + geom_point(aes(1,1), colour="white") +
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter = ggplot(my.data, aes(dr, ssp))+ #, color = dyn
  geom_point(alpha = 0.5)
#scale_fill_manual(values = c("red", "blue"))
hist_right = ggplot(my.data, aes(ssp))+
  geom_histogram(aes(y = ..density..))+
  coord_flip()
gg = grid.arrange(hist_top, empty, scatter, hist_right,
                  ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
ggsave(filename = paste0("SSP_DR_sim_", paste0(ppp, collapse = "_"), ".pdf"),
       plot = gg)


# # FIXED
# try = lapply(1:3, function(i) opt_pbd_sim_cpp2(1, 0.3, 1, 100, 0, 0.5, taxa = 100, ntry = 1000))
# Ls = lapply(try, "[[", 1)
# step1 = lapply(Ls, twoL)
# source("get_phylo.R")
# step2 = lapply(step1, function(x) lapply(x, get2phylo))
# phy = mapply(combineL, phy1 = lapply(step2, "[[", 1), phy2 = lapply(step2, "[[", 2), SIMPLIFY = FALSE)
# plot(phy[[3]]);axisPhylo()


source("~/Desktop/gitHub/protracted_sp/SSP_DR/R/auxiliary_functions.R")
source("~/Desktop/gitHub/protracted_sp/analysis/R/get_phylo.R")
nn = 100
comb.data = data.frame(dr = NA, ssp = NA, dyn = NA, age = NA)
for(i in 1:nn){
  load(paste0("~/Desktop/gitHub/protracted_sp/analysis/data/SSP_DR_sim/SSP_DR_sim_trans_01_", i, ".RData"))
  this_phy = lapply(lapply(lapply(aa, "[[", 1), get_phylo), drop.fossil)
  my.data$age = unlist(lapply(this_phy, get_sp_age))
  comb.data = rbind(comb.data, my.data)
}; comb.data = comb.data[-1 , ]
ggplot(comb.data, aes(dr, ssp)) + geom_point() +
  labs(x = "Diversification rate (DR)", y = "Subspecies richness") + 
  labs(title = "Simulated data") + 
  theme(title = element_text(size = 28,face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=14), axis.title=element_text(size=22))
ggsave("protracted_sp/SSP_DR/SSP_DR_sim.png")
sim.dr.age = ggplot(comb.data, aes(age, log(dr))) + geom_point()
sim.dr.age.log = ggplot(comb.data, aes(log(age), log(dr))) + geom_point()
grid.arrange(sim.dr.age, sim.dr.age.log)
ggsave(grid.arrange(sim.dr.age, sim.dr.age.log), filename = "protracted_sp/SSP_DR/data/species_age_sim.png")


