

# setwd("/Users/mauro/Desktop/gitHub")

library(reshape2)
library(ape)
library(ggplot2)
library(dplyr)
library(MASS)
library(visreg)
library(phytools)



# get 1 of Jetz1s trees
trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = trees[[1]]
plot(jetz, show.tip.label = FALSE); axisPhylo()

# estimate the DR
# source("SSP_DR/get_DR.R")
# dr = get_DR(jetz)
# save(dr, file = "SSP_DR/DR_estimate.RData")
load("SSP_DR/DR_estimate.RData")
dr2plot = data.frame("DR" = dr, "log" = log(dr))
ggplot(dr2plot, aes(x = DR)) + geom_histogram(bins = 100) + xlim(0, 4.5)
ggplot(dr2plot, aes(x = log)) + geom_histogram(aes(y = ..density..), bins = 100) +
  stat_function(geom = "point", n = length(dr), fun = dnorm, args = list(mean(log(dr)), sd(log(dr))))


botero = read.csv("protracted_sp/analysis/data/Botero14.csv", as.is = TRUE)
# select only birds (ie, drop mammals)
birds = botero %>% 
  filter(TAXON == "birds") %>%
  arrange(SPECIES)
summary(birds)
# estimate the parameter of the poisson
(lambda = mean(birds$SUBSPECIES))
var(birds$SUBSPECIES) # it's completely different
# plot the histogram of #ssp and the predicted poisson distribution
ggplot(birds, aes(SUBSPECIES)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  stat_function(geom = "point", n = max(birds$SUBSPECIES), fun = dpois, args = list(lambda))

# select only the species with data
present = names(dr) %in% birds$SPECIES
dat = data.frame("sp" = names(dr)[present],
                 "DR" = dr[present],
                 "ssp" = numeric(sum(present)), 
                 "lat" = numeric(sum(present)),
                 stringsAsFactors = FALSE)
dat = dat[sort(dat$sp), ]
# there's 1 sp present in Botero dataset that doesn't appear in the phylogeny "Pseudohirundo_griseopyga"
birds = birds[-which(birds$SPECIES == "Pseudohirundo_griseopyga"), ]
dat$ssp = birds$SUBSPECIES
dat$lat = birds$CENTROID
summary(dat)
# write.csv(dat, row.names = FALSE, file = "SSP_DR/DR_ssp_centroid.csv")
dat=read.csv(file = "SSP_DR/DR_ssp_centroid.csv", row.names = NULL)
# plot DATA
ggplot(dat, aes(DR, ssp, color = lat)) + geom_point() + scale_colour_gradient2()
ggsave("SSP_DR/DR_ssp_centroid.pdf")
ggplot(dat, aes(x=lat, y=DR)) + geom_point() + geom_smooth(method = "lm") 
ggplot(dat, aes(x=lat, y=ssp)) + geom_point() + geom_smooth(method = "lm")


######## estimate the effect of #ssp on DR
# fit the model with POISSON
model = glm(ssp ~ DR, data = dat, family = poisson(link = log))
summary(model)
anova(model)
plot(model)
visreg(model)
ggplot(dat, aes(DR, ssp)) + geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) 


# fit the model with QUASIpoisson
model2 = glm(ssp ~ DR, data = dat, family = quasipoisson)
summary(model2)
anova(model2)
plot(model2)
visreg(model2)

# plot the histogram of #ssp and the predicted NEGATIVE BINOMIAL distributions
nbinom = fitdistr(birds$SUBSPECIES, densfun = "negative binomial")
my.points = data.frame("xxx" = 1:max(birds$SUBSPECIES),
                       "yyy" = dnbinom(1:max(birds$SUBSPECIES), size = nbinom$estimate[1], mu = nbinom$estimate[2]))
ggplot(birds, aes(SUBSPECIES)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  stat_function(aes(color = "Poisson"), geom = "point", n = max(birds$SUBSPECIES), fun = dpois, args = list(lambda)) + 
  # this produces NaNs....
  #stat_function(aes(color = "Negat Binomial"), geom = "point", n = max(birds$SUBSPECIES),
  #              fun = dnbinom, args = list(nbinom$estimate[1], nbinom$estimate[2])) + 
  geom_point(data = my.points, aes(x = xxx, y = yyy, color = "Negat Binomial")) +
  guides(color = guide_legend("Model Type"))
# fit the model with NEGATIVE BINOMIAL
model3 = glm.nb(ssp ~ DR, data = dat)
summary(model3)
visreg(model3)
anova(model3)
plot(model3)



all.models = data.frame("xxx" = rep(dat$DR, 3),
                        "model" = rep(c("Poisson", "QuasiPois", "NegBinom"), each = length(dat$DR)),
                        "yyy" = log(c(model$fitted.values, model2$fitted.values, model3$fitted.values)))
ggplot(dat, aes(DR, ssp)) + geom_point() + scale_y_continuous(trans = "log") +
  geom_point(data = all.models, aes(xxx, yyy, color = model)) 





#### Now considering the phylogenetic dependency of the observations
pruned = drop.tip(jetz, tip = c(jetz$tip.label[!present], "Pseudohirundo_griseopyga"))
my.pic = data.frame("DR" = pic(dat$DR, pruned),
                    "ssp" = pic(dat$ssp, pruned))
summary(my.pic)
# estimate PHYLOGENETIC SIGNAL
physig = list("dr" = list(), "ssp" = list())
tophysig = list(dat$DR, dat$ssp)
tophysig = lapply(tophysig, function(x){names(x) = dat$sp;x})
physig[[1]] = phylosig(pruned, tophysig[[1]], method = "lambda")
physig[[2]] = phylosig(pruned, tophysig[[2]], method = "lambda")
physig
# let's look at the NEW variables
pic2plot = data.frame("var" = rep(c("dr", "ssp"), each = nrow(my.pic)),
                      "value" = c(my.pic$DR, my.pic$ssp))
stats = data.frame("group" = c("dr", "ssp"),
                   "mean" = apply(my.pic, 2, mean),
                  "sd" = apply(my.pic, 2, sd))
ggplot(pic2plot, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 200) +
  xlim(-1, 1) +
  with(stats[stats$group == "dr", ],
       stat_function(data = pic2plot[pic2plot$var == "dr",], fun = dnorm, args = list(mean = mean, sd = sd))) + 
  with(stats[stats$group == "ssp", ],
       stat_function(data = pic2plot[pic2plot$var == "ssp",], fun = dnorm, args = list(mean = mean, sd = sd))) +
  facet_grid(var ~ ., scales = "free")

# The model
ggplot(my.pic, aes(DR, ssp)) + geom_point() + geom_smooth(method = "lm") +
model.pic = lm(ssp ~ DR, data = my.pic)
summary(model.pic)
anova(model.pic)




source("SSP_DR/ggplotRegression.R")
##### MODELS
ggplotRegression(model3)
ggsave(filename = "SSP_DR/model_negativeBinomial.pdf")
ggplotRegression(model.pic)
ggsave(filename = "SSP_DR/model_pic.pdf")
###################

###### DATA
ggplot(dr2plot, aes(x = log)) + geom_histogram(aes(y = ..density..), bins = 100) +
  stat_function(geom = "point", n = length(dr), fun = dnorm, args = list(mean(log(dr)), sd(log(dr))))
ggsave(filename = "SSP_DR/DR.pdf")
ggplot(birds, aes(SUBSPECIES)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  stat_function(aes(color = "Poisson"), geom = "point", n = max(birds$SUBSPECIES), fun = dpois, args = list(lambda)) + 
  geom_point(data = my.points, aes(x = xxx, y = yyy, color = "Negat Binomial")) +
  guides(color = guide_legend("Model Type"))
ggsave(filename = "SSP_DR/ssp.pdf")
ggplot(dat, aes(x = DR, y = ssp)) +
  geom_point()
ggsave(filename = "SSP_DR/ssp_DR.pdf")
ggplot(pic2plot, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 200) +
  xlim(-1, 1) +
  with(stats[stats$group == "dr", ],
       stat_function(data = pic2plot[pic2plot$var == "dr",], fun = dnorm, args = list(mean = mean, sd = sd))) + 
  with(stats[stats$group == "ssp", ],
       stat_function(data = pic2plot[pic2plot$var == "ssp",], fun = dnorm, args = list(mean = mean, sd = sd))) +
  facet_grid(var ~ ., scales = "free")
ggsave(filename = "SSP_DR/pic.pdf")
###################

###### DR vs. SSP family level
load(file = "SSP_DR/Reference_species_per_Family.RData")
dat$Family = character(nrow(dat))
for(x in refBirdLife){
  ind = na.exclude(match(x$Scientific.name, dat$sp))
  if(!is.null(ind) & length(ind) > 0){
    dat$Family[ind] = as.character(x$Family.name[which(!is.na(ind))])
  }
}
datsummary = dat %>%
  filter(Family != "") %>%
  group_by(Family) %>%
  summarise(ssp.mean = mean(ssp),
            ssp.sd = sd(ssp),
            ssp.se = sd(ssp)/sqrt(length(ssp)),
            dr.mean = mean(DR),
            dr.sd = sd(DR),
            dr.se = sd(DR)/sqrt(length(DR)))
datsummary$ssp.sd[is.na(datsummary$ssp.sd)] = 0
datsummary$ssp.se[is.na(datsummary$ssp.se)] = 0
datsummary$dr.sd[is.na(datsummary$dr.sd)] = 0
datsummary$dr.se[is.na(datsummary$dr.se)] = 0
library(ggforce)
ggplot(datsummary, aes(x = dr.mean, y = ssp.mean)) + 
  geom_errorbarh(aes(xmin = dr.mean-dr.se, xmax = dr.mean+dr.se)) +
  geom_errorbar(aes(ymin = ssp.mean-ssp.se, ymax = ssp.mean+ssp.se), width=0) +
  geom_point() +
  facet_zoom(y = ssp.mean < 5, zoom.size = 0.5)
ggsave(filename = "SSP_DR/ssp_DR_Family_SE.pdf")
ggplot(datsummary, aes(x = dr.mean, y = ssp.mean)) + 
  geom_point()
ggsave(filename = "SSP_DR/ssp_DR_Family.pdf")
# using pic as the phylogenetic mean
