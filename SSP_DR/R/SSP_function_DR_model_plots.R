


# setwd("/Users/mauro/Desktop/gitHub")

library(gridExtra)
library(ggplot2)
library(ape)
library(nlme)
library(visreg)
library(MASS)
library(dplyr)
library(reshape2)
library(ggtree)
library(ggExtra)
library(phytools)


my.theme = theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold")) 

# get 1 of Jetz1s trees
trees = read.tree("protracted_sp/analysis/data/Jetz_etal/EricsonStage2_0001_1000/AllBirdsEricson1_reduced.tre")
jetz = trees[[1]]

botero = read.csv("protracted_sp/analysis/data/Botero14.csv", as.is = TRUE)
# select only birds (ie, drop mammals)
birds = botero %>% 
  filter(TAXON == "birds") %>%
  arrange(SPECIES)
summary(birds)
load("protracted_sp/SSP_DR/data/DR_estimate.RData")
present = names(dr) %in% birds$SPECIES
dat=read.csv(file = "protracted_sp/SSP_DR/data/Updated_Botero_dat.csv", row.names = 1, as.is = TRUE)
summary(dat)
head(dat)
load(file = "protracted_sp/SSP_DR/data/pruned.RData")
my.pic = data.frame("DR" = pic(dat$DR, pruned),
                    "ssp" = pic(dat$ssp, pruned))
summary(my.pic)
load("protracted_sp/analysis/data/physig_results.RData")
physig
load("protracted_sp/SSP_DR/data/datsummary.RData")
summary(datsummary)
load("protracted_sp/SSP_DR/data/datgenus.RData")
summary(datgenus)
load("protracted_sp/SSP_DR/data/datfamily.RData")
summary(datfamily)
pic2plot = data.frame("var" = rep(c("dr", "ssp"), each = nrow(my.pic)),
                      "value" = c(my.pic$DR, my.pic$ssp))
stats = data.frame("group" = c("dr", "ssp"),
                   "mean" = apply(my.pic, 2, mean),
                   "sd" = apply(my.pic, 2, sd))



#########################################################################
########################### THE DATA ####################################
#########################################################################
# TREE and SSP
ggssp = sapply(dat$ssp, function(x)ifelse(is.na(x), 0, x))
facet_plot(ggtree(jetz), panel = 'Subspecies richness (log)', data = data.frame(sp=dat$sp, ssp=log(dat$ssp+1)), geom = geom_segment, aes(x=0, xend = ssp, y=y, yend=y)) + theme_tree2()
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_log.png")
facet_plot(ggtree(jetz), panel = 'Subspecies richness', data = data.frame(sp=dat$sp, ssp=dat$ssp), geom = geom_segment, aes(x=0, xend = ssp, y=y, yend=y)) + theme_tree2()
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP.png")
#colloting the groups
groupInfo = split(dat$sp, as.factor(dat$Order))
jpoints = groupOTU(jetz, groupInfo)
jpoints$tip.label = sapply(jpoints$tip.label, function(x)paste0(rep(".",ggssp[which(dat$sp==x)]), collapse=""))
ggtree(jpoints, layout = "fan") + geom_tiplab2()
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_fan.png")
ggtree(jpoints, aes(color = group), layout = "fan") + geom_tiplab2()
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_fan_col.png",  width = 10, height = 10, units = "in")
ggtree(jpoints, aes(color = group), layout = "fan") + geom_tiplab2() +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent")) # bg of the plot
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_fan_col_transp.png",  width = 10, height = 10, units = "in", bg = "transparent")
ggtree(jpoints, aes(color = group), layout = "fan") + 
  theme(legend.position="left", legend.title.align = 0.5) + guides(fill = guide_legend(title = "Order"))
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_col_legend.png",  width = 7, height = 7, units = "in")
# labelling the ORDER
(mrcAncestor = sapply(groupInfo, function(x) ifelse(length(x)>1, MRCA(x, obj = jetz), NA)))
ggjetz = ggtree(jetz)
aaa = na.exclude(mrcAncestor[which(names(mrcAncestor) %in% names(table(dat$Order)[which(table(dat$Order)>200)]))]) # only the groups with 100+ tips
bs = 2; fs = 5; off = rep(seq(0,16,by=4), 2); ot = 1
ggj = ggjetz; for(i in 1:length(aaa)){
  ggj = ggj + geom_cladelabel(aaa[i], label = names(aaa)[i], barsize = bs, color = "black", angle = 270, fontsize = fs, hjust = 0.5, offset.text = ot, offset = off[i])
}; ggj
# Using only the species present in BOTERO
bot = dat %>% filter(dat$sp %in% birds$SPECIES)
groupInfo2 = split(bot$sp, as.factor(bot$Order))
jetz2 = drop.tip(jetz, tip = jetz$tip.label[which(!jetz$tip.label %in% bot$sp)])
jpoints2 = groupOTU(jetz2, groupInfo2)
ggssp2 = sapply(bot$ssp, function(x)ifelse(is.na(x), 0, x))
jpoints2$tip.label = sapply(jpoints2$tip.label, function(x)paste0(rep(".",ggssp2[which(bot$sp==x)]), collapse=""))
ggtree(jpoints2, aes(color = group), layout = "fan") + geom_tiplab2(size = 3)
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_fan_col_BOTERO.png",  width = 10, height = 10, units = "in")
ggtree(jpoints2, aes(color = group), layout = "fan") + theme(legend.position="left", legend.title.align = 0.5)
ggsave(filename = "protracted_sp/SSP_DR/output/jetz_SSP_col_legend_BOTERO.png",  width = 7, height = 7, units = "in")

# triangles
backbone = data.frame(tip.label = sapply(unique(dat$Order), function(x) dat$sp[which(dat$Order == x)[1]]),
                      clade.label = unique(dat$Order),
                      N = sapply(unique(dat$Order), function(x) sum(dat$Order == x)),
                      stringsAsFactors = FALSE)
jetz.single = drop.tip(phy = jetz, tip = dat$sp[!dat$sp%in%backbone$tip.label])
backbone$depth = sapply(backbone$tip.label,function(x,y) 0.5*y$edge.length[which(y$edge[,2]==which(y$tip.label==x))],y=jetz.single)

# histogram DR
ggplot(dat, aes(x = DR)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  stat_function(geom = "point", n = length(dat$DR), fun = dnorm, args = list(mean(dat$DR, na.rm=T), sd(dat$DR, na.rm=T)))
ggsave(filename = "protracted_sp/SSP_DR/output/DR.pdf")
# plot the histogram of #ssp and the predicted NEGATIVE BINOMIAL distributions
nbinom = fitdistr(birds$SUBSPECIES, densfun = "negative binomial")
my.points = data.frame("xxx" = 1:max(birds$SUBSPECIES),
                       "yyy" = dnbinom(1:max(birds$SUBSPECIES), size = nbinom$estimate[1], mu = nbinom$estimate[2]))
# estimate the parameter of the poisson
(lambda = mean(birds$SUBSPECIES))
var(birds$SUBSPECIES) # it's completely different
ggplot(birds, aes(SUBSPECIES)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  stat_function(aes(color = "Poisson"), geom = "point", n = max(birds$SUBSPECIES), fun = dpois, args = list(lambda)) + 
  # this produces NaNs....
  #stat_function(aes(color = "Negat Binomial"), geom = "point", n = max(birds$SUBSPECIES),
  #              fun = dnbinom, args = list(nbinom$estimate[1], nbinom$estimate[2])) + 
  geom_point(data = my.points, aes(x = xxx, y = yyy, color = "Negat Binomial")) +
  guides(color = guide_legend("Model Type"))
ggsave(filename = "protracted_sp/SSP_DR/output/ssp_models.pdf")
# histogram SSP
ggplot(birds, aes(SUBSPECIES)) + geom_histogram(aes(y = ..density..), binwidth = 1)
ggsave(filename = "protracted_sp/SSP_DR/output/ssp.pdf")
# SSP as a function of DR
scatter = ggplot(dat, aes(x = exp(DR), y = ssp)) + geom_point() +
  labs(x = "Diversification rate (DR)", y = "Subspecies richness") + xlim(0,4) + #my.theme
  theme(axis.text=element_text(size=15), axis.title=element_text(size=24,face="bold")) 
ggsave(scatter, filename = "protracted_sp/SSP_DR/output/ssp_DR.png")
ggMarginal(scatter, type = "histogram", xparams = list(binwidth = .1), yparams = list(binwidth = 1)) 
ggsave(plot = ggMarginal(scatter, type = "histogram", xparams = list(binwidth = .1), yparams = list(binwidth = 1)), filename = "protracted_sp/SSP_DR/output/ssp_DR_marginal.png")
# histogram of AGE
age1log = ggplot(dat, aes(age)) + scale_x_log10() + geom_histogram(aes(y = ..density..), bins = 50)
age1 = ggplot(dat, aes(age)) + geom_histogram(aes(y = ..density..), binwidth = 1)
# DR as a function of AGE
age2log = ggplot(dat, aes(x = age, y = DR)) + scale_x_log10() + geom_smooth(method = "lm") + geom_point()
age2 = ggplot(dat, aes(x = age, y = DR)) + geom_point()
# SSP as a function of AGE
age3log = ggplot(dat, aes(x = age, y = log(ssp))) + scale_x_log10() + geom_smooth(method = "lm") + geom_point()
age3 = ggplot(dat, aes(x = age, y = log(ssp))) + geom_point()
grid.arrange(age1, age2, age3)
ggsave(grid.arrange(age1, age2, age3), filename = "protracted_sp/SSP_DR/output/species_age.pdf")
grid.arrange(age1log, age2log, age3log)
ggsave(grid.arrange(age1log, age2log, age3log), filename = "protracted_sp/SSP_DR/output/species_age_LOG.pdf")
ggsave(age2log, filename = "protracted_sp/SSP_DR/output/species_age_DR.pdf")
# PIC for ssp and DR
ggplot(pic2plot, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 200) +
  xlim(-1, 1) +
  with(stats[stats$group == "dr", ],
       stat_function(data = pic2plot[pic2plot$var == "dr",], fun = dnorm, args = list(mean = mean, sd = sd))) + 
  with(stats[stats$group == "ssp", ],
       stat_function(data = pic2plot[pic2plot$var == "ssp",], fun = dnorm, args = list(mean = mean, sd = sd))) +
  facet_grid(var ~ ., scales = "free")
# PIC - SSP as a function of DR
ggplot(pic2plot, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 200) +
  xlim(-1, 1) +
  with(stats[stats$group == "dr", ],
       stat_function(data = pic2plot[pic2plot$var == "dr",], fun = dnorm, args = list(mean = mean, sd = sd))) + 
  with(stats[stats$group == "ssp", ],
       stat_function(data = pic2plot[pic2plot$var == "ssp",], fun = dnorm, args = list(mean = mean, sd = sd))) +
  facet_grid(var ~ ., scales = "free")
ggsave(filename = "protracted_sp/SSP_DR/output/pic.pdf")
# LATITUDE
ggplot(dat, aes(log(DR), ssp, color = lat)) + geom_point() + scale_colour_gradient2()
ggsave("protracted_sp/SSP_DR/output/DR_ssp_centroid.pdf")
ggplot(dat, aes(x=lat, y=log(DR))) + geom_point() + geom_smooth(method = "lm") 
ggplot(dat, aes(x=lat, y=ssp)) + geom_point() + geom_smooth(method = "lm")
# Rates of diversification
rates = melt(datsummary, measure.vars = c(16:20))
ggplot(rates, aes(value, y = ..density..)) + geom_histogram(bins = 50) + facet_grid(variable ~ ., scales = "free")
hist.b = ggplot(datsummary, aes(b, y = ..density..)) + geom_histogram()
hist.bdLik = ggplot(datsummary, aes(bdLik, y = ..density..)) + geom_histogram()
hist.ssp.rate = ggplot(datsummary, aes(ssp.rate, y = ..density..)) + geom_histogram()
hist.ssp.lam = ggplot(datsummary, aes(ssp.lambda, y = ..density..)) + geom_histogram()
grid.arrange(hist.b, hist.bdLik)
ggsave(grid.arrange(hist.b, hist.bdLik), filename = "protracted_sp/SSP_DR/output/Divers_rate_family.pdf")
grid.arrange(hist.ssp.rate, hist.ssp.lam)
ggsave(grid.arrange(hist.ssp.rate, hist.ssp.lam), filename = "protracted_sp/SSP_DR/output/Divers_rate_SSP.pdf")
###################




#########################################################################
######################### THE MODELS ####################################
#########################################################################
######## estimate the effect of #ssp on DR
# fit the model with POISSON
model = glm(ssp ~ DR, data = dat, family = poisson(link = log))
summary(model)
anova(model)
visreg(model)
ggplot(dat, aes(DR, ssp)) + geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) 
# fit the model with QUASIpoisson
model2 = glm(ssp ~ DR, data = dat, family = quasipoisson)
summary(model2)
anova(model2)
visreg(model2)
# fit the model with NEGATIVE BINOMIAL
model3 = glm.nb(ssp ~ DR, data = dat)
summary(model3)
visreg(model3)
anova(model3)


all.models = data.frame("xxx" = rep(dat$DR, 3),
                        "model" = rep(c("Poisson", "QuasiPois", "NegBinom"), each = length(dat$DR)),
                        "yyy" = log(c(model$fitted.values, model2$fitted.values, model3$fitted.values)))
ggplot(dat, aes(DR, ssp)) + geom_point() + scale_y_continuous(trans = "log") +
  geom_point(data = all.models, aes(xxx, yyy, color = model)) 


####### SSP vs Age
mod.age.log = lm(log(ssp) ~ log(age), data = dat)
summary(mod.age.log)
anova(mod.age.log)
visreg(mod.age.log)
dat2gls = na.exclude(dat[ , c("sp", "ssp", "age")])
jetz2gls = drop.tip(jetz, tip = which(!jetz$tip.label %in% dat2gls$sp))
#mod.age.gls = gls(log(ssp) ~ log(age), data = dat2gls, correlation = corPagel(value = 1, phy = jetz2gls, fixed = FALSE))
# 0.09486835 
mod.age.gls = gls(log(ssp) ~ log(age), data = dat2gls, correlation = corPagel(value = 0.09486835, phy = jetz2gls, fixed = TRUE))
summary(mod.age.gls)
anova(mod.age.gls)
visreg(mod.age.gls)



######## PIC
ggplot(my.pic, aes(DR, ssp)) + geom_point() + geom_smooth(method = "lm")
model.pic = lm(ssp ~ DR - 1, data = my.pic)
summary(model.pic)
anova(model.pic)
visreg(model.pic)
togls = list(ssp = dat$ssp, DR = dat$DR)
togls = lapply(togls, function(x) {names(x) = dat$sp;x})
model.gls = gls(ssp ~ DR, data = togls, correlation = corPagel(1, pruned, fixed = FALSE))
summary(model.gls)
anova(model.gls)
visreg(model.gls)








#########################################################################
##################### THE MODELS and PLOTS ##############################
#########################################################################
source("protracted_sp/SSP_DR/auxiliary_functions.R")
##### MODELS
ggplotRegression(model3)
ggsave(filename = "protracted_sp/SSP_DRoutput//model_negativeBinomial.pdf")
ggplotRegression(model.pic)
ggsave(filename = "protracted_sp/SSP_DR/output/model_pic.pdf")


## FAMILY level analysis
library(ggforce)
ggplot(datsummary, aes(x = dr.mean, y = ssp.mean)) + 
  geom_errorbarh(aes(xmin = dr.mean-dr.sd, xmax = dr.mean+dr.sd)) +
  geom_errorbar(aes(ymin = ssp.mean-ssp.sd, ymax = ssp.mean+ssp.sd), width=0) +
  geom_point() #+ facet_zoom(y = ssp.mean < 5, zoom.size = 0.5)
ggsave(filename = "protracted_sp/SSP_DR/output/ssp_DR_Family_SD.pdf")
ggplot(datsummary, aes(x = dr.mean, y = ssp.mean)) + 
  geom_point() + geom_smooth(method = "lm")
mod.ssp.dr = lm(ssp.mean ~ dr.mean, data = datsummary)
summary(mod.ssp.dr)
anova(mod.ssp.dr)
ggsave(filename = "protracted_sp/SSP_DR/output/ssp_DR_Family.pdf")
ggplot(datsummary, aes(x = dr.mean, y = log(ssp.mean))) + 
  geom_point() + geom_smooth(method = "lm")
mod.ssp.dr.log = lm(log(ssp.mean) ~ dr.mean, data = datsummary)
summary(mod.ssp.dr.log)
anova(mod.ssp.dr.log)
########## Family level rates
### SUBSPECIATION
# MODEL and ssp/DR
mod.bdLik.ssp = lm(ssp.mean ~ bdLik, data = datsummary)
anova(mod.bdLik.ssp)
mod.bdLik = lm(dr.mean ~ bdLik, data = datsummary)
anova(mod.bdLik)
# PLOT
bdLik.ssp = ggplot(datsummary, aes(x = bdLik, y = ssp.mean)) + geom_point()
bdLik.dr = ggplot(datsummary, aes(x = bdLik, y = dr.mean)) + geom_point() 
grid.arrange(bdLik.ssp, bdLik.dr)
ggsave(grid.arrange(bdLik.ssp, bdLik.dr), filename = "protracted_sp/SSP_DR/output/ssp_Family_bd_DR.pdf")
bdLik.ssp + labs(x = "Diversification rate (Family level)", y = "Average subspecies richness") + ylim(0,8) + my.theme
ggsave(filename = "protracted_sp/SSP_DR/output/ssp_Family_bd.png")
#######
### SPECIATION and ssp/DR/PBD
# MODEL
mod.b.ssp = lm(ssp.mean ~ b, data = datsummary)
anova(mod.b.ssp)
mod.b = lm(dr.mean ~ b, data = datsummary)
anova(mod.b)
summary(mod.b)
mod.pbd.spe = lm(pbd.specia ~ b, data = datsummary)
anova(mod.pbd.spe)
mod.pbd.la1 = lm(pbd.la1 ~ b, data = datsummary)
anova(mod.pbd.la1)
mod.pbd.b = lm(pbd.b ~ b, data = datsummary)
anova(mod.pbd.b)
# PLOT
b.ssp = ggplot(datsummary, aes(x = b, y = ssp.mean)) + geom_point()
b.dr = ggplot(datsummary, aes(x = b, y = dr.mean)) + geom_point() + geom_smooth(method = "lm")
grid.arrange(b.ssp, b.dr)
ggsave(grid.arrange(b.ssp, b.dr), filename = "protracted_sp/SSP_DR/output/ssp_Family_b.pdf")
ggplot(datsummary, aes(x = pbd.specia, y = ssp.mean)) + geom_point()
#########
### diversification and PBD
# MODEL
mod.pbd.bdLik = lm(pbd.specia ~ bdLik, data = datsummary)
anova(mod.pbd.bdLik)
mod.pbd.la1.bdLik = lm(pbd.la1 ~ bdLik, data = datsummary)
anova(mod.pbd.la1.bdLik)
mod.pbd.b.bdLik = lm(pbd.b ~ bdLik, data = datsummary)
anova(mod.pbd.b.bdLik)
# PLOT
pbd2plot = datsummary %>% select(bdLik, pbd.specia, pbd.la1, pbd.b) %>% melt(id.vars = "bdLik")
ggplot(pbd2plot, aes(bdLik, value)) + geom_point() + geom_smooth(method = "lm") + facet_grid(variable ~ ., scales = "free")
### diversification and subspeciation rates
# MODEL
mod.phill = lm(ssp.lambda ~ bdLik, data = datsummary)
anova(mod.phill)
# PLOT
lambda.bd = ggplot(datsummary, aes(x = bdLik, y = ssp.lambda)) + geom_point()  + ylim(0,4)
#######
# speciation and subspeciation rates
# MODEL
mod.phillB = lm(ssp.lambda ~ b, data = datsummary)
summary(mod.phillB)
anova(mod.phillB)
# PLOT
lambda.b = ggplot(datsummary, aes(x = b, y = ssp.lambda)) + 
  geom_point() + geom_smooth(method = "lm") + ylim(0,4)
######
# SSP as a function of ALL metrics
comp = melt(datsummary[ , c("Family", "Nsp", "ssp.mean", "dr.mean", "age.mean", "centroid", "b", "bdLik", "ssp.lambda", "ssp.rate")],
            id.vars = c("Family", "Nsp", "ssp.mean"))
comp$value[comp$N < 3] = NA
ggplot(comp, aes(x = value, y = ssp.mean)) + lims(y = c(0,9)) +
  geom_point() + facet_grid(. ~ variable, scales = "free") + geom_smooth(method = "lm")
ggsave(filename = "protracted_sp/SSP_DR/output/ssp_Family_ALL.pdf")
####################################

############ check if richness is nested
### SSP per species within GENUS
# MODEL
ggplot(datgenus, aes(log(ssp.mean))) + geom_histogram()
ggplot(datgenus, aes(Nsp)) + geom_histogram()
cor.test(log(datgenus$ssp.mean), datgenus$Nsp, method = "kendall")
mod.ssp.gen = lm(log(ssp.mean) ~ Nsp, data = datgenus)
summary(mod.ssp.gen)
anova(mod.ssp.gen)
visreg(mod.ssp.gen)
# PLOT
ssp.sp.gen = ggplot(datgenus, aes(Nsp, ssp.mean, group = Nsp)) + geom_count() + scale_size_area() + 
  scale_size_continuous(name  ="# genus", breaks = c(1,10,100,200)) +
  labs(x = "Average species richness (LOG)", y = "Average subspecies richness (LOG)") +
  scale_y_log10() + scale_x_log10() + my.theme
ssp.sp.gen10 = ggplot(datgenus, aes(Nsp, ssp.mean, group = Nsp)) + geom_boxplot() + 
  labs(x = "Average species richness", y = "Average subspecies richness") +
  xlim(0,10.5) + ylim(0,10)
grid.arrange(ssp.sp.gen, ssp.sp.gen10)
ggsave(grid.arrange(ssp.sp.gen, ssp.sp.gen10), filename = "protracted_sp/SSP_DR/output/Nested_ssp_sp_genus.png")
ggplot(datgenus, aes(Nsp, ssp.mean, group = Nsp)) + geom_jitter(width = .1, height = .1) +
  scale_size_continuous(name  ="# genus", breaks = c(1,10,100,200)) +
  labs(x = "Average species richness (LOG)", y = "Average subspecies richness (LOG)") +
  #  labs(title = "Species richness within genus vs. Subspecies richness within species") +
  scale_y_log10() + scale_x_log10() + 
  theme(title = element_text(size = 24), axis.text=element_text(size=14), axis.title=element_text(size=22,face="bold")) 
ggsave(filename = "protracted_sp/SSP_DR/output/Nested_ssp_sp_genusFULL.png")
#######
### species per genera within FAMILY
# MODEL
ggplot(datfamily, aes(log(sp.mean))) + geom_histogram()
ggplot(datfamily, aes(Ngen)) + geom_histogram()
cor.test(log(datfamily$sp.mean), datfamily$Ngen, method = "kendall")
mod.sp.fam = lm(log(sp.mean) ~ Ngen, data = datfamily)
summary(mod.sp.fam)
anova(mod.sp.fam)
visreg(mod.sp.fam)
# PLOT
sp.gen.fam = ggplot(datfamily, aes(Ngen, sp.mean, group = Ngen)) + geom_count() + scale_x_log10() + scale_y_log10()
sp.gen.fam10 = ggplot(datfamily, aes(Ngen, sp.mean, group = Ngen)) + geom_boxplot() + xlim(0,10.5) + ylim(0,10)
grid.arrange(sp.gen.fam10, sp.gen.fam)
ggsave(grid.arrange(sp.gen.fam10, sp.gen.fam), filename = "protracted_sp/SSP_DR/output/Nested_sp_gen_fam.pdf")
########
### SSP per species within FAMILY
# MODEL
ggplot(datfamily, aes(log(ssp.mean))) + geom_histogram()
ggplot(datfamily, aes(Nsp)) + geom_histogram()
cor.test(log(datfamily$ssp.mean), datfamily$Nsp, method = "kendall")
mod.ssp.fam = lm(log(ssp.mean) ~ Nsp, data = datfamily)
summary(mod.ssp.fam)
anova(mod.ssp.fam)
visreg(mod.ssp.fam)
# PLOT
ggplot(datfamily, aes(Nsp, ssp.mean, group = Nsp)) + geom_count() + scale_x_log10()
ggplot(datfamily, aes(Nsp, ssp.mean, group = Nsp)) + geom_boxplot() + xlim(0,11)





# possible causes
other.theme = theme(axis.text=element_text(size=12), axis.title=element_text(size=.1)) 
my.lab = labs(x = "", y = "")
cause1 = b.ssp + ylim(0,8) + xlim(0.0925,0.1025) + other.theme + my.lab #labs(x = "Speciation rate (Family level)", y = "Average subspecies richness")
cause2 = bdLik.ssp + ylim(0,8) + other.theme + my.lab #labs(x = "Diversification rate (Family level)", y = "Average subspecies richness")
cause3 = lambda.b + xlim(0.0925,0.1025) + other.theme + my.lab #labs(x = "Speciation rate (Family level)", y = "Subspecies diversification rate")
cause4 = lambda.bd + other.theme + my.lab #labs(x = "Diversification rate (Family level)", y = "Subspecies diversification rate")
ggsave(filename = "protracted_sp/SSP_DR/output/possible4causes.png", plot = 
         grid.arrange(cause1, cause2, cause3, cause4, 
                      top = grid::textGrob("Possible explanations", gp=grid::gpar(fontsize=28,font=2)),
                      bottom = grid::textGrob("Speciation rate (Family)                             Diversification rate (Family)", gp=grid::gpar(fontsize=14,font="bold")),
                      left = grid::textGrob("Subspecies diversification rate          Average subspecies richness", gp=grid::gpar(fontsize=13,font="bold"), rot = 90))
)




################# MATT
summary(datsummary)
mod.full = lm(log(ssp.median) ~ abs(centroid.median) * dr.median * order, data = datsummary, na.action = "na.omit")
anova(mod.full)
summary(mod.full)
mod.add = lm(log(ssp.median) ~ abs(centroid.median) + dr.median + order, data = datsummary, na.action = "na.omit")
anova(mod.add)
summary(mod.add)







######### PBD example
set.seed(661); tree = rbdtree(.5,0,5); (tree = drop.tip(tree, c(1, 11:14, 10, 18)))
plot(tree)

tree0 = drop.tip(tree, c(2, 5:8, 10:12)); tree0$tip.label = paste0("sp.",LETTERS[4:1])
ec0 = c("black", "green", "orange", "black", "blue", "red")
tc0 = c("green", "orange", "blue", "red")

tree1 = tree; tree1$tip.label = c("sp.D", paste0("spC.",2:1), paste0("spB.",5:1), paste0("spA.",4:1))
ec1 = c("black", "green", rep("orange", 3), "black", rep("blue", 9), rep("red", 7))
tc1 = c("green", rep("orange", 2), rep("blue", 5), rep("red", 4))

tree2 = tree; tree2$tip.label = c("sp.D", paste0("spC.",2:1), paste0("spB.",5:1), paste0("spA.",4:2), "spE")
ec2 = c("black", "green", rep("orange", 3), "black", rep("blue", 9), "black", rep("red", 5), "purple")
tc2 = c("green", rep("orange", 2), rep("blue", 5), rep("red", 3), "purple")


ew = 5
lo = .15
ct = 3
NOmargin = FALSE
stip = TRUE
xxx = 5.05
png("protracted_sp/SSP_DR/output/PBD_example_3phy.png", width = 1000, height = 1200)
par(mfrow=c(3,1), mar = c(2, 0, 2, 0) + 0.1)
plot.phylo(tree0, show.tip.label = TRUE, no.margin = NOmargin, label.offset = lo, edge.width = ew, cex = ct, tip.color = tc0, edge.color = ec0)
plot.phylo(tree1, show.tip.label = stip, no.margin = NOmargin, label.offset = lo, edge.width = ew, cex = ct, tip.color = rep("white", 12), edge.color = ec1)
segments(x0 = xxx, y0 = 8.7, x1 = xxx, y1 = 12.3, lwd = ew, col = "red")
segments(x0 = xxx, y0 = 3.8, x1 = xxx, y1 = 8.2, lwd = ew, col = "blue")
segments(x0 = xxx, y0 = 1.8, x1 = xxx, y1 = 3.3, lwd = ew, col = "orange")
text(x = 5.25, y = c(10.5, 6.1, 2.57, 1.1), labels = paste0("sp.",LETTERS[1:4]), col = rev(tc0), cex = ct)
plot.phylo(tree2, show.tip.label = stip, no.margin = NOmargin, label.offset = lo, edge.width = ew, cex = ct, tip.color = rep("white", 12), edge.color = ec2)
segments(x0 = xxx, y0 = 8.7, x1 = xxx, y1 = 11.2, lwd = ew, col = "red")
segments(x0 = xxx, y0 = 3.8, x1 = xxx, y1 = 8.2, lwd = ew, col = "blue")
segments(x0 = xxx, y0 = 1.8, x1 = xxx, y1 = 3.3, lwd = ew, col = "orange")
text(x = 5.25, y = c(11.9, 9.95, 6.1, 2.57, 1.1), labels = c("sp.E", paste0("sp.",LETTERS[1:4])), col = c("purple", rev(tc0)), cex = ct)
dev.off()



