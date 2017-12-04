
library(dplyr)
library(ggplot2)
library(stats4)



# load subspecies data
botero = read.csv("analysis/data/Botero14.csv")
head(botero)
# select only birds (ie, drop mammals)
birds = botero %>% 
  filter(TAXON == "birds")
# estimate the number of ssp for different regions of the World
ssp = birds %>%
  group_by(LAT.RANGE) %>% 
  summarise(avg = mean(SUBSPECIES, na.rm = TRUE),
            sd = sd(SUBSPECIES, na.rm = TRUE),
            median = median(SUBSPECIES, na.rm = TRUE))
# get separate data frames for each region
mixed = birds %>%
  filter(LAT.RANGE == "Mixed") %>%
  select(SPECIES,SUBSPECIES)
tropical = birds %>%
  filter(LAT.RANGE == "Tropical") %>%
  select(SPECIES,SUBSPECIES)
temperate = birds %>%
  filter(LAT.RANGE == "Temperate") %>%
  select(SPECIES,SUBSPECIES)
#levels(birds$LAT.RANGE) = paste(levels(birds$LAT.RANGE), table(birds$LAT.RANGE))
gg = ggplot(birds, aes(SUBSPECIES)) + geom_histogram() + facet_grid(. ~ LAT.RANGE) 
gg
ggsave("analysis/output/subspecies_per_region_Botero14.pdf", gg)

# using MLE to determine the parameter of the prior
# it is the same as estimate analytic
### estimate the parameter for the LOG NORMAL prior
b1names = list("Mixed", "Tropical", "Temperate", c("Mixed", "Tropical"), c("Mixed", "Temperate"))
prior4b1.lnorm = sapply(b1names, function(nm){
  LLlnorm = function(meanlog, sdlog){
    out = suppressWarnings(dlnorm(birds$SUBSPECIES[birds$LAT.RANGE %in% nm], meanlog, sdlog, log = TRUE))
    -sum(out)
  }
  par = mle(LLlnorm, start = list(meanlog = 1, sdlog = 1))
  coef(par)
})
colnames(prior4b1.lnorm) = sapply(b1names, function(nm)paste0(nm, collapse = "+"))
prior4b1.lnorm
### estimate the parameter for the EXPONENCIAL prior
prior4b1.exp = sapply(b1names, function(nm){
  LLexp = function(r) {
    out = suppressWarnings(dexp(x = birds$SUBSPECIES[birds$LAT.RANGE %in% nm], rate = r, log = TRUE))
    -sum(out)
  }
  par = mle(LLexp, start = list(r = .1))
  out = structure(coef(par), names = paste0(nm, collapse = "+"))
})
prior4b1.exp
# hist(mixed$SUBSPECIES, breaks = 50, freq = FALSE)
# curve(dlnorm(x, coef(par)), 0, 40, add=TRUE)

save(birds, b1names, prior4b1.exp, prior4b1.lnorm, file = "analysis/data/Botero14.RData")


par2plot = data.frame(region = names(prior4b1.exp), rate = prior4b1.exp, t(prior4b1.lnorm), row.names = NULL)
iii = 1
hist(birds$SUBSPECIES[which(birds$LAT.RANGE %in% b1names[[iii]])], breaks = 60, probability = TRUE, col = "black")
curve(dexp(x, rate = par2plot$rate[iii], log = FALSE), add = TRUE, col = "red", lwd = 1.5)
curve(dlnorm(x, meanlog = par2plot$meanlog[iii], sdlog = par2plot$sdlog[iii], log = FALSE), add = TRUE, col = "green", lwd = 1.5)


#LatHarsh = cor(abs(birds$CENTROID), log(birds$ENV.HARSHNESS - min(birds$ENV.HARSHNESS) + 1) )
LatHarsh = cor(abs(birds$CENTROID), birds$ENV.HARSHNESS)
#ggLatHarsh = ggplot(birds, aes(abs(CENTROID), log(ENV.HARSHNESS - min(ENV.HARSHNESS) + 1))) +
ggLatHarsh = ggplot(birds, aes(abs(CENTROID), ENV.HARSHNESS)) +
  geom_point() +
  #geom_text(aes(x = 50, y = -1.5, label = paste0("cor = ", round(LatHarsh, 3))), size = 10) +
  xlab("abs(Latitude)") + ylab("Harshness") +
  theme(axis.title.y = element_text(size=30, face="bold"),
        axis.title.x = element_text(size=30, face="bold"))
ggLatHarsh
ggsave(filename = "analysis/output/Botero14_Latitude_Harshness.pdf", plot = ggLatHarsh+geom_text(aes(x = 50, y = -1.5, label = paste0("cor = ", round(LatHarsh, 3))), size = 10) )
ggsave(filename = "analysis/output/Botero14_Latitude_Harshness.png", plot = ggLatHarsh+geom_text(aes(x = 50, y = -1.5, label = paste0("cor = ", round(LatHarsh, 3))), size = 10) )
ggsave(filename = "analysis/output/Botero14_Latitude_Harshness_NOcorr.png", plot = ggLatHarsh)
