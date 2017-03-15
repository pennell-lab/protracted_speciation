
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
gg = ggplot(birds, aes(SUBSPECIES)) + geom_histogram() + facet_grid(. ~ LAT.RANGE)
gg
ggsave("analysis/output/subspecies_per_region_Botero14.pdf", gg)
# using MLE to determine the parameter of the prior
# it is the same as estimate analytic
# LLexp <- function(rate) {
#   out = dexp(mixed$SUBSPECIES, rate, log = TRUE)
#   -sum(out)
# }
# mle(LLexp, start = list(rate = 1))

### estimate the parameter for the LOG NORMAL prior
b1names = list("Mixed", "Tropical", "Temperate", c("Mixed", "Tropical"), c("Mixed", "Temperate"))
prior4b1 = sapply(b1names, function(nm){
  LLlnorm = function(meanlog, sdlog){
    out = suppressWarnings(dlnorm(birds$SUBSPECIES[birds$LAT.RANGE %in% nm], meanlog, sdlog, log = TRUE))
    -sum(out)
  }
  par = mle(LLlnorm, start = list(meanlog = 1, sdlog = 1))
  coef(par)
})
colnames(prior4b1) = c("Mixed", "Tropical", "Temperate", "Mixed+Tropical", "Mixed+Temperate")
prior4b1
# hist(mixed$SUBSPECIES, breaks = 50, freq = FALSE)
# curve(dlnorm(x, coef(par)), 0, 40, add=TRUE)

save(birds, b1names, prior4b1, file = "analysis/data/Botero14.RData")
