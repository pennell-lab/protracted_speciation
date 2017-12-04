

# PLOTS
# Subclade clean harshness 3 parameters analysis

library(parallel)
library(geiger)
library(coda)
library(dplyr)
library(ggplot2)

# setwd("~/Desktop/gitHub/protracted_sp/")
# source("analysis/R/auxiliary_functions.R") 
file_path = "analysis/data/subclade_harshness/" # path to the files

# name of the files
files = list.files(path = file_path, pattern = "*_PBD_Bayes.txt")
files = files[grepl("_harshnessClean3par_", files)]
# creates a data frame with the info about the runs - eg, category and chain.
cntrl = do.call(rbind, sapply(files, strsplit, split = "_"))
cntrl = data.frame(type = cntrl[ , 6],
                   phy = cntrl[ , 5],
                   category = cntrl[ , 7],
                   est = cntrl[ , 9],
                   chain = cntrl[ , 11],
                   row.names = NULL)

# Ncores = 3
# # read all the files, already skips the burn-in to accelarate the process
# raw = mclapply(paste0(file_path, files), FUN = read.csv, sep = "\t", skip = ceiling(5000000 * 0.3), mc.cores = Ncores)
# # assigns the correct colnames
# raw = lapply(raw, function(x){colnames(x) = c("logLik", "prior", "posterior", "b", "mu1", "la1", "mu2", "accepted", "proposed.par", "converged"); x})
# # discards the in between values; ie, keeps only the real samples from posterior
# coda.reads = lapply(lapply(raw, function(x)x[ , 1:7]), mcmc, thin = 1000)
# # creates clean data frames, with only the real samples and the 4 parameters!
# ind = seq(1, 3500001, by = 1000)
# reads = do.call(rbind, lapply(raw, function(x)x[ind, 1:7]))
# nr = length(ind)
# reads$type = rep(cntrl$type, each = nr)
# reads$phy = rep(cntrl$phy, each = nr)
# reads$category = rep(cntrl$category, each = nr)
# reads$est = rep(cntrl$est, each = nr)
# reads$chain = rep(cntrl$chain, each = nr)
# 
# save(coda.reads, file = "analysis/data/subclade_harshness/subclade_harshnessClean3par.RData")
# write.csv(file = "analysis/data/subclade_harshness/subclade_harshnessClean3par_reads.csv", reads, row.names = FALSE)

load("analysis/data/subclade_harshness/subclade_harshnessClean3par.RData")
reads = read.csv(file = "analysis/data/subclade_harshness/subclade_harshnessClean3par_reads.csv", header = TRUE)


# general plot parameters
lab.size = 30
lab.face = "bold"
tick.text.size = 19
# gets the colours for the plot, mimeking the ggplot default
# gets the hue
h = seq(15, 375, length = 2 + 1) # number of categories plus 1
col.values = hcl(h = h, l = 65, c = 100)[2:1]
# subspecies formation rate
ggb = ggplot(reads, aes(b, fill = category)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  xlab("Rate") + ylab("Density") + ggtitle("Speciation-initiation") +
  labs(fill = "Harshness") +
  guides(fill=FALSE) + # remove legends
  scale_fill_manual(values = col.values) +
  theme(plot.title = element_text(size=lab.size, face=lab.face),
        axis.title.y = element_text(size=lab.size, face=lab.face),
        axis.title.x = element_text(size=lab.size, face=lab.face),
        axis.text = element_text(size=tick.text.size))
ggb
ggsave(filename = "analysis/output/Botero14_HarshnessClean3par_b.png", plot = ggb)
# time to speciation
ggla1 = ggplot(reads, aes(la1, fill = category)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  xlab("Rate") + ylab("") + ggtitle("Speciation-completion") +
  labs(fill = "Harshness") +
  guides(fill=FALSE) + # remove legends
  scale_fill_manual(values = col.values) +
  xlim(0, 10) +
  theme(plot.title = element_text(size=lab.size, face=lab.face),
        axis.title.y = element_text(size=lab.size, face=lab.face),
        axis.title.x = element_text(size=lab.size, face=lab.face),
        axis.text = element_text(size=tick.text.size))
ggla1
ggsave(filename = "analysis/output/Botero14_HarshnessClean3par_la1.png", plot = ggla1)
# species extinction rate
ggmu1 = ggplot(reads, aes(mu1, fill = category)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col.values)
ggmu1
# subspecies extinction rate
ggmu2 = ggplot(reads, aes(mu2, fill = category)) + 
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col.values)
ggmu2

# species diversification rate
ggSPdiv = ggplot(reads, aes(b - mu1, fill = category)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  xlab("Rate") + ylab("Density") + ggtitle("Species diversification") +
  labs(fill = "Harshness") +
  guides(fill=FALSE) + # remove legends
  scale_fill_manual(values = col.values) +
  theme(plot.title = element_text(size=lab.size, face=lab.face),
        axis.title.y = element_text(size=lab.size, face=lab.face),
        axis.title.x = element_text(size=lab.size, face=lab.face),
        axis.text = element_text(size=tick.text.size))
ggSPdiv
ggsave(filename = "analysis/output/Botero14_HarshnessClean3par_SPdivers.png", plot = ggSPdiv)
# subspecies diversification rate
ggSUBdiv = ggplot(reads, aes((1/la1) - mu2, fill = category)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  xlab("Rate") + ylab("") + ggtitle("SUBspecies diversification") +
  labs(fill = "Harshness") +
  guides(fill=FALSE) + # remove legends
  scale_fill_manual(values = col.values) +
  theme(plot.title = element_text(size=lab.size, face=lab.face),
        axis.title.y = element_text(size=lab.size, face=lab.face),
        axis.title.x = element_text(size=lab.size, face=lab.face),
        axis.text = element_text(size=tick.text.size))
ggSUBdiv
ggsave(filename = "analysis/output/Botero14_HarshnessClean3par_SUBdivers.png", plot = ggSUBdiv)


# Just to get the legend in a separate plot
legend = get_legend(ggmu1 + labs(fill = "Harshness") )
plot_grid(legend)
ggsave(filename = "analysis/output/Botero14_HarshnessClean3par_legend.png")
