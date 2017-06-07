


path2read = "~/Desktop/gitHub/protracted_sp/analysis/R/"
path2data = "~/Desktop/gitHub/protracted_sp/analysis/data/"
path2plot = "~/Desktop/gitHub/protracted_sp/analysis/output/"
source(paste0(path2read, "plot_functions.R"))


par = get_param(output = path2read, out_post = "posterior/", initial.patt = "Bayes_simulation_study_")
head(par)
str(par)
#res = read.BPBD(chains = 2, folder = "~/Desktop/gitHub/protracted_sp/analysis/R/posterior/", par = par)



res = read.BPBD(chains = 2, output = path2read, out_post = "posterior/", initial.patt = "Bayes_simulation_study_")
save(res, file = paste0(path2data, "simulation_study_results.RData"))
# check ACCEPTANCE RATIO and convergence
accp = check_MCMC(res$samples, "acc")
gg0 = ggplot(melt(accp), aes(value)) + geom_histogram() + facet_grid(. ~ variable) + geom_vline(aes(xintercept = 0.3), colour = "red")
ggsave(filename = paste0(path2plot, "Simulation_study_accpRatio.pdf"), plot = gg0)

(ZERO = check_MCMC(res$samples, "all"))

# checking the "my.mcmc" object
plot(res$mcmc[[1]])
# check CONVERGENCE
library(coda)
gelman = lapply(res$mcmc, gelman.diag)
conver = do.call(rbind, lapply(gelman, function(x) setNames(c(x[[2]], x[[1]][ , 1]), c("mult","b","mu1","la1","mu2") )))
gg = ggplot(melt(conver), aes(value)) + geom_histogram() + facet_grid(. ~ Var2) + geom_vline(aes(xintercept = 1.1), colour = "red")
ggsave(filename = paste0(path2plot, "Simulation_study_convergence.pdf"), plot = gg)

ml = mclapply(X = branches, FUN = function(b) pbd_ML(brts = b, initparsopt = c(1,0.5,3,2), exteq = 0, verbose = FALSE), mc.cores = 3)
orig.par = par()
pdf(paste0(path2plot, "Simulation_study.pdf"))
for(w in 1:length(res$samples)){
  this_par = as.numeric(strsplit(names(res$samples)[w], split = "_")[[1]])
  layout(matrix(1:2, ncol = 2), widths = c(3,1))
  par(mai=(c(.5,.8,.8,0)))
  boxplot(res$samples[[w]][ , c("b", "mu1", "mu2")])#, log = "y")
  points(x = 1:3, y = this_par[c(1,4,5)], col = "red")
  points(x = 1:3, y = ml[[w]][c(1,2,4)], col = "green")
  title(main = paste("b =", this_par[1], "mu1 =", this_par[4], "mu2 =", this_par[5], "la1 =", this_par[2], "phylog", this_par[6]))
  par(mai=(c(.5,0,.8,.8)))
  boxplot(res$samples[[w]][ , "la1"], ylab="", yaxt="n", xlab="la1")#, log = "y")
  axis(4)
  points(x = 1, y = this_par[2], col = "red")
  points(x = 1, y = ml[[w]][3], col = "green")
}
dev.off()
par(orig.par)






all_par = data.frame(t(sapply(names(res$samples), function(x)as.numeric(strsplit(x, split = "_")[[1]]))))
colnames(all_par) = c("b1", "la1", "b2", "mu1", "mu2", "phylogeny")
