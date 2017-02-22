library(ape)
library(PBD)
source("analysis/R/opt_pbd_sim.cpp")

library(Rcpp)
library('inline')


pars = c(0.3, 2, 0, 0.1, 0.1)
age = 10
soc = 2
set.seed(666)
(mm = opt_pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 5, soc = 2))
set.seed(660)
pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 3, soc = 1, plotit = FALSE)

library(microbenchmark)
compare <- microbenchmark(opt_pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 5, soc = 1),
                          pbd_sim(pars = c(0.3, 2, 0, 0.1, 0.1), age = 3, soc = 1, plotit = FALSE),
                          times = 1000)




opt_pbd_sim = function (pars, age, soc = 2, ntry = 1) 
  # ntry - controls the number of times the simulation will try to generate a valid phylogeny
{
  # la1 = pars[1]
  # la2 = pars[2]
  # la3 = pars[3]
  # mu1 = pars[4]
  # mu2 = pars[5]
  # need to be ordered accordingly 
  pars = pars[c(1, 4, 2, 3, 5)]
  
  if(soc == 1){
    L = pbd_loop()
  } else{
    L
    
    else if(i == 2){
      int id = id1 + 1;
      int Sid = Sid1;
      Rcpp::NumericVector sg = NumericVector::create(NA_REAL);
      Rcpp::NumericVector si = NumericVector::create(-id);
      Rcpp::NumericMatrix L(1, 6);
      L(0, 0) = id;
      L(0, 1) = 1;
      L(0, 2) = t[0];
      L(0, 3) = -1;
      L(0, 4) = -1;
      L(0, 5) = 1;
    }
  }
  
  
  
  
  
  if (i == 1) {
    if ((control$Ng + control$Ni) > 0) {
      i = i + 1
      L1 = control$L
      control$id1 = control$id
      control$Sid1 = control$Sid
      control$si1 = control$si
      control$sg1 = control$sg
    }
  } else { # if (i == 2)
    if (PBD::checkgood(control$L, control$si, control$sg) == 1) {
      i = i + 1
      L2 = control$L
      si2 = control$si
      sg2 = control$sg
    }
  }
  
  if(i > soc){
    control$L = L1
    if (soc == 2) {
      control$L = rbind(L1, L2)
    }
    L0 = control$L
    absL = control$L
    absL[, 2] = abs(control$L[, 2])
    
    rtree = try(read.tree(text = PBD::detphy(absL, age)))
    if(is.null(rtree) | class(rtree) == "try-error"){
      if(!is.binary.phylo(rtree)){
        ntry = ntry-1
        if(ntry > 0){
          i = 1
        } else{
          return(NULL)
        }
      }
    } else {
      break()
    }
  }
}

tree = as.phylo(rtree)
control$L[, 3:5][which(control$L[, 3:5] == -1)] = age + 1
control$L[, 3:5] = age - control$L[, 3:5]
control$L = control$L[order(control$L[, 1]), ]
Ltreeslist = list(tree = tree, 
                  L = control$L)
return(Ltreeslist)
}
