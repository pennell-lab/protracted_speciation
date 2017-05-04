
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <iostream>
#include <armadillo>

using namespace Rcpp;
using namespace std;
using namespace arma;


//ordered.pars = pars[c(1, 4, 2, 3, 5)]
// [[Rcpp::export]]
NumericVector estimateProbs(NumericVector pars, int Ng, int Ni1, int Ni2){
  NumericVector out(7);
  out[0] = pars[0] * Ng;
  out[1] = pars[1] * Ng;
  out[2] = pars[2] * Ni1;
  out[3] = pars[3] * Ni2;
  out[4] = pars[4] * (Ni1 + Ni2);
  out[5] = pars[5] * Ni1;
  out[6] = pars[5] * Ni2;
  
  return out;
}



// [[Rcpp::export]]
NumericMatrix addRow7(NumericMatrix mat, NumericVector vec){
  int dim = mat.nrow();
  NumericMatrix out(dim+1, 7);
  for(int i = 0; i < dim; i++){
    out(i, _) = mat(i, _);
  }
  out(dim, _) = vec;
  
  return(wrap(out));
}





// [[Rcpp::export]]
int get_which(NumericVector vec){
  int ll = vec.size();
  if(ll == 2){
    return vec[1];
  } else{
    NumericVector pp(ll, 1.0);
    pp[0] = 0;
    NumericVector out = Rcpp::RcppArmadillo::sample(vec, 1, FALSE, pp);
    return out[0];
  }
}



// [[Rcpp::export]]
int get_id(NumericVector vec){
  int ll = vec.size();
  if(ll == 2){
    return 1;
  } else{
    NumericVector pp(ll, 1.0);
    pp[0] = 0;
    IntegerVector ss = seq_along(vec);
    ss = Rcpp::RcppArmadillo::sample(ss, 1, FALSE, pp);
    return (ss[0] - 1);
  }
}

// [[Rcpp::export]]
int get_inverse(double a){
  int out = (int)a / 2;
  out = abs(out - 2);
  
  return(out);
}

// [[Rcpp::export]]
int check_dynamics(double dyn, double tr){
  if(Rcpp::runif(1)[0] < tr){
    return(get_inverse(dyn));
  } else{
    return(dyn);
  }
}



// [[Rcpp::export]]
Rcpp::NumericMatrix pbdLoop_taxa2Var(Rcpp::NumericVector pars, double trans, int taxa, double ntry){
  
  
  double t, denom;
  int id, Sid; //int Sid1 = 0;
  int event, Ng, Ni1, Ni2, dynamics;
  Rcpp::NumericVector sg, si1, si2, vec, probs;
  
  int parent, iddie, todie, idcomplete, tocomplete;
  Rcpp::NumericVector indices = NumericVector::create(1, 2, 3, 4, 5, 6, 7);
  Rcpp::NumericMatrix L(1, 7);
  bool bb = TRUE;
  
  while(ntry > 0){
    id = 1;
    Sid = 1;
    sg = NumericVector::create(NA_REAL, 1);
    si1 = NumericVector::create(NA_REAL);
    si2 = NumericVector::create(NA_REAL);
    Ng = sg.size() - 1;
    Ni1 = si1.size() - 1;
    Ni2 = si2.size() - 1;
    Rcpp::NumericMatrix L(1, 7);
    L(0, 0) = 1;
    L(0, 1) = 0;
    L(0, 2) = 0;
    L(0, 3) = 0;
    L(0, 4) = -1;
    L(0, 5) = 1;
    L(0, 6) = 1;
    probs = estimateProbs(pars, Ng, Ni1, Ni2);
    denom = sum(probs);
    probs = probs / denom;
    t = 0;
    
    
    while(bb) {
      probs = estimateProbs(pars, Ng, Ni1, Ni2);
      denom = sum(probs);
      probs = probs / denom;
      t = t - log(Rcpp::runif(1)[0]) / denom;
      event = Rcpp::RcppArmadillo::sample(indices, 1, FALSE, probs)[0];
      switch(event){
        case (1): { 
          parent = get_id(sg);
          id += 1;
          dynamics = check_dynamics(L(parent-1, 6), trans);
          vec = NumericVector::create(id, sg(parent), t, -1, -1, L(abs(sg(parent)) - 1, 5), dynamics);
          L = addRow7(L, vec);
          if(dynamics == 1){
            si1.insert(si1.end(), -id);
            Ni1 += 1;
          } else {
            si2.insert(si2.end(), -id);
            Ni2 += 1;
          }
          break;
        } 
        case (2): {
          iddie = get_id(sg);
          todie = sg(iddie);
          L(todie - 1, 4) = t;
          sg.erase(sg.begin() + iddie);
          Ng -= 1;
          break;
        } 
        case (3): {
          if(Ng == taxa){
            bb = FALSE;
            break;
          }
          idcomplete = get_id(si1);
          tocomplete = abs(si1(idcomplete));
          L(tocomplete - 1, 3) = t;
          Sid += 1;
          L(tocomplete - 1, 5) = Sid;
          sg.insert(sg.end(), tocomplete);
          si1.erase(si1.begin() + idcomplete);
          Ng += 1;
          Ni1 -= 1;
          break;
        }
        case (4): {
          if(Ng == taxa){
            bb = FALSE;
            break;
          }
          idcomplete = get_id(si2);
          tocomplete = abs(si2(idcomplete));
          L(tocomplete - 1, 3) = t;
          Sid += 1;
          L(tocomplete - 1, 5) = Sid;
          sg.insert(sg.end(), tocomplete);
          si2.erase(si2.begin() + idcomplete);
          Ng += 1;
          Ni2 -= 1;
          break;
        }
        case (5): {
          //parent = get_which(si);
          id += 1;
          vec = NumericVector::create(id, parent, t, -1, -1, L(abs(parent) - 1, 5), NA_REAL);
          L = addRow7(L, vec);
          //si.insert(si.end(), -id);
          //Ni += 1;
          break;
        }
        case (6): {
          iddie = get_id(si1);
          todie = abs(si1[iddie]);
          L(todie - 1, 4) = t;
          si1.erase(si1.begin() + iddie);
          Ni1 -= 1;
          break;
        }
        case (7): {
          iddie = get_id(si2);
          todie = abs(si2[iddie]);
          L(todie - 1, 4) = t;
          si2.erase(si2.begin() + iddie);
          Ni2 -= 1;
          break;
        }
      }
      
      if((Ng + Ni1 + Ni2) == 0){
        bb = FALSE; // if all sp die
      }
    }
    
    ntry -= 1;
    bb = TRUE; 
    
    if(Ng == taxa){
      L(0, 3) = t;
      return L;
    } 
    
  }
  
  return L;
}




