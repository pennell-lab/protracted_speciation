
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
NumericVector estimateProbs(NumericVector pars, int Ng, int Ni){
  NumericVector out(5);
  out[0] = pars[0] * Ng;
  out[1] = pars[1] * Ng;
  out[2] = pars[2] * Ni;
  out[3] = pars[3] * Ni;
  out[4] = pars[4] * Ni;
  
  return out;
}



// [[Rcpp::export]]
NumericMatrix addRow(NumericMatrix mat, NumericVector vec){
  int dim = mat.nrow();
  NumericMatrix out(dim+1, 6);
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
Rcpp::NumericMatrix pbdLoop_taxa(Rcpp::NumericVector pars, int taxa, double ntry){
  
  
  double t, denom;
  int id, Sid; //int Sid1 = 0;
  int event, Ng, Ni;
  Rcpp::NumericVector sg, si, vec, probs;
  
  int parent, iddie, todie, idcomplete, tocomplete;
  Rcpp::NumericVector indices = NumericVector::create(1, 2, 3, 4, 5);
  Rcpp::NumericMatrix L(1, 6);
  bool bb = TRUE;
  
  while(ntry > 0){
    id = 1;
    Sid = 1;
    sg = NumericVector::create(NA_REAL, id);
    si = NumericVector::create(NA_REAL);
    Ng = sg.size() - 1;
    Ni = si.size() - 1;
    Rcpp::NumericMatrix L(1, 6);
    L(0, 0) = id;
    L(0, 1) = 0;
    L(0, 2) = 0;
    L(0, 3) = 0;
    L(0, 4) = -1;
    L(0, 5) = 1;
    probs = estimateProbs(pars, Ng, Ni);
    denom = sum(probs);
    probs = probs / denom;
    t = 0;
    
    
    
    while(bb) {
      probs = estimateProbs(pars, Ng, Ni);
      denom = sum(probs);
      probs = probs / denom;
      t = t - log(Rcpp::runif(1)[0]) / denom;
      event = Rcpp::RcppArmadillo::sample(indices, 1, FALSE, probs)[0];
      
      switch(event){
      case (1): { 
        parent = get_which(sg);
        id += 1;
        vec = NumericVector::create(id, parent, t, -1, -1, L(abs(parent) - 1, 5));
        L = addRow(L, vec);
        si.insert(si.end(), -id);
        Ni += 1;
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
        idcomplete = get_id(si);
        tocomplete = abs(si(idcomplete));
        L(tocomplete - 1, 3) = t;
        Sid += 1;
        L(tocomplete - 1, 5) = Sid;
        sg.insert(sg.end(), tocomplete);
        si.erase(si.begin() + idcomplete);
        Ng += 1;
        Ni -= 1;
        break;
      } 
      case (4): {
        parent = get_which(si);
        id += 1;
        vec = NumericVector::create(id, parent, t, -1, -1, L(abs(parent) - 1, 5));
        L = addRow(L, vec);
        si.insert(si.end(), -id);
        Ni += 1;
        break;
      }
      case (5): {
        iddie = get_id(si);
        todie = abs(si[iddie]);
        L(todie - 1, 4) = t;
        si.erase(si.begin() + iddie);
        Ni -= 1;
        break;
      }
      }
      
      if((Ng + Ni) == 0){
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




