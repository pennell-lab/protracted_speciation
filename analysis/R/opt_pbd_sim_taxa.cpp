
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
Rcpp::NumericMatrix pbdLoop_taxa(Rcpp::NumericVector pars, int taxa){
  double t = 0;
  
  int id1 = 0;
  int id = id1 + 1;
  //int Sid1 = 0;
  int Sid = 1;
  Rcpp::NumericVector sg = NumericVector::create(NA_REAL, id);
  Rcpp::NumericVector si = NumericVector::create(NA_REAL);
  Rcpp::NumericMatrix L(1, 6);
  L(0, 0) = id;
  L(0, 1) = 0;
  L(0, 2) = -1 * exp(-10);
  L(0, 3) = 0;
  L(0, 4) = -1;
  L(0, 5) = 1;
  
  
  int Ng = sg.size() - 1;
  int Ni = si.size() - 1;
  Rcpp::NumericVector probs = estimateProbs(pars, Ng, Ni);
  double denom = sum(probs);
  probs = probs / denom;
  t = t - log(Rcpp::runif(1)[0]) / denom;
  
  
  int parent, iddie, todie, idcomplete, tocomplete, event;
  Rcpp::NumericVector indices = NumericVector::create(1, 2, 3, 4, 5);
  Rcpp::NumericVector vec;
  event = Rcpp::RcppArmadillo::sample(indices, 1, FALSE, probs)[0];
  while(1) {
    if (event == 1) {
      parent = get_which(sg);
      id += 1;
      vec = NumericVector::create(id, parent, t, -1, -1, L(abs(parent) - id1 - 1, 5));
      L = addRow(L, vec);
      si.insert(si.end(), -id);
      Ni += 1;
    } else if (event == 2) {
      iddie = get_id(sg);
      todie = sg(iddie);
      L(todie - id1 - 1, 4) = t;
      sg.erase(sg.begin() + iddie);
      Ng -= 1;
    } else if (event == 3) {
      if(Ng == taxa){
        break;
      }
      idcomplete = get_id(si);
      tocomplete = abs(si(idcomplete));
      L(tocomplete - id1 - 1, 3) = t;
      Sid += 1;
      L(tocomplete - id1 - 1, 5) = Sid;
      sg.insert(sg.end(), tocomplete);
      si.erase(si.begin() + idcomplete);
      Ng += 1;
      Ni -= 1;
    } else if (event == 4) {
      parent = get_which(si);
      id += 1;
      vec = NumericVector::create(id, parent, t, -1, -1, L(abs(parent) - id1 - 1, 5));
      L = addRow(L, vec);
      si.insert(si.end(), -id);
      Ni += 1;
    } else {
      iddie = get_id(si);
      todie = abs(si[iddie]);
      L(todie - id1 - 1, 4) = t;
      si.erase(si.begin() + iddie);
      Ni -= 1;
    }
    
    probs = estimateProbs(pars, Ng, Ni);
    denom = sum(probs);
    probs = probs / denom;
    t = t - log(Rcpp::runif(1)[0]) / denom;
    event = Rcpp::RcppArmadillo::sample(indices, 1, FALSE, probs)[0];
  }
  
  L(0, 3) = t;
  
  return L;
}

