
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
  for(int i = 0; i < 2; i++){
    out[i] = pars[i] * Ng;
  }
  for(int i = 2; i < 5; i++){
    out[i] = pars[i] * Ni;
  }
  
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
    for(int i = 0; i < ll; i++){
      vec[i] = i;
    }
    NumericVector out = Rcpp::RcppArmadillo::sample(vec, 1, FALSE, pp);
    return out[0];
  }
}





// [[Rcpp::export]]
Rcpp::NumericMatrix pbdLoop(Rcpp::NumericVector pars, double age){
  NumericVector t = 0;
  
  int id1 = 0;
  int id = id1 + 1;
  int Sid1 = 0;
  int Sid = 1;
  Rcpp::NumericVector sg = NumericVector::create(NA_REAL, id);
  Rcpp::NumericVector si = NumericVector::create(NA_REAL);
  Rcpp::NumericMatrix L(1, 6);
  L(0, 0) = id;
  L(0, 1) = 0;
  L(0, 2) = -1 * exp(-10);
  L(0, 3) = t[0];
  L(0, 4) = -1;
  L(0, 5) = 1;
  
  
  int Ng = sg.size() - 1;
  int Ni = si.size() - 1;
  Rcpp::NumericVector probs = estimateProbs(pars, Ng, Ni);
  double denom = sum(probs);
  probs = probs / denom;
  t = t[0] - log(Rcpp::runif(1)) / denom;
  Rcpp::NumericVector indices = NumericVector::create(1, 2, 3, 4, 5);
  
  
  while(t[0] <= age) {
    Rcout <<
      "probs = " << probs << '\n';
    
    
    Rcpp::NumericVector event = Rcpp::RcppArmadillo::sample(indices, 1, FALSE, probs);
    Rcout <<
      "time --> " << t[0] << "\t event = " << event[0] << '\n';
    if (event[0] == 1) {
      double parent = get_which(sg);
      id += 1;
      Rcpp::NumericVector vec = NumericVector::create(id, parent, t[0], -1, -1, L(abs(parent) - id1, 6));
      L = addRow(L, vec);
      si.insert(si.end(), -id);
      Ni += 1;
    } else if (event[0] == 2) {
      int iddie = get_id(sg);
      int todie = sg(iddie);
      L(todie - id1, 5) = t[0];
      sg.erase(sg.begin() + iddie);
      Ng -= 1;
    } else if (event[0] == 3) {
      int idcomplete = get_id(si);
      int tocomplete = abs(si(idcomplete));
      L(tocomplete - id1, 4) = t[0];
      Sid += 1;
      L(tocomplete - id1, 6) = Sid;
      sg.insert(sg.end(), tocomplete);
      si.erase(si.begin() + idcomplete);
      Ng += 1;
      Ni -= 1;
    } else if (event[0] == 4) {
      int parent = get_which(si);
      id += 1;
      Rcpp::NumericVector vec = NumericVector::create(id, parent, t[0], -1, -1, L(abs(parent) - id1, 6));
      L = addRow(L, vec);
      si.insert(si.end(), -id);
      Ni += 1;
    } else {
      int iddie = get_id(si);
      int todie = abs(si[iddie]);
      L(todie - id1, 5) = t[0];
      si.erase(si.begin() + iddie);
      Ni -= 1;
    }
    
    probs = estimateProbs(pars, Ng, Ni);
    denom = sum(probs);
    probs = probs / denom;
    t = t[0] - log(Rcpp::runif(1)) / denom;
  }
  
  return L;
}

