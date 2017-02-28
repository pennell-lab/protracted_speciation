
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


