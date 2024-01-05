// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec test(arma::vec x) {
    arma::vec x_unique = unique(x);
    int J = x_unique.size();
    arma::vec x_ranks = arma::regspace(0, J-1);
    int I = x.size();
    for(int i = 0; i < I; i++){
        arma::uvec ind = find(x_unique == x(i));
        x(i) = x_ranks(ind(0)) / J;
    }
    x = x / max(x);
  return x;
}

// [[Rcpp::export]]
arma::mat quantile_normalize_test(arma::mat X) {
  int P = X.n_cols;
  for(int i = 0; i < P; i++){
    X.col(i) = test(X.col(i));
  }
  
  return(X);
}
