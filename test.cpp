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

double newtonsMethod(double a, double b, double c, double d, double lowerbound, double initialGuess) {
  double x0 = initialGuess;
  double x1;
  double tolerance = 1e-6;
  int maxIterations = 1000;
  
  for (int i = 0; i < maxIterations; ++i) {
    double f = cubicFunction(a, b, c, d, x0);
    double fPrime = cubicDerivative(a, b, c, x0);
    
    if (fabs(fPrime) < std::numeric_limits<double>::epsilon()) {
      Rcpp::Rcout << "Derivative too small." << std::endl;
      return x0;
    }
    
    x1 = x0 - f / fPrime;
    
    if(x1 < lowerbound){
      return lowerbound;
    }
    if (fabs(x1 - x0) <= tolerance) {
      return x1;
    }
    
    x0 = x1;
  }
  
  Rcpp::Rcout << "Max iterations reached." << std::endl;
  return x0;
}
double cubicFunction(double a, double b, double c, double d, double x) {
  return a*x*x*x + b*x*x + c*x + d;
}

double cubicDerivative(double a, double b, double c, double x) {
  return 3*a*x*x + 2*b*x + c;
}

arma::vec find_theta_center(const std::vector<Node*>& forest, const arma::vec& Y, 
                            const arma::vec& weights, const arma::mat& X, Hypers& hypers,
                            arma::vec theta_center){
  arma::vec theta_new = theta_center;
  arma::vec theta_grid = arma::regspace(0, M_PI, 0.05);
  double llik_old = loglik_Theta(forest, theta_center, X, Y, weights, hypers);
  double llik = llik_old;
  double theta_tmp = 0.0;
  for(unsigned int i = 0; i < theta_new.size(); i++){
    theta_tmp = theta_center(i);
    // Rcpp::Rcout << "theta center is used to be" << theta_center(i);
    for(unsigned int j = 0; j < theta_grid.size(); j++){
      theta_new(i) = theta_grid(j);
      llik = loglik_Theta(forest, theta_new, X, Y, weights, hypers);
      if(llik > llik_old){
        llik_old = llik;
        theta_center(i) = theta_new(i);
      }
      // Rcpp::Rcout << "theta center is changed to" << theta_center(i);
    }
    // theta_new(i) = theta_tmp;
    //If independent from others and update one theta when others all remain.
    theta_new(i) = theta_center(i);
  }
  
  return theta_center;
}

arma::vec theta_proposal_beta(arma::vec theta, arma::vec d, const double& c){
  arma::vec theta_new = d;
  arma::vec randnum;
  double x,y;
  for(unsigned int j = 0; j < theta_new.size(); j++){
    randnum = arma::randg(1, arma::distr_param(c, 1.0));
    x = randnum(0);
    randnum = arma::randg(1, arma::distr_param(d(j),1.0));
    y = randnum(0);
    theta_new(j) = x / (x + y) * M_PI;
  }
  return(theta_new);
}

arma::vec find_beta_d(arma::vec theta_center, const double& c){
  arma::vec d = theta_center;
  for(unsigned int i = 0; i < theta_center.size(); i++){
    d(i) = c/theta_center(i) * M_PI - c;
  }
  return d;
}