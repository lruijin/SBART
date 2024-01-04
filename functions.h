#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <RcppArmadillo.h>

int sample_class(const arma::vec& probs) {
  double U = R::unif_rand();
  double foo = 0.0;
  int K = probs.size();

  // Sample
  for(int k = 0; k < K; k++) {
    foo += probs(k);
    if(U < foo) {
      return(k);
    }
  }
  return K - 1;
}

int sample_class(int n) {
  double U = R::unif_rand();
  double p = 1.0 / ((double)n);
  double foo = 0.0;

  for(int k = 0; k < n; k++) {
    foo += p;
    if(U < foo) {
      return k;
    }
  }
  return n - 1;
}

double logit(double x) {
  return log(x) - log(1.0-x);
}

double expit(double x) {
  return 1.0 / (1.0 + exp(-x));
}

// [[Rcpp::export]]
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision) {
  arma::vec z = arma::zeros<arma::vec>(mean.size());
  for(unsigned int i = 0; i < mean.size(); i++) {
    z(i) = norm_rand();
  }
  arma::mat Sigma = inv_sympd(Precision);
  arma::mat L = chol(Sigma, "lower");
  arma::vec h = mean + L * z;
  /* arma::mat R = chol(Precision); */
  /* arma::vec h = solve(R,z) + mean; */
  return h;
}

// [[Rcpp::export]]
arma::mat choll(const arma::mat& Sigma) {
  return chol(Sigma);
}

double log_sum_exp(const arma::vec& x) {
  double M = x.max();
  return M + log(sum(exp(x - M)));
}

/* Truncated normal proposal*/
// [[Rcpp::export]]
double truncated_normal_arma(double mu, double sigma, double lower, double upper) {
  double x;
  do {
    x = arma::randn() * sigma + mu; // Generate a normal random variable
  } while (x < lower || x > upper);
  return x;
}

// [[Rcpp::export]]
double log_truncated_normal_pdf(double x, double mu, double sigma, double lower, double upper) {
  double cdf_upper = arma::normcdf(upper, mu, sigma);
  double cdf_lower = arma::normcdf(lower, mu, sigma);
  return log(arma::normpdf(x, mu, sigma) / (cdf_upper - cdf_lower));
}

#endif

/* Quantile Normalize indices */
arma::vec trank(arma::vec x) {
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
arma::mat quantile_normalize_bart(arma::mat X) {
  int P = X.n_cols;
  for(int i = 0; i < P; i++){
    X.col(i) = trank(X.col(i));
  }
  
  return(X);
}

