#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat att(arma::mat &X, arma::colvec &y, int unit_exclude, arma::mat &unit_weights, arma::mat &time_weights);
  