#include <RcppArmadillo.h>
#include <math.h>
#include <string>

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat unit_weights(arma::mat &Y00, arma::mat &Y10, double dzeta, double rho = 1, double tol = 1e-8, std::string weight_type = "SDID");
  