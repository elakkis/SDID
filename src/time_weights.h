#include <RcppArmadillo.h>
#include <math.h>
#include <string>
#include "admm_solver.h"

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat time_weights(arma::mat &Y00, arma::mat &Y01, double dzeta, double rho = 1, double tol = 1e-8, std::string weight_type = "SDID");
  