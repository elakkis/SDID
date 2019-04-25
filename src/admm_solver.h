#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat admm_solver(arma::mat &Dmat, arma::colvec &dvec, double rho, int max_iter = 10000,double tol = 1e-8);
