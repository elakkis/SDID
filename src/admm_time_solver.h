#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat admm_time_solver(arma::mat &Dmat, arma::colvec &dvec,double rho, int max_iter,double tol);