#include <RcppArmadillo.h>
#include <math.h>
#include <string>
#include "admm_solver.h"
#include "unit_weights.h"
// [[Rcpp::depends("RcppArmadillo")]]
//' The function find optimal unit weights
//' @param Y00 - T0 x N0 matrix of control unit outcomes before the treatment
//' @param Y10 - T0 x N - N0 matrix of treated unit outcomes before the treatment
//' @param dzeta - Ridge penalty parameter
//' @param rho - parameter used in ADMM
//' @param tol - Convergence tolerance for ADMM algorithm
//' @param weight_type - type of weight to be computed. Supported options are "SDID" and "Kernel"
//' @return N x 1 matrix of time weights for the estimation
//' @export
// [[Rcpp::export]]
arma::mat unit_weights(arma::mat &Y00, arma::mat &Y10, double dzeta, double rho, double tol, std::string weight_type)
{
  int T0 = Y00.n_rows;
  int N = Y00.n_cols + 1;
  arma::mat unit_weights(N - 1, 1);
  
  arma::mat eye_mat(N - 1, N - 1);
  eye_mat.eye();

  arma::mat V(N - 1, N - 1);

  V = Y00.t() * Y00/T0 + dzeta*eye_mat;
  //arma::colvec dvec = - Y10.t() * Y00/T0;
  arma::colvec dvec = - Y00.t() * Y10/T0;
  

  //depending on which weights the user has specified:
  if (weight_type == "SDID")
  {

    //prepare the input for ADMM
    //unit_weights = admm_solver(V, dvec, rho, 10000, 1e-8);
    
    
    arma::mat result(N, 1);
    result.submat(0, 0, N - 2, 0) = admm_solver(V, dvec, rho, 10000, 1e-8);
    result(N - 1, 0) = 1/(N - (N - 1));
    
    return result;
    //return unit_weights;
  }
  else if (weight_type == "Kernel")
  {
    exit(1);
  }
  else
  {
    Rcpp::Rcout << "Unknown weight type: " << weight_type << "\n";
    Rcpp::Rcout << "Please specify one of the supported weights \n";
  }
  
  exit(1);
}