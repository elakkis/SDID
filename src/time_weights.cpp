#include "time_weights.h"
#include "admm_time_solver.h"
// [[Rcpp::depends("RcppArmadillo")]]
//' The function find optimal time weights
//' @param Y00 - T0 x N0 matrix of control unit outcomes before the treatment
//' @param Y01 - T-T0 x N0 matrix of control unit outcomes after the treatment
//' @param dzeta - Ridge penalty parameter
//' @param rho - parameter used in ADMM
//' @param tol - Convergence tolerance for ADMM algorithm
//' @param weight_type - type of weight to be computed. Supported options are "SDID" and "Kernel"
//' @return T + 1 x 1 matrix of time weights for the estimation
//' @export
// [[Rcpp::export]]
arma::mat time_weights(arma::mat &Y00, arma::mat &Y01, double dzeta, double rho, double tol, std::string weight_type)
{
  int T0 = Y00.n_rows;
  int N = Y00.n_cols + 1;
  
  int T = T0 + Y01.n_rows;
  arma::mat result(T, 1);
  
  arma::mat V(T0 + 1, T0 + 1);
  V.zeros();
  
  arma::rowvec tmp_vec(N - 1);
  tmp_vec.ones();
  arma::mat eye_mat(T0, T0);
  eye_mat.eye();
  
  V(0, 0) = 1;
  V.submat(1, 1, T0, T0) = Y00 * Y00.t()/(N - 1) + dzeta*eye_mat;
  V.submat(1, 0, T0, 0) = Y00 * tmp_vec.t()/(N - 1);
  V.submat(0, 1, 0, T0) = tmp_vec * Y00.t()/(N - 1);
  //V = V + dzeta*eye_mat;

  //need to compute the average of Y01 over time (average for each column)
  arma::colvec Y01_av(N - 1);
  
  tmp_vec = mean(Y01, 0);
  Y01_av = tmp_vec.t();

  arma::colvec dvec(T0 + 1);
  dvec(0) = -sum(Y01_av)/(N - 1);

  dvec.subvec(1, T0) = - Y00 * Y01_av/(N - 1);

  //depending on which weights the user has specified:
  if (weight_type == "SDID")
  {
    //prepare the input for ADMM
    //ignore lambda 0 here
    result.submat(0, 0, T0 - 1, 0) = admm_time_solver(V, dvec, rho, 10000, tol).submat(1, 0, T0, 0);
    
    for (int i = T0; i < T; i++)
    {
      result(i, 0) = 1.0/(T - T0);
    }
    return result;
  }
  else if (weight_type == "Kernel")
  {
    //add kernel weights here and option for kernel functions
    exit(1);
  }
  else
  {
    Rcpp::Rcout << "Unknown weight type: " << weight_type << "\n";
    Rcpp::Rcout << "Please specify one of the supported weights \n";
  }
  
  exit(1);
}