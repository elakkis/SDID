#include "admm_solver.h"
// [[Rcpp::depends("RcppArmadillo")]]
//' Solves the  quadratic programming with linear equality 
//' and box inequality constraints (non-negative x that sum up to one) using ADMM algorithm.
//' @param Dmat : A (n x n) matrix representing  matrix involved in the quadratic term. Must be positive semidefinite.
//' @param dvec : A size n vector representing involved in the linear term of
//' @param rho : A scalar value specifying the augmented Lagrange multiplier rho for ADMM algorithm
//' @param max_iter : A integer indicating maximum number of iterations for ADMM algorithm.
//' @param tol : A scalar indicating the relative tolerance to determine the convergence of ADMM algorithm.
//' @return x : A size n vector containing optimal primal parameters from ADMM
//' @export
// [[Rcpp::export]]
arma::mat admm_solver(arma::mat &Dmat, arma::colvec &dvec,double rho, int max_iter,double tol)
{
  int n = Dmat.n_rows;
  int m = 1;

  //set up the constraints for weights
  arma::mat Amat(1, n);
  Amat.ones();

  arma::mat x(n, 1);
  arma::mat z(n, 1);
  arma::mat u(n, 1);
  
  x.zeros();
  z.zeros();
  u.zeros();
  

  arma::mat matr(m + n, m + n);
  matr.zeros();
  
  arma::mat tmp_eye(n, n, arma::fill::eye);
  matr.submat(0, 0, n - 1, n - 1) = Dmat + rho*tmp_eye;
  matr.submat(n, 0, n + m - 1, n - 1) = Amat;
  matr.submat(0, n, n - 1, n + m - 1) = Amat.t();
  

  matr.submat(n, n, n + m - 1, n + m - 1) = 0;
  

  //maybe we can get a speed gain here?
  //arma::mat m_inv = matr.i();
  arma::mat intercept(m + n, 1, arma::fill::zeros);
  intercept.submat(n, 0, n + m - 1, 0) = -1;
  

  double err =  tol*10;
  int i_iter = 0;
  for (i_iter = 0; i_iter < max_iter; i_iter++)
  {
    intercept.submat(0, 0, n - 1, 0) = dvec - rho*(z - u);

    //arma::mat tmp = -m_inv * intercept;
    arma::mat tmp = - solve(matr, intercept);
    x = tmp.submat(0, 0, n - 1, 0);
    z = x + u;
    
    //check that z belongs to the feasible region
    for (int i = 0; i < n; i++)
    {
      if (z(i, 0) < 0)
      {
        z(i, 0) = 0;
      }
      else if (z(i, 0) > 1)
      {
        z(i, 0) = 1;
      }
    }
    
    u = u + x - z;

    err = arma::norm(x - z, 1)/sqrt(pow(arma::norm(x, 1), 2) + pow(arma::norm(z, 1), 2));
    if (err <= tol)
    {
      break;
    }
  }
  
  if (i_iter == max_iter)
  {
    Rcpp::Rcout << "The algorithm has exceeded the maximum number of iterations \n";
  }
  
  return x;
}
