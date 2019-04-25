#include "att.h"
// [[Rcpp::depends("RcppArmadillo")]]
//' The function estimate the treatment effect
//' @param X - N*T x k data matrix
//' @param y - N*T vector of outcomes
//' @param unit_exclude - number of unit to exclude from the estimation. To include all units set this to -1
//' @param unit_weights - N x 1 matrix of unit weights
//' @param time_weights - T + 1 x 1 matrix of time weights
//' @param weight_type - type of weight to be computed. Supported options are "SDID" and "Kernel"
//' @export
//' @return att: Average Treatment Effect
// [[Rcpp::export]]
arma::mat att(arma::mat &X, arma::colvec &y, int unit_exclude, arma::mat &unit_weights, arma::mat &time_weights)
{
  //ordering here would be unit 1 over time, then unit 2 over time, etc
  
  //first construct the weighting matrix
  int N = unit_weights.n_rows;
  int T = X.n_rows/N;
  int k = X.n_cols;
  
  arma::vec W(N*T);

  for (int i = 0; i < N*T; i++)
  {
      W(i) = unit_weights(i/T)*time_weights(1 + i % T);
  }
  
  //the first element here would be the treatment effect
  //need to specify ordering in X
  arma::mat att(k, 1);
  
  arma::mat term_1(k, k);
  arma::mat term_2(k, 1);
  term_1.zeros();
  term_2.zeros();
  
  
  int i_unit;
  for (int i = 0; i < N*T; i++)
  {
    //if this is the unit we are excluding - continue
    i_unit = i/T;
    if (i_unit != unit_exclude)
    {
      term_1 = term_1 + W(i) * X.row(i).t() * X.row(i);
      term_2 = term_2 + W(i) * X.row(i).t() * y(i);
    }
  }
  
  if (unit_exclude >= 0)
  {
    term_1.shed_row(unit_exclude + 1);
    term_1.shed_col(unit_exclude + 1);
    term_2.shed_row(unit_exclude + 1);
  }

    
  att = solve(term_1, term_2);
  return att.submat(0, 0, 0, 0);
}