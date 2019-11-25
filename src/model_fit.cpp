// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// This is a Rcpp function used for least-square matrix operation.

// [[Rcpp::export]]
List model_fit(const arma::mat& X, const arma::colvec& y) {

  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);    // estimate coefficients for model y ~ X

  arma::colvec fitted_value=X*coef;         // estimated y_hat

  arma::colvec res  = y - X*coef;           // residuals


  double MSE = std::inner_product(res.begin(), res.end(), res.begin(), 0.00)/(n-k);  // mean squared error

  arma::colvec std_err = arma::sqrt(MSE * arma::diagvec(arma::pinv(arma::trans(X)*X))); // std.errors of coefficients

  arma::mat estimated_variance = MSE * arma::pinv(arma::trans(X)*X); // estimated variance and covariance for coefficients

  arma::mat x_inverse = arma::pinv(arma::trans(X)*X); // tanspose(x)%*%x

  return List::create(Named("beta_estimates") = coef,
                      Named("fitted_value") = fitted_value,
                      Named("residual") = res,
                      Named("MSE") = MSE,
                      Named("beta_SE")  = std_err,
                      Named("estimated_variance") = estimated_variance,
                      Named("x_inverse")  = x_inverse);
}
