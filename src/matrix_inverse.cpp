// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List model_fit(const arma::mat& X, const arma::colvec& y) {

int n = X.n_rows, k = X.n_cols;

arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
arma::colvec fitted_value=X*coef;
arma::colvec res  = y - X*coef;           // residuals

// std.errors of coefficients
double MSE = std::inner_product(res.begin(), res.end(), res.begin(), 0.00)/(n-k);

arma::colvec std_err = arma::sqrt(MSE * arma::diagvec(arma::pinv(arma::trans(X)*X)));

arma::mat estimated_variance = MSE * arma::pinv(arma::trans(X)*X);

arma::mat x_inverse = arma::pinv(arma::trans(X)*X);

return List::create(Named("beta_estimates") = coef,
                    Named("fitted_value") = fitted_value,
                    Named("residual") = res,
                    Named("MSE") = MSE,
                    Named("beta_SE")  = std_err,
                    Named("estimated_variance") = estimated_variance,
                    Named("x_inverse")  = x_inverse);
}
