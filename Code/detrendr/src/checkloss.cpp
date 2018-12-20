#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
using namespace Rcpp;
using namespace arma;

//' @useDynLib detrendr
//' @importFrom Rcpp evalCpp

//' Check loss function
//' \code{check_loss} Evaluates check loss function
//' @param r vector of residuals
//' @param tau quantile level must be in [0,1]
//' @export
//'
// [[Rcpp::export]]
arma::mat checkloss(arma::mat r,
                  arma::vec tau){
  int nT = tau.n_elem;
  if (r.n_cols != nT) {
    throw(Rcpp::exception("Number of columns in y must be same as length of tau"));
  }

  arma::mat f = zeros<arma::vec>(nT);
  for (int j = 0; j < nT; j++){
    for (int i = 0; i < r.n_elem; i++){
      if (r(i,j) > 0){
        f(i,j) = tau(j)*r(i,j);
      } else {
        f(i,j) = (tau(j)-1)*r(i,j);
      }
    }
  }
  return f;
}
