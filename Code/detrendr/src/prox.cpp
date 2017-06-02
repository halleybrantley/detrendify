#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
using namespace Rcpp;
using namespace arma;

//' @useDynLib detrendr
//' @importFrom Rcpp evalCpp
//' 
//' @title
//' \code{prox_quantile} computes the proximal mapping of the check function.
//'
//' @param w input
//' @param tau quantile parameter
//' @param alpha scale parameter
//' @examples
//' set.seed(12345)
//' n <- 1e3
//' w <- seq(-3, 3, length.out=n)
//' tau <- 0.5
//' alpha <- 2
//' prox_out <- prox_quantile(w, tau, alpha)
//' plot(w, prox_out, type='l', main=expression(paste(tau," = ")))
//'
//' tau <- 0.05
//' alpha <- 2
//' prox_out <- prox_quantile(w, tau, alpha)
//' plot(w, prox_out, type='l', main=expression(paste(tau," = ")))
//' @export
// [[Rcpp::export]]
arma::vec prox_quantileC(arma::vec w,
                        double tau,
                        double alpha){
  int n = w.n_elem;
  arma::vec prox_out = zeros<vec>(n);
  double threshold1 = tau*alpha;
  double threshold2 = -(1 - tau)*alpha;

  for (int i=0; i<n; i++){
    if (w(i) > threshold1) {
      prox_out(i) = w(i) - threshold1;
    } else if (w(i) < threshold2){
      prox_out(i) = w(i) - threshold2;
    }
  }
  return prox_out;
}
