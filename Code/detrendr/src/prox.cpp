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
arma::vec prox_quantile(arma::vec w,
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


//' @title
//' Proximal mapping of f_1
//' 
//' \code{prox_f1} computes the proximal mapping of the average quantile loss
//'
//' @param theta input
//' @param y response
//' @param tau quantile parameter
//' @param step step-size
//' @examples
//' @export
// [[Rcpp::export]]
arma::vec prox_f1(arma::vec theta, 
                  arma::vec y, 
                  double tau = 0.05, 
                  double step = 1.0){
  int n = theta.n_elem;
  arma::vec w = y - theta;
  double alpha = step/n;
  return y - prox_quantile(w, tau, alpha);
}

//' @title
//' Proximal mapping of f_2
//' 
//' \code{prox_f2} computes the proximal mapping of the L1 penalty
//' 
//' @param eta input
//' @param lambda regularization parameter
//' @param step step-size
//' @export
//' @examples
//' set.seed(12345)
//' n <- 1e3
//' eta <- seq(-3, 3, length.out=n)
//' lambda <- 1
//' prox_out <- prox_f2(eta, lambda)
//' plot(eta, prox_out, type = 'l')
//' abline(0,1)
// [[Rcpp::export]]
arma::vec prox_f2(arma::vec eta, 
                  double lambda, 
                  double step = 1){
  return prox_quantile(eta, 0.5, 2*step*lambda);
}