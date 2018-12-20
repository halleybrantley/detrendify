#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include "getDk.hpp"
#include "prox.hpp"
using namespace Rcpp;
using namespace arma;

//'
//' One step of Spingarn's algorithm
//'
//' \code{spingarn_one_step} updates theta and eta in place
//' @param theta input 1
//' @param eta input 2
//' @param Vdiff vector 
//' @param y response
//' @param D differencing matrix
//' @param cholM upper cholesky of  (I + DtD)
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export
// [[Rcpp::export]]
void spingarn_one_step(arma::vec& theta,
                       arma::vec& eta,
                       arma::vec y,
                       arma::sp_mat D,
                       arma::sp_mat cholM,
                       double lambda,
                       double tau = 0.05,
                       double step = 1,
                       int k = 3){
  arma::vec theta_old = theta;
  arma::vec eta_old = eta;
  prox(theta, eta, y, lambda, tau, step);
  arma::vec thetaMid = 2*theta-theta_old;
  arma::vec etaMid = 2*eta-eta_old;
  project_V(thetaMid, etaMid, D, cholM, k);
  theta = theta_old + 1*(thetaMid - theta);
  eta = eta_old + 1*(etaMid - eta);
}

//'
//' Multiple steps of Spingarn's algorithm
//'
//' \code{spingarn_one_step}
//' @param theta input 1
//' @param eta input 2
//' @param y response
//' @param k order of derivative
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @param numberIter number of iterations
//' @export
//' @examples
//' set.seed(12345)
//' n <- 4e3
//' x <- seq(1/n, 1, length.out=n)
//' f <- 2*(x + 2)^2 + 3*cos(3*pi*x)
//' tau <- 5e4
//' g1 <- 40*exp(-tau*(x-0.3)^2)
//' g2 <- 30*exp(-3*tau*(x-0.5)^2)
//' g3 <- 37*exp(-tau*(x-0.7)^2)
//' g4 <- 45*exp(-tau/10*(x-0.55)^2)
//' y <- f + g1 + g2 + g3 + g4 + rnorm(n)
//' plot(y~x, type="l")
//' k <- 3
//' D <- get_Dk(n, k)
//' M <- diag(n) + crossprod(D)
//' cholM <- as.matrix(chol(M))
//' lambda <- 10
//' tau <- 0.01
//' step <- 1
//' numberIter <- 100
//' multi_step <- spingarn_multi_iter(theta, eta, y, n, k, lambda,
//' tau, step, numberIter)
//' theta <- multi_step[[1]]
//' theta_last <- prox_f1(theta, y, tau)
//' plot(x,f,type='l',col='blue', ylim=c(min(y),max(y)), lwd=3)
//' points(x,y,pch=16)
//' lines(x,theta_last,col='red', lwd=3)
// [[Rcpp::export]]
Rcpp::List spingarn_multi_step(arma::vec theta,
                               arma::vec eta,
                               arma::vec y,
                               arma::sp_mat D,
                               arma::sp_mat cholM,
                               double lambda,
                               double tau = 0.05,
                               double step = 1,
                               double numberIter=1,
                               int k=3){
  arma::vec Vdiff = zeros<vec>(1);
  arma::vec theta_cp = theta;
  arma::vec eta_cp = eta;
  
  for (int i = 1; i < numberIter; i++){
    spingarn_one_step(theta_cp, eta_cp, y, D, cholM,
                      lambda, tau, step, k);
    if (i % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
  }
  
  prox(theta_cp, eta_cp, y, lambda, tau, step);
  return Rcpp::List::create(Named("theta")=theta_cp,
                            _["eta"]=eta_cp);
}


