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

//' @title
//' Proximal mapping
//' 
//' \code{prox} computes the block separable proximal mapping, changes theta
//' and eta in place
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export 
//[[Rcpp::export]]
void prox(arma::vec& theta, 
                arma::vec& eta, 
                arma::vec y, 
                double lambda, 
                double tau = 0.05, 
                double step = 1.0){
  theta = prox_f1(theta, y, tau, step);
  eta = prox_f2(eta, lambda, step);
}

//' @title
//' Proximal Mapping Test
//' 
//' \code{prox_test} computes the block separable proximal mapping. 
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export 
//[[Rcpp::export]]
Rcpp::List prox_test(arma::vec theta, 
                     arma::vec eta, 
                     arma::vec y, 
                     double lambda, 
                     double tau = 0.05, 
                     double step = 1.0){
  prox(theta, eta, y, lambda, tau, step);
  return Rcpp::List::create(theta = theta, eta = eta);
}

//' @title 
//' 
//' \code{get_D1} computes in the discrete derivative matrix.
//' 
//' @param n length of input
arma::sp_mat get_D1(int n){
  int numberNonZero = 2*(n-1);
  arma::vec values = ones<vec>(numberNonZero);
  values.subvec(n-1, 2*(n-1)-1) = -1*values.subvec(n-1, 2*(n-1)-1);
  
  arma::umat locs = repmat(linspace<urowvec>(0,n-2,n-1),2,2);
  locs.submat(1, n-1, 1, numberNonZero-1) = locs.submat(1, n-1, 1, 
              numberNonZero-1) + 1;
  
  arma::sp_mat D1 = arma::sp_mat(locs, values);

  //Rcout << "D1" << std::endl << D1 << std::endl;
  
  return D1;
}

//' @title
//' kth order difference matrix
//' 
//' \code{get_Dkn} computes the discrete kth derivative matrix
//' 
//' @param n length of input
//' @param k order of the derivative
//' 
//[[Rcpp::export]]
arma::sp_mat get_Dk(int n, 
                     int k){
  arma::sp_mat D = get_D1(n);
  for (int i=2; i < k+1; i++){
    D = get_D1(n-i+1)*D;
  }
  // Rcout << "D_k" << std::endl << D << std::endl;
  return D;
}


//' @title 
//' Project onto subspace (updates values of theta and eta in place)
//' 
//' \code{project_V} projects (theta, eta) onto the subspace eta = D%*%theta
//' 
//' @param theta first input
//' @param eta second input
//' @param D differencing matrix
//' @param cholM upper triangular cholesky of  I + DtD
//[[Rcpp::export]]
void project_V(arma::vec& theta,
                     arma::vec& eta,
                     arma::sp_mat D, 
                     arma::mat cholM){
  arma::vec DtEta = vectorise(D.t()*eta);
  theta = solve(trimatu(cholM), solve(trimatl(cholM.t()), theta + DtEta));
  eta = D*theta;
}


//' @title
//' One step of Spingarn's algorithm
//' 
//' \code{spingarn_one_step}
//' @param theta input 1
//' @param eta input 2
//' @param y response
//' @param D differencing matrix
//' @param cholM upper cholesky of  (I + DtD)
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size

//[[Rcpp::export]]
Rcpp::List spingarn_one_step(arma::vec theta, 
                             arma::vec eta, 
                             arma::vec y, 
                             arma::sp_mat D, 
                             arma::mat cholM,
                             double lambda, 
                             double tau = 0.05, 
                             double step = 1){
  arma::vec theta_old = theta;
  arma::vec eta_old = eta;
  prox(theta, eta, y, lambda, tau, step);
  arma::vec thetaMid = 2*theta-theta_old;
  arma::vec etaMid = 2*eta-eta_old;
  project_V(thetaMid, etaMid, D, cholM);

  theta = theta_old + 1.9*(thetaMid - theta);
  eta = eta_old + 1.9*(etaMid - eta);
  return Rcpp::List::create(theta=theta, eta=eta);
}

//' @title
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
//[[Rcpp::export]]
Rcpp::List spingarn_multiple(arma::vec theta,
                             arma::vec eta,
                             arma::vec y,
                             arma::sp_mat D, 
                             arma::mat cholM,
                             double lambda,
                             double tau = 0.05,
                             double step = 1,
                             double numberIter=1){

  for (int i = 0; i < numberIter; i++){
    spingarn_one_step(theta, eta, y, D, cholM, lambda, tau, step);
  }

  return Rcpp::List::create(theta=theta, eta=eta);
}
