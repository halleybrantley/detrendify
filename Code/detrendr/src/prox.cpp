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
//' Function to test differencing matrix, returns single element of Dk(n)
//' 
//' @param n length of input
//' @param k order of the derivative
//' @param row row of element to return
//' @param col column of element to return
//' @export
//[[Rcpp::export]]
double test_Dk(int n, 
               int k, 
               int row,
               int col){
  arma::sp_mat Dk = get_Dk(n,k);
  return Dk(row, col);
}

//' @title 
//' Project onto subspace (updates values of theta and eta in place)
//' 
//' \code{project_V} projects (theta, eta) onto the subspace eta = D%*%theta
//' 
//' @param theta first input
//' @param eta second input
//' @param D differencing matrix
//' @param M = I + DtD
void project_V(arma::vec& theta,
                     arma::vec& eta,
                     arma::sp_mat D, 
                     arma::sp_mat M){
  arma::vec DtEta = vectorise(D.t()*eta);
  theta = spsolve(M, theta + DtEta, "lapack");
  eta = D*theta;
}

//' @title 
//' Test Cpp project_V function
//' 
//' \code{test_project_V}
//' @param theta input 1
//' @param eta input 2
//' @param k order of derivative
//[[Rcpp::export]]
Rcpp::List test_project_V(arma::vec theta, 
                          arma::vec eta, 
                          int k){
  int n = theta_n.elem;
  arma::sp_mat D = get_Dk(n, k);
  arma::sp_mat M = speye<sp_mat>(n,n) + D.t()*D;
  project_V(theta, eta, D, M);
  return Rcpp::List::create(theta=theta, eta=eta);
}

//' @title
//' One step of Spingarn's algorithm
//' 
//' \code{spingarn_one_step}
//' @param theta input 1
//' @param eta input 2
//' @param y response
//' @param D differencing matrix
//' @param M = I + DtD
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
void spingarn_one_step(arma::vec& theta, 
                             arma::vec& eta, 
                             arma::vec y, 
                             arma::sp_mat D, 
                             arma::sp_mat M,
                             double lambda, 
                             double tau = 0.05, 
                             double step = 1){
  arma::vec theta_old = theta;
  arma::vec eta_old = eta;
  prox(theta, eta, y, lambda, tau, step);
  arma::vec thetaMid = 2*theta-theta_old;
  arma::vec etaMid = 2*eta-eta_old;
  project_V(thetaMid, etaMid, D, M);
  theta = theta_old + 1.9*(thetaMid - theta);
  eta = eta_old + 1.9*(etaMid - eta);
}

//' @title
//' One step of Spingarn's algorithm
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
//' n <- 1e3
//' x <- seq(1/n, 1, length.out=n)
//' f <- 2*(x + 2)^2 + 3*cos(3*pi*x)
//' tau <- 1e4
//' g <-100*exp(-tau*(x-0.5)^2)
//' y <- f + g + rnorm(n)
//' k <- 3
//' D <- get_Dk_R(n, k)
//' M <- diag(n) + crossprod(D)
//' lambda <- 10
//' tau <- 0.01
//' step <- 1
//' numberIter <- 100
//' one_step <- spingarn_multi_iter(theta, eta, y, n, k, lambda, 
//' tau, step, numberIter)
//' theta <- one_step[[1]]
//' theta_last <- prox_f1(theta, y, tau)
//' plot(x,f,type='l',col='blue', ylim=c(min(y),max(y)), lwd=3)
//' points(x,y,pch=16)
//' lines(x,theta_last,col='red', lwd=3)
//[[Rcpp::export]]
Rcpp::List spingarn_multi_iter(arma::vec theta, 
                             arma::vec eta, 
                             arma::vec y, 
                             int k,
                             double lambda, 
                             double tau = 0.05, 
                             double step = 1, 
                             double numberIter=1){
  int n = y.n_elem;
  arma::sp_mat D = get_Dk(n, k);
  arma::sp_mat M = speye<sp_mat>(n,n) + D.t()*D;
  
  for (int i = 0; i < numberIter; i++){
    spingarn_one_step(theta, eta, y, D, M, lambda, tau, step);
  }
  
  return Rcpp::List::create(theta=theta, eta=eta);
}  
