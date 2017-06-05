#include <RcppArmadillo.h>
#define ARMA_USE_SUPERLU 1
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
//' \code{prox} computes the block separable proximal mapping. 
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export 
//[[Rcpp::export]]
Rcpp::List prox(arma::vec theta, 
                arma::vec eta, 
                arma::vec y, 
                double lambda, 
                double tau = 0.05, 
                double step = 1.0){
  return Rcpp::List::create(theta = prox_f1(theta, y, tau, step),
                            eta = prox_f2(eta, lambda, step));
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
//' Project onto subspace
//' 
//' \code{project_V} projects (theta, eta) onto the subspace eta = D%*%theta
//' 
//' @param theta first input
//' @param eta second input
//' @param D differencing matrix

// //[[Rcpp::export]]
// Rcpp::List project_V(arma::vec theta, 
//                      arma::vec eta, 
//                      arma::sp_mat D){
//   int n = D.n_cols;
//   arma::sp_mat M = speye<sp_mat>(n,n) + D.t()*D;
//   theta = spSolve(M, theta + D.t()*eta);
//   eta = D*theta;
//   return Rcpp::List::create(theta=theta, eta=eta);
// }

//[[Rcpp::export]]
Rcpp::List test_project_V(arma::vec theta, 
                          arma::vec eta, 
                          int n, 
                          int k){
  arma::sp_mat D = get_Dk(n, k);
  arma::sp_mat M = speye<sp_mat>(n,n) + D.t()*D;
  arma::vec DtEta = vectorise(D.t()*eta);
  //Rcout << "DtEta" << std::endl << DtEta << std::endl;
  theta = spsolve(M, theta + DtEta);
  arma::vec eta2 = vectorise(D*theta);
  return Rcpp::List::create(theta=theta, eta=eta);
}
