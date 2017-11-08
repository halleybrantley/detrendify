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
double check_loss(arma::vec r, 
                  double tau){
  double f = 0.0;
  for (int i = 0; i < r.n_elem; i++){
    if (r(i) > 0){
      f += tau*r(i);
    } else {
      f += (tau-1)*r(i);
    }
  }
  return f;
}


//' Banded Cholesky Solve
//' 
//' \code{chol_solve} Solves a linear system cholM%*%x=b 
//' when cholM is a sparse banded cholesky.
//' 
//' @param cholM sparse banded cholseky decomposition of discrete derivative 
//' matrix of order k
//' @param b dense solution vector
//' @param k order of discrete derivative matrix
//' @param upper boolean indicator of whether cholM is upper or lower triangular
//' @export
//' 
// [[Rcpp::export]]
arma::vec chol_solve(arma::sp_mat cholM, 
                     arma::vec b, 
                     int k,
                     bool upper = true){
  int n = cholM.n_cols;
  arma::uvec ind;
  arma::vec x=b;
  
  if(upper){
    ind = regspace<uvec>(n-1, 0);
  } else {
    ind = regspace<uvec>(0, n-1);
  }
  
  for (int i=0; i < k+1; i++){
    for (int j=0; j < i; j++){
      x(ind(i)) -= cholM(ind(i), ind(j))*x(ind(j));   
    }
    x(ind(i)) /= cholM(ind(i), ind(i));
  }
  
  for (int i=k+1; i < n; i++){
    for (int j = i-k; j < i; j++){
      x(ind(i)) -= cholM(ind(i), ind(j))*x(ind(j)); 
    }
    x(ind(i)) /= cholM(ind(i), ind(i));
  }
  return x;
}


//' Proximal Mapping
//' 
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



//' Proximal mapping of f_1
//' 
//' \code{prox_f1} computes the proximal mapping of the average quantile loss
//'
//' @param theta input
//' @param y response
//' @param tau quantile parameter
//' @param step step-size
//' @export
// [[Rcpp::export]]
arma::vec prox_f1(arma::vec theta, 
                  arma::vec y, 
                  double tau = 0.05, 
                  double step = 1.0){
  int n = theta.n_elem;
  arma::vec w = y - theta;
  return y - prox_quantile(w, tau, step);
}


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


//' Proximal mapping
//' 
//' \code{prox} computes the block separable proximal mapping, changes theta
//' and eta in place
//' 
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export 
// [[Rcpp::export]]
void prox(arma::vec& theta, 
                arma::vec& eta, 
                arma::vec y, 
                double lambda, 
                double tau = 0.05, 
                double step = 1.0){
  theta = prox_f1(theta, y, tau, step);
  eta = prox_f2(eta, lambda, step);
}

//' Proximal Mapping Test
//' 
//' \code{prox_test} computes the block separable proximal mapping.
//' Returns values of theta and eta. 
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export 
// [[Rcpp::export]]
Rcpp::List prox_test(arma::vec theta, 
                     arma::vec eta, 
                     arma::vec y, 
                     double lambda, 
                     double tau = 0.05, 
                     double step = 1.0){
  prox(theta, eta, y, lambda, tau, step);
  return Rcpp::List::create(theta = theta, eta = eta);
}

//' Discrete derivative matrix
//' 
//' \code{get_D1} computes the sparse discrete derivative matrix.
//' 
//' @param n length of input
//' @examples
//' n <- 5
//' D1 <- get_D1(n)
//' @export
// [[Rcpp::export]]
arma::sp_mat get_D1(int n){
  int numberNonZero = 2*(n-1);
  arma::vec values = ones<vec>(numberNonZero);
  values.subvec(n-1, 2*(n-1)-1) = -1*values.subvec(n-1, 2*(n-1)-1);
  arma::umat locs = repmat(linspace<urowvec>(0,n-2,n-1),2,2);
  locs.submat(1, n-1, 1, numberNonZero-1) = locs.submat(1, n-1, 1, 
              numberNonZero-1) + 1;
  arma::sp_mat D1 = arma::sp_mat(locs, values);
  return D1;
}

//' kth order sparse difference matrix
//' 
//' \code{get_Dkn} computes the sparse discrete kth derivative matrix
//' 
//' @param n length of input
//' @param k order of the derivative
//' @export
// [[Rcpp::export]]
arma::sp_mat get_Dk(int n, 
                     int k){
  arma::sp_mat D = get_D1(n);
  for (int i=2; i < k+1; i++){
    D = get_D1(n-i+1)*D;
  }
  // Rcout << "D_k" << std::endl << D << std::endl;
  return D;
}

//' Project onto subspace 
//' 
//' \code{project_V} projects (theta, eta) onto the subspace eta = D*theta.
//' Updates values of theta and eta in place.
//' 
//' @param theta first input
//' @param eta second input
//' @param D differencing matrix
//' @param cholM upper triangular cholesky decomposition of  I + DtD
//' @export
// [[Rcpp::export]]
void project_V(arma::vec& theta,
                     arma::vec& eta,
                     arma::sp_mat D, 
                     arma::sp_mat cholM, 
                     int k){
  arma::vec DtEta = vectorise(D.t()*eta);
  theta = chol_solve(cholM, chol_solve(cholM.t(), theta + DtEta, k, false), 
                     k, true);
  eta = D*theta;
}


//' 
//' One step of Spingarn's algorithm
//' 
//' \code{spingarn_one_step} updates theta and eta in place
//' @param theta input 1
//' @param eta input 2
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
  arma::vec theta_cp = theta;
  arma::vec theta_cur = theta_cp;
  arma::vec eta_cp = eta;
  arma::vec f1 = zeros<vec>(numberIter);
  arma::vec f2 = zeros<vec>(numberIter);
  arma::vec fmean = zeros<vec>(numberIter);
  f1(0) = check_loss(y-prox_f1(theta_cp, y, tau, step), tau);
  f2(0) = lambda*norm(prox_f2(eta_cp, lambda, step),1);
  fmean(0) =  f1(0)+f2(0);
  int endIter = numberIter;
  int avgNum = 200;
  int n = theta.n_elem;
  if (lambda/n > 200){
    avgNum = lambda/n;
  }
  for (int i = 1; i < numberIter; i++){
    spingarn_one_step(theta_cp, eta_cp, y, D, cholM, 
                      lambda, tau, step, k);
    f1(i) = norm(theta_cp);
    // f1(i) = check_loss(y-prox_f1(theta_cp, y, tau, step), tau);
    // f2(i) = lambda*norm(prox_f2(eta_cp, lambda, step),1);
    // fmean(i) = (fmean(i-1)*(avgNum-1) + f1(i) + f2(i))/(avgNum);
    // if(i > avgNum){
    //   if (fmean(i)-fmean(i-1) > 0) {
    //     endIter = i;
    //     break;
    //   }
    // }
    if (i % 100 == 0){
      Rcpp::checkUserInterrupt(); 
    }
  }
  
  return Rcpp::List::create(theta_cp, eta_cp, fmean, f1, f2, endIter);
}
