// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppEigen)]]
// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

//' @title
//' Perform cholesky decomposition and solve linear system
//' @param X sparse matrix
//' @param yy vector 
//' @export
// [[Rcpp::export]]
Eigen::VectorXd chol_solve_eigen(Eigen::MappedSparseMatrix<double>  X,
                                 Eigen::VectorXd yy){
  const Eigen::SparseMatrix<double> At((X).adjoint());
  const Eigen::VectorXd Aty(At * yy);
  const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > Ch(At * At.adjoint());
  return Ch.solve(Aty);
}

// [[Rcpp::export]]
SEXP chol_eigen(Eigen::MappedSparseMatrix<double> X){
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Lower, NaturalOrdering<int> >	cholOfX(X);	
  return wrap(cholOfX.matrixU());
}