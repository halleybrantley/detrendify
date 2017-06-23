// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// chol_solve
arma::vec chol_solve(arma::sp_mat cholM, arma::vec b, int k, bool upper);
RcppExport SEXP detrendr_chol_solve(SEXP cholMSEXP, SEXP bSEXP, SEXP kSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type cholM(cholMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(chol_solve(cholM, b, k, upper));
    return rcpp_result_gen;
END_RCPP
}
// prox_quantile
arma::vec prox_quantile(arma::vec w, double tau, double alpha);
RcppExport SEXP detrendr_prox_quantile(SEXP wSEXP, SEXP tauSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_quantile(w, tau, alpha));
    return rcpp_result_gen;
END_RCPP
}
// prox_f1
arma::vec prox_f1(arma::vec theta, arma::vec y, double tau, double step);
RcppExport SEXP detrendr_prox_f1(SEXP thetaSEXP, SEXP ySEXP, SEXP tauSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_f1(theta, y, tau, step));
    return rcpp_result_gen;
END_RCPP
}
// prox_f2
arma::vec prox_f2(arma::vec eta, double lambda, double step);
RcppExport SEXP detrendr_prox_f2(SEXP etaSEXP, SEXP lambdaSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_f2(eta, lambda, step));
    return rcpp_result_gen;
END_RCPP
}
// prox
void prox(arma::vec& theta, arma::vec& eta, arma::vec y, double lambda, double tau, double step);
RcppExport SEXP detrendr_prox(SEXP thetaSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    prox(theta, eta, y, lambda, tau, step);
    return R_NilValue;
END_RCPP
}
// prox_test
Rcpp::List prox_test(arma::vec theta, arma::vec eta, arma::vec y, double lambda, double tau, double step);
RcppExport SEXP detrendr_prox_test(SEXP thetaSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_test(theta, eta, y, lambda, tau, step));
    return rcpp_result_gen;
END_RCPP
}
// get_D1
arma::sp_mat get_D1(int n);
RcppExport SEXP detrendr_get_D1(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_D1(n));
    return rcpp_result_gen;
END_RCPP
}
// get_Dk
arma::sp_mat get_Dk(int n, int k);
RcppExport SEXP detrendr_get_Dk(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(get_Dk(n, k));
    return rcpp_result_gen;
END_RCPP
}
// project_V
void project_V(arma::vec& theta, arma::vec& eta, arma::sp_mat D, arma::sp_mat cholM, int k);
RcppExport SEXP detrendr_project_V(SEXP thetaSEXP, SEXP etaSEXP, SEXP DSEXP, SEXP cholMSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type cholM(cholMSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    project_V(theta, eta, D, cholM, k);
    return R_NilValue;
END_RCPP
}
// spingarn_one_step
void spingarn_one_step(arma::vec& theta, arma::vec& eta, arma::vec y, arma::sp_mat D, arma::sp_mat cholM, double lambda, double tau, double step, int k);
RcppExport SEXP detrendr_spingarn_one_step(SEXP thetaSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP DSEXP, SEXP cholMSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP stepSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type cholM(cholMSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    spingarn_one_step(theta, eta, y, D, cholM, lambda, tau, step, k);
    return R_NilValue;
END_RCPP
}
// spingarn_multi_step
Rcpp::List spingarn_multi_step(arma::vec theta, arma::vec eta, arma::vec y, arma::sp_mat D, arma::sp_mat cholM, double lambda, double tau, double step, double numberIter, int k);
RcppExport SEXP detrendr_spingarn_multi_step(SEXP thetaSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP DSEXP, SEXP cholMSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP stepSEXP, SEXP numberIterSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type cholM(cholMSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< double >::type numberIter(numberIterSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(spingarn_multi_step(theta, eta, y, D, cholM, lambda, tau, step, numberIter, k));
    return rcpp_result_gen;
END_RCPP
}
// chol_solve_eigen
Eigen::VectorXd chol_solve_eigen(Eigen::MappedSparseMatrix<double> X, Eigen::VectorXd yy);
RcppExport SEXP detrendr_chol_solve_eigen(SEXP XSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(chol_solve_eigen(X, yy));
    return rcpp_result_gen;
END_RCPP
}
// chol_eigen
SEXP chol_eigen(Eigen::MappedSparseMatrix<double> X);
RcppExport SEXP detrendr_chol_eigen(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(chol_eigen(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"detrendr_chol_solve", (DL_FUNC) &detrendr_chol_solve, 4},
    {"detrendr_prox_quantile", (DL_FUNC) &detrendr_prox_quantile, 3},
    {"detrendr_prox_f1", (DL_FUNC) &detrendr_prox_f1, 4},
    {"detrendr_prox_f2", (DL_FUNC) &detrendr_prox_f2, 3},
    {"detrendr_prox", (DL_FUNC) &detrendr_prox, 6},
    {"detrendr_prox_test", (DL_FUNC) &detrendr_prox_test, 6},
    {"detrendr_get_D1", (DL_FUNC) &detrendr_get_D1, 1},
    {"detrendr_get_Dk", (DL_FUNC) &detrendr_get_Dk, 2},
    {"detrendr_project_V", (DL_FUNC) &detrendr_project_V, 5},
    {"detrendr_spingarn_one_step", (DL_FUNC) &detrendr_spingarn_one_step, 9},
    {"detrendr_spingarn_multi_step", (DL_FUNC) &detrendr_spingarn_multi_step, 10},
    {"detrendr_chol_solve_eigen", (DL_FUNC) &detrendr_chol_solve_eigen, 2},
    {"detrendr_chol_eigen", (DL_FUNC) &detrendr_chol_eigen, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_detrendr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
