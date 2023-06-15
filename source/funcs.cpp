// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/LU>

using namespace Eigen;


// [[Rcpp::export]]
SEXP eigenMultSolveSp(SparseMatrix<double> X, SparseMatrix<double> XX, VectorXd Y){
  MatrixXd XtX = (XX.transpose() * X).pruned();
  VectorXd XtY = XX.transpose() * Y;
  
  VectorXd C = XtX.lu().solve(XtY);
  
  return Rcpp::wrap(C);
}
