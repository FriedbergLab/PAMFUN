#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix Jaccard_cpp(const arma::mat& x, const arma::mat& y, Rcpp::CharacterVector rownames_x, Rcpp::CharacterVector rownames_y) {
  
  if (arma::any(arma::vectorise(x) > 1.0) || arma::any(arma::vectorise(y) > 1.0)) {
    Rcpp::warning("Matrices contain values greater than 1. Please binarize matrices before running Jaccard");
  }
  
  arma::mat intersection = x * y.t();
  
  arma::vec union_counts_x = arma::sum(x, 1);
  arma::vec union_counts_y = arma::sum(y, 1);
  
  arma::mat A = arma::repmat(union_counts_x, 1, intersection.n_cols);
  arma::mat B = arma::repmat(union_counts_y.t(), intersection.n_rows, 1);
  
  arma::mat jaccard_matrix = intersection / (A + B - intersection);
  
  // Create an Rcpp NumericMatrix to hold the result and set its row and column names
  Rcpp::NumericMatrix result = Rcpp::wrap(jaccard_matrix);
  rownames(result) = rownames_x;
  colnames(result) = rownames_y;
  
  return result;
}
