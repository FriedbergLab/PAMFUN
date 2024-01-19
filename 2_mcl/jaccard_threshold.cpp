#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::DataFrame Jaccard_cpp(const arma::mat& x, const arma::mat& y, 
                            Rcpp::CharacterVector rownames_x, Rcpp::CharacterVector rownames_y,
                            double threshold) {
  
  if (arma::any(arma::vectorise(x) > 1.0) || arma::any(arma::vectorise(y) > 1.0)) {
    Rcpp::warning("Matrices contain values greater than 1. Please binarize matrices before running Jaccard");
  }
  
  arma::mat intersection = x * y.t();
  
  arma::vec union_counts_x = arma::sum(x, 1);
  arma::vec union_counts_y = arma::sum(y, 1);
  
  arma::mat A = arma::repmat(union_counts_x, 1, intersection.n_cols);
  arma::mat B = arma::repmat(union_counts_y.t(), intersection.n_rows, 1);
  
  arma::mat jaccard_matrix = intersection / (A + B - intersection);

  // Filter the Jaccard similarity matrix based on the threshold and store the results
  std::vector<std::string> row_names;
  std::vector<std::string> col_names;
  std::vector<double> values;

  for(unsigned int i = 0; i < jaccard_matrix.n_rows; ++i) {
    for(unsigned int j = 0; j < jaccard_matrix.n_cols; ++j) {
      if(jaccard_matrix(i,j) >= threshold) {
        row_names.push_back(Rcpp::as<std::string>(rownames_x[i]));
        col_names.push_back(Rcpp::as<std::string>(rownames_y[j]));
        values.push_back(jaccard_matrix(i,j));
      }
    }
  }

  // Create a DataFrame to store the results
  return Rcpp::DataFrame::create(Rcpp::Named("protein1") = row_names,
                                 Rcpp::Named("protein2") = col_names,
                                 Rcpp::Named("sim") = values);
}
