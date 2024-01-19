#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List precision_recall(LogicalVector true_labels, NumericVector probs) {
  int n = probs.size();
  
  if (n != true_labels.size()) {
    stop("probs and true_labels must have the same length");
  }

  // Check the range of probabilities
  for (int i = 0; i < n; i++) {
    if (probs[i] < 0 || probs[i] > 1) {
      stop("All probability values must be between 0 and 1.");
    }
  }

  int num_thresholds = 101;  // For thresholds from 0 to 1 in increments of 0.01
  NumericVector thresholds(num_thresholds);
  NumericVector precision(num_thresholds);
  NumericVector recall(num_thresholds);
  NumericVector specificity(num_thresholds);  // Added for specificity

  for (int t = 0; t < num_thresholds; t++) {
    double threshold = t * 0.01;
    thresholds[t] = threshold;
    
    int TP = 0, FP = 0, FN = 0, TN = 0;

    for (int i = 0; i < n; i++) {
      bool pred = probs[i] >= threshold;
      if (pred && true_labels[i]) TP++;
      else if (pred && !true_labels[i]) FP++;
      else if (!pred && true_labels[i]) FN++;
      else TN++;
    }

    precision[t] = TP == 0 ? 0 : (double)TP / (TP + FP);
    recall[t] = TP == 0 ? 0 : (double)TP / (TP + FN);
    specificity[t] = TN + FP == 0 ? 0 : (double)TN / (TN + FP);  // Calculating specificity
  }

  return List::create(
    _["thresholds"] = thresholds,
    _["precision"] = precision,
    _["recall"] = recall,
    _["specificity"] = specificity  // Added for specificity
  );
}
