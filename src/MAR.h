#include <RcppArmadillo.h> // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// [[Rcpp::depends(RcppArmadillo)]]

class MAR{ // Multivariate autoregressive process class
  // Model Y_t = m + A * Y_{t - 1} + E,  E is a N(0, Sigma) vector
private:
  // Attributes
  unsigned int dimension; // Dimension of the processed
  arma::rowvec m; // intercept dimension * 1
  arma::mat A; // Dynamics matrix dimension * dimension
  arma::mat Sigma; // Innovation var-covariance matrix dimension * dimension
  arma::rowvec mu0;
  arma::mat Sigma0;
public:
  // Constructor
  // Getters
  double get_log_likelihood(arma::mat observations){
    unsigned int n = observations.n_rows - 1;
    arma::mat Y = observations.rows(1, n);
    arma::mat X = arma::join_cols()
  }
}
