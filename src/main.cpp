#include <Utils.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  return Utils::dmvnrm_arma_fast(x, mean, sigma, logd);
}

