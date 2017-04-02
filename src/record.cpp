#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP CCOUNT(SEXP X, SEXP sequence) {
  BEGIN_RCPP
  mat X_internal = as<mat>(X);
  uvec sequence_internal = as<uvec>(sequence);

  mat sM = X_internal.submat(sequence_internal.subvec(0, sequence_internal.size() - 2),
                      sequence_internal.subvec(1, sequence_internal.size() - 1));
  return wrap(accu(sM.diag()));
  END_RCPP
}
