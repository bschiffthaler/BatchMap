#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP CCOUNT(SEXP X, SEXP sequence) {
  BEGIN_RCPP
  mat _X = as<mat>(X);
  uvec _sequence = as<uvec>(sequence);

  mat _sM = _X.submat(_sequence.subvec(0, _sequence.size() - 2),
                      _sequence.subvec(1, _sequence.size() - 1));
  return wrap(accu(_sM.diag()));
  END_RCPP
}
