#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
SEXP GET_RF_MAT_NO_LOD(SEXP seqnum, SEXP nmrk, SEXP CC,
                       SEXP CR, SEXP RC, SEXP RR, SEXP minLOD,
                       SEXP maxRF)
{
  uvec _seqnum = as<uvec>(seqnum);
  mat _cc = as<mat>(CC);
  mat _cr = as<mat>(CR);
  mat _rc = as<mat>(RC);
  mat _rr = as<mat>(RR);
  uword _nmrk = as<uword>(nmrk);
  double _minLOD = as<double>(minLOD);
  double _maxRF = as<double>(maxRF);

  mat r(_nmrk, _nmrk);
  std::fill( r.begin(), r.end(), NumericVector::get_na() ) ;

  for(uword i = 0; i < _nmrk - 1; i++)
  {
    for(uword j = i+1; j < _nmrk; j++)
    {
      uword k1 = _seqnum(i) - 1;
      uword k2 = _seqnum(j) - 1;
      if(k1 > k2){
        uword tmp = k1;
        k1 = k2; k2 = tmp;
      }
      vec rfs(4); vec lods(4);
      rfs(0) = _cc(k2,k1); lods(0) = _cc(k1,k2);
      rfs(1) = _cr(k2,k1); lods(1) = _cr(k1,k2);
      rfs(2) = _rc(k2,k1); lods(2) = _rc(k1,k2);
      rfs(3) = _rr(k2,k1); lods(3) = _rr(k1,k2);
      uword good = 0;
      double rfatlomax = 0;
      double lomax = 0;
      for(uword x = 0; x < 4; x++)
      {
        if(rfs(x) <= _maxRF && lods(x) >= _minLOD)
        {
          good++;
          if(lods(x) > lomax){
            lomax = lods(x);
            rfatlomax = rfs(x);
          }
        }
      }
      if(good == 0) continue;
      if(rfatlomax > 0.5) rfatlomax = 0.5;
      r(i,j) = rfatlomax;
      r(j,i) = rfatlomax;
    }
  }
  return wrap(r);
}

