// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GET_RF_MAT_NO_LOD
SEXP GET_RF_MAT_NO_LOD(SEXP seqnum, SEXP nmrk, SEXP CC, SEXP CR, SEXP RC, SEXP RR, SEXP minLOD, SEXP maxRF);
RcppExport SEXP _BatchMap_GET_RF_MAT_NO_LOD(SEXP seqnumSEXP, SEXP nmrkSEXP, SEXP CCSEXP, SEXP CRSEXP, SEXP RCSEXP, SEXP RRSEXP, SEXP minLODSEXP, SEXP maxRFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type seqnum(seqnumSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nmrk(nmrkSEXP);
    Rcpp::traits::input_parameter< SEXP >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< SEXP >::type CR(CRSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RC(RCSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RR(RRSEXP);
    Rcpp::traits::input_parameter< SEXP >::type minLOD(minLODSEXP);
    Rcpp::traits::input_parameter< SEXP >::type maxRF(maxRFSEXP);
    rcpp_result_gen = Rcpp::wrap(GET_RF_MAT_NO_LOD(seqnum, nmrk, CC, CR, RC, RR, minLOD, maxRF));
    return rcpp_result_gen;
END_RCPP
}
// READ_OUTCROSS
SEXP READ_OUTCROSS(SEXP file);
RcppExport SEXP _BatchMap_READ_OUTCROSS(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(READ_OUTCROSS(file));
    return rcpp_result_gen;
END_RCPP
}
// CCOUNT
SEXP CCOUNT(SEXP X, SEXP sequence);
RcppExport SEXP _BatchMap_CCOUNT(SEXP XSEXP, SEXP sequenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sequence(sequenceSEXP);
    rcpp_result_gen = Rcpp::wrap(CCOUNT(X, sequence));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP est_hmm_out(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_rf_out_wrap(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP get_bins(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BatchMap_GET_RF_MAT_NO_LOD", (DL_FUNC) &_BatchMap_GET_RF_MAT_NO_LOD, 8},
    {"_BatchMap_READ_OUTCROSS", (DL_FUNC) &_BatchMap_READ_OUTCROSS, 1},
    {"_BatchMap_CCOUNT", (DL_FUNC) &_BatchMap_CCOUNT, 2},
    {"est_hmm_out",     (DL_FUNC) &est_hmm_out,     6},
    {"est_rf_out_wrap", (DL_FUNC) &est_rf_out_wrap, 5},
    {"get_bins",        (DL_FUNC) &get_bins,        3},
    {NULL, NULL, 0}
};

RcppExport void R_init_BatchMap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
