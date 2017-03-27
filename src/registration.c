#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BatchMap_CCOUNT(SEXP, SEXP);
extern SEXP BatchMap_GET_RF_MAT_NO_LOD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BatchMap_READ_OUTCROSS(SEXP);
extern SEXP est_hmm_out(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_rf_out_wrap(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_bins(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BatchMap_CCOUNT",            (DL_FUNC) &BatchMap_CCOUNT,            2},
    {"BatchMap_GET_RF_MAT_NO_LOD", (DL_FUNC) &BatchMap_GET_RF_MAT_NO_LOD, 8},
    {"BatchMap_READ_OUTCROSS",     (DL_FUNC) &BatchMap_READ_OUTCROSS,     1},
    {"est_hmm_out",                (DL_FUNC) &est_hmm_out,                6},
    {"est_rf_out_wrap",            (DL_FUNC) &est_rf_out_wrap,            5},
    {"get_bins",                   (DL_FUNC) &get_bins,                   3},
    {NULL, NULL, 0}
};

void R_init_BatchMap(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
