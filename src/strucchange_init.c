#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP recresid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"recresid", (DL_FUNC) &recresid, 13},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

