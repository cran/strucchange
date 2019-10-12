#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>


/* (1) Helper functions. */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt, names;
  PROTECT(elmt = R_NilValue);
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  for(int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }

  UNPROTECT(2);

  return elmt;
}


/* (2) Fallback evaluation. */
SEXP eval_fallback(SEXP fallback, SEXP r,
  SEXP fm, SEXP betar, SEXP check, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang5(fallback, r, fm, betar, check));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


/* (3) Recursive residuals. */
SEXP recresid(SEXP START, SEXP END, SEXP X1, SEXP XR,
  SEXP FR, SEXP BETAR, SEXP RVAL, SEXP X, SEXP Y,
  SEXP CHECK, SEXP FALLBACK, SEXP FM, SEXP RHO)
{
  int start = INTEGER(START)[0] - 1;
  int end = INTEGER(END)[0];
  int k = nrows(X1);
  int n = nrows(X);

  SEXP XT, XT2, XT3, FB, R, CHECK2, RVAL2;
  PROTECT(RVAL2 = duplicate(RVAL));
  PROTECT(CHECK2 = duplicate(CHECK));
  PROTECT(XT = duplicate(X1));
  PROTECT(XT2 = duplicate(X1));
  PROTECT(XT3 = duplicate(X1));
  PROTECT(R = duplicate(START));
  PROTECT(FB = eval_fallback(FALLBACK, R, FM, BETAR, CHECK2, RHO));

  double *X1ptr = REAL(X1);
  double *XRptr = REAL(XR);
  double *BETARptr = REAL(BETAR);
  double *RVALptr = REAL(RVAL2);
  double *Xptr = REAL(X);
  double *Yptr = REAL(Y);

  double *XTptr = REAL(XT);
  double *XT2ptr = REAL(XT2);
  double *XT3ptr = REAL(XT3);

  double *FB_X1ptr = 0;
  double *FB_betarptr = 0;

  double sum = 0.0;
  double sum2 = 0.0;
  double fr = REAL(FR)[0];
  double fr2 = pow(fr, 0.5);
  double fr0 = fr;
  double fr20 = fr2;

  int r, i, j, l;

  for(r = start; r < end; r++) {
    for(i = 0; i < k; i++) {
      for(j = 0; j < k; j++) {
        XTptr[i + j * k] = 0.0;
        XT2ptr[i + j * k] = 0.0;
        for(l = 0; l < k; l++) {
          XTptr[i + j * k] += X1ptr[i + l * k] * XRptr[l] * XRptr[j];
        }
      }
      for(j = 0; j < k; j++) {
        for(l = 0; l < k; l++) {
          XT2ptr[i + j * k] += XTptr[i + l * k] * X1ptr[l + j * k];
        }
        XT3ptr[i + j * k] = X1ptr[i + j * k] - (XT2ptr[i + j * k] / fr0);
        BETARptr[i] += XT3ptr[i + j * k] * XRptr[j] * RVALptr[r - start] * fr20;
      }
    }
    if(LOGICAL(CHECK2)[0]) {
      INTEGER(R)[0] = r + 1;
      FB = eval_fallback(FALLBACK, R, FM, BETAR, CHECK2, RHO);
      FM = getListElement(FB, "fm");
      LOGICAL(CHECK2)[0] = LOGICAL(getListElement(FB, "check"))[0];

      FB_X1ptr = REAL(getListElement(FB, "X1"));
      FB_betarptr = REAL(getListElement(FB, "betar"));

      for(i = 0; i < k; i++) {
        for(j = 0; j < k; j++) {
          XT3ptr[i + j * k] = FB_X1ptr[i + j * k];
        }
        BETARptr[i] = FB_betarptr[i];
      }
    }
    fr = 0.0;
    sum2 = 0.0;
    for(i = 0; i < k; i++) {
      sum = 0.0;
      for(j = 0; j < k; j ++) {
        sum += Xptr[r + j * n] * XT3ptr[j + i * k];
        X1ptr[j + i * k] = XT3ptr[j + i * k];
      }
      fr += sum * Xptr[r + i * n];
      sum2 += Xptr[r + i * n] * BETARptr[i];
      XRptr[i] = Xptr[r + i * n];
    }
    fr += 1.0;
    fr2 = pow(fr, 0.5);
    RVALptr[r - start + 1] = (Yptr[r] - sum2) / fr2;
  }

  UNPROTECT(7);

  return RVAL2;
}

