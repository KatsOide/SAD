#include <random_driver.h>

#include <sim/sad_api.h>
#include <float.h>
#include <math.h>

/* Support NAN by using DBL_QNAN for non-C99 compiler */
#if !defined(NAN)       && defined(DBL_QNAN)
#define NAN DBL_QNAN
#endif

/* Fortran compatible function for src/tgauss.f */
real8 tgetgcut_(void) {
  return random_get_gcut();
}

void tsetgcut_(void) {
  random_set_gcut(rgetgl1("GCUT"));
}

real8 tran_(void) {
  double prn = NAN;

  /* Uniform distribution [0, 1) */
  random_generate(RANDOM_UNIFORM0, 1, &prn);
  return prn;
}

real8 tgauss_(void) {
  double prn = NAN;

  /* Gauss distribution with GCUT */
  random_generate(RANDOM_GAUSS, 1, &prn);
  return prn;
}

void tran_array_(real8 *prn, const integer4 *n) {
  /* Uniform distribution [0, 1) */
  random_generate(RANDOM_UNIFORM0, *n, prn);
}

void tgauss_array_(real8 *prn, const integer4 *n) {
  /* Gauss distribution with GCUT */
  random_generate(RANDOM_GAUSS, *n, prn);
}

/* End of File */
