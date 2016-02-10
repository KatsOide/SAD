#include <sim/sad_f2c.h>

#include <math.h>

/* atanh for REAL(8) from Fortran2008 standard */
real8 atanh_(real8* x) {
  return atanh(*x);
}

/* End of File */
