#include <sim/sad_f2c.h>

#include <math.h>

/* hypot for REAL(8) from Fortran2008 standard */
real8 hypot_(real8* x, real8* y) {
  return hypot(*x, *y);
}

/* End of File */
