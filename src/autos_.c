#include <sim/sad_f2c.h>
#include <sim/TFCBK.h>

#include <stdio.h>
#include <math.h>

#include "gdtoaimp.h"

#define	COMPAT_SAD_AUTOS

void autos_(character result, ftnlen result_len, real8 *x) {
  char dtoa_buffer[result_len + 1];
  const char *buffer = dtoa_buffer;

#ifdef	COMPAT_SAD_AUTOS
  const char minus_infinity[] = "-INF";
  if(isinf(*x)) {
    buffer = (*x < 0) ? minus_infinity : (minus_infinity + 1);
  } else {
    g_dfmt(dtoa_buffer, x, 0, sizeof(dtoa_buffer));
  }
#else
  g_dfmt(dtoa_buffer, x, 0, sizeof(dtoa_buffer));
#endif

  /* Copy into Fortran character buffer */
  while(result_len > 0 && *buffer) {
    *result++ = *buffer++;
    --result_len;
  }
  while(result_len > 0) {
    *result++ = ' ';
    --result_len;
  }
}

void autoa_(character result, ftnlen result_len, real8 *x) {
  char dtoa_buffer[result_len + 1];
  const char *buffer = dtoa_buffer;

#ifdef	COMPAT_SAD_AUTOS
  const char minus_infinity[] = "-INF";
  if(isinf(*x)) {
    buffer = (*x < 0) ? minus_infinity : (minus_infinity + 1);
  } else if(isnan(*x)) {
    buffer = "NaN";
  } else {
    snprintf(dtoa_buffer, result_len + 1, "%a", *x);
  }
#else
  snprintf(dtoa_buffer, result_len + 1, "%a", *x);
#endif

  /* Copy into Fortran character buffer */
  while(result_len > 0 && *buffer) {
    *result++ = *buffer++;
    --result_len;
  }
  while(result_len > 0) {
    *result++ = ' ';
    --result_len;
  }
}

/* End of File */
