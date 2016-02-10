#include <sim/sad_f2c.h>

#include <stdlib.h>
#include <stdio.h>

static FILE *npin = NULL, *npout = NULL;

void nfputm_(integer4 *lfno,
	     real8 *xa, real8 *ya, real8 *xxa, real8 *xya, real8 *yya) {
  char *pout;

  if(npout == NULL) {
    pout = getenv("SADPOUT");
    fprintf(stderr, " opening pipe:%s for write\n", pout);
    npout = fopen(pout, "a");
    fprintf(stderr, " %s is opened for write\n", pout);
    if(npout == NULL) {
      fprintf(stderr, " putm:unable to open fifo %s\n", pout);
      exit(1);
    }
  }

  fprintf(npout,
	  "%15.7lg,%15.7lg,%15.7lg,%15.7lg,%15.7lg\n",
	  *xa, *ya, *xxa, *xya, *yya);
  fflush(npout);

  fprintf(stderr,
	  "Putting data:%15.7lg,%15.7lg,%15.7lg,%15.7lg,%15.7lg \n",
	  *xa, *ya, *xxa, *xya, *yya); 
  fflush(stderr);
}

void nfgetm_(integer4 *lfno,
	     real8 *xa, real8 *ya, real8 *xxa, real8 *xya, real8 *yya) {
  char *pin;

  if(npin == NULL) {
    pin = getenv("SADPIN");
    npin = fopen(pin, "r");
    if(npin == NULL) {
      fprintf(stderr, " getm:unable to open fifo %s\n", pin);
      exit(1);
    }
  }

  fscanf(npin, "%lg,%lg,%lg,%lg,%lg",  xa, ya, xxa, xya, yya);

  fprintf(stderr,
	  "data gotten from :%15.7lg,%15.7lg,%15.7lg,%15.7lg,%15.7lg \n",
	  *xa, *ya, *xxa, *xya, *yya);
  fflush(stderr); 
}

/* End of File */
