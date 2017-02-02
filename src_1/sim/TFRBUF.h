#ifndef _TFRBUF_H_
#define _TFRBUF_H_

/* sim for inc/TFRBUF.inc */
#include <sim/sad_f2c.h>

enum {
  irbinit          =  1, irbopen     =  2, irbclose     =  3,
  irbreadrecord    =  4, irbreadbuf  =  5, irbmovepoint =  6,
  irbbor           =  7, irbgetpoint =  8, irbreset     =  9,
  irbreadrecordbuf = 10, irbeor2bor  = 11, irbsetinp    = 12,
  irbcloseinp      = 13, irbsetbuf   = 14, irbibuf      = 15
} SAD_RBUF_COMMANDS;

#endif /* _TFRBUF_H_ */
