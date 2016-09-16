#ifndef _MACCODE_H_
#define _MACCODE_H_

/* sim for inc/MACCODE.inc */
#include <sim/sad_f2c.h>

enum {
  icNULL   =   0, icDRFT   =   1,
  icBEND   =   2, icQUAD   =   4, icSEXT   =   6,
  icOCTU   =   8, icDECA   =  10, icDODECA =  12,
  icUND    =  18, icWIG    =  19, icSOL    =  20,
  icST     =  21, icMULT   =  22, icTEST   =  30,
  icCAVI   =  31, icTCAV   =  32, icMAP    =  33,
  icINS    =  34, icCOORD  =  35, icBEAM   =  36,
  icPROT   =  37, icSPCH   =  38,
  icMARK   =  41, icMONI   =  42, icAPRT   =  43,
  icMXEL   =  99, icLINE   = 100, icCELL   = 100,

  icRSVD = 258,
  icDEF  = icRSVD + 2, icACT  = icDEF  + 2,
  icPROC = icACT  + 2, icVAR  = icPROC + 2,
  icKWRD = icVAR  + 2, icUNIT = icKWRD + 2,
  icRAND = icUNIT + 2, icENV  = icRAND + 2,
  icFLAG = icENV  + 2, icGLI  = icFLAG + 4,
  icGLL  = icGLI  + 4, icGLR  = icGLI  + 8,
  icGraf = icGLL  + 2, icPART = icGraf + 2
} SAD_MACCODE_T;

#endif /* _MACCODE_H_ */
