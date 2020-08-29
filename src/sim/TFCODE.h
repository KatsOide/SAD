#ifndef _TFCODE_H_
#define _TFCODE_H_

/* sim for inc/TFCODE.inc */
#include <sim/sad_f2c.h>

typedef enum {
  ntfoper = 0, ntffun = 0,
  ntfreal = 1, ntflist = 3, ntflistr = 4, ntfstkseq = 5, ntfstk = 6,
  ntfstring = 101,
  ntfsymbol = 201, ntfdef = 201,
  ntfpat = 203, ntfarg = 205
} SAD_OBJECT_TYPES;

typedef enum {
  nfunlengyh = 20,
  nfundo     = 29,
  nfunmodule = 34,
  nfunblock  = 35,
  nfunif     = 39,
  nfununeval = 132,
  nfunwith   = 141
} SAD_FUNCTION_TYPES;

typedef enum {
  mtfnull  = 0, mtfneg   = 1, mtfplus  = 3, mtfminus = 4,
  mtfmult  = 5, mtfdiv   = 6, mtftimes = 5,
  mtfrevpower = 7, mtfpower    = 8, 
mtfgreater  = 9,
  mtfless     = 10, mtfgeq      = 11, mtfleq      = 12,
  mtfequal    = 13,  mtfunequal  = 14, 
  mtfand = 15, mtfor  = 16, mtfnot = 17, 
  mtfsame     = 18, mtfunsame   = 19,
  mtfconcat = 20,
  mtfleftbra    = 21, mtfrightbra   = 22,
  mtfleftbrace  = 23, mtfrightbrace = 24, mtflist = 23,
  mtfsetdelayed = 25, mtfset        = 26,
  mtfcomplex = 27,
  mtfleftparen  = 28, mtfrightparen = 29,
  mtfcomma = 30, mtfcomp  = 31, mtffun   = 32, mtfcolon = 33,
  mtfrule            = 34, mtfruledelayed     = 35,
  mtfreplace         = 36, mtfreplacerepeated = 37,
  mtfupset           = 38, mtfupsetdelayed    = 39,
  mtfunset = 40, mtfpattest = 41, mtfflag = 42,
  mtfslot = 43, mtfslotseq = 44, mtfdot = 45, mtfalt = 46,
  mtfmap = 47, mtfmapall = 48, mtfapply = 49, mtfrepeated = 50,
  mtfrepeatednull = 51, mtfinequality = 52, mtfaddto = 53,
  mtfsubtractfrom = 54, mtftimesby = 55, mtfdivideby = 56,
  mtfincrement = 57, mtfdecrement = 58,
  mtfpart = 59, mtfatt = 60, mtfmessagename = 61, mtftagset = 62,
  mtfleftcomment  = 63, mtfrightcomment = 64,
  mtfhold = 65, mtfend = 66,
  mtfnopc = 66 /* Last SAD_OPERATOR_TYPES */
} SAD_OPERATOR_TYPES;

typedef enum {
  lsimplepat   =  1, lsimplepatlist =   2, larglist       =   4,
  lconstlist   =  8, lnoconstlist   =  16, lnopatarg      =  32,
  lmemberlist  = 64, lnodefsymbol   = 128, lnoseqlist     = 256,
  lnonreallist = 512
} SAD_LIST_FLAGS;

typedef enum {
  knopatarg      =   lnonreallist*2,
  kpatarg   =    lnonreallist*4, kconstlist  =    lnonreallist*8,
  knoconstlist   =   lnonreallist*16, kfixedarg =   lnonreallist*32, 
  knofixedarg =   lnonreallist*64,
  kallnofixedarg =  lnonreallist*128, kconstarg =  lnonreallist*256, 
  knoconstarg =  lnonreallist*512,
  kseqarg        = lnonreallist*1024, knoseqarg = lnonreallist*2048
} SAD_ARG_FLAGS;

typedef enum {
  iattrholdfirst =  1, iattrholdrest  =  2,
  iattrconstant  =  4, iattrimmediate =  8,
  iattrorderless = 16, iattrdynamic   = 32,
  iattrholdall   =  3 /* holdfirst | holdrest */
} SAD_ATTRIBUTE_FLAGS;

static const integer8 ktfnull  =0xfff0000000000000;
static const integer8 ktfother =0xfff2000000000000;
static const integer8 ktfnr    =0x7ff2000000000000;
static const integer8 ktfoper  =0xfff6000000000000;
static const integer8 ktfref   =0xfffa000000000000;
static const integer8 ktfobj   =0x7ff2000000000000;
static const integer8 ktflist  =0x7ff2000000000000;
static const integer8 ktfpat   =0x7ff6000000000000;
static const integer8 ktfstring=0x7ffa000000000000;
static const integer8 ktfsymbol=0x7ffe000000000000;
static const integer8 ktomask  =0xfff2000000000000;
static const integer8 ktrmask  =0x7ff2000000000000;
static const integer8 ktfmask  =0xfffe000000000000;
static const integer8 ktamask  =0x0001ffffffffffff;
static const integer8 ktftrue  =0x3ff0000000000000;
static const integer8 ktffalse =0x0000000000000000;
static const integer8 ktfnan   =0xfff8000000000000;

static const real8 rtfnull = 0.0e0;

#endif /* _TFCODE_H_ */
