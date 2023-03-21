#ifndef _TFCBK_H_
#define _TFCBK_H_

/* sim for inc/TFCBK.inc */
#include <sim/sad_f2c.h>

static const integer4 maxgeneration = 1L << 30;
static const integer4 maxlevele     = 4096;
static const integer4 nsymhash      = 2047;
static const integer4 maxlbuf       = 32768;
static const integer4 nslots        = 32;

/* tffvp common block */
/*
typedef struct {
  integer8
  __itfcontroot, __itfcontext, __itfcontextpath, __itflocal,
    __kxeof, __kxfailed, __iaxhold, __iaximmediate,
    __iaxline, __kxliteral, __kxnull, __kxnulll, __kxnulls,
    __iaxout, __iaxpriority, __iaxschar, 
    __iaxslotnull, __iaxslotpart, __iaxslotseqnull,
    __kxvect, __kxvect1, __iavpw,
    __kerror, __ierrorf, __ierrorgen, __ierrorprint, __ierrorth,
    __ierrorexp, __ifunbase, __initmessage, __levelcompile;
  real8 __dinfinity, __dnotanumber;
  integer4  
  __levele, __levelp, __lgeneration, __ltrace,
    __modethrow, __iordless;
} tffvp_t;
extern tffvp_t tffvp_;
*/
extern integer8 __tfcbk_MOD_itfcontroot;
extern integer8 __tfcbk_MOD_itfcontext;
extern integer8 __tfcbk_MOD_itfcontextpath;
extern integer8 __tfcbk_MOD_itflocal;

extern integer8 __tfcbk_MOD_kxeof;
extern integer8 __tfcbk_MOD_kxfailed;
extern integer8 __tfcbk_MOD_iaxhold;
extern integer8 __tfcbk_MOD_iaximmediate;

extern integer8 __tfcbk_MOD_iaxline;
extern integer8 __tfcbk_MOD_kxliteral;
extern integer8 __tfcbk_MOD_kxnull;
extern integer8 __tfcbk_MOD_kxnulll;
extern integer8 __tfcbk_MOD_kxnulls;

extern integer8 __tfcbk_MOD_iaxout;
extern integer8 __tfcbk_MOD_iaxpriority;
extern integer8 __tfcbk_MOD_iaxschar;

extern integer8 __tfcbk_MOD_iaxslotnull;
extern integer8 __tfcbk_MOD_iaxslotpart;
extern integer8 __tfcbk_MOD_iaxslotseqnull;

extern integer8 __tfcbk_MOD_kxvect;
extern integer8 __tfcbk_MOD_kxvect1;
extern integer8 __tfcbk_MOD_iavpw;

extern integer8 __tfcbk_MOD_kerror;
extern integer8 __tfcbk_MOD_ierrorf;
extern integer8 __tfcbk_MOD_ierrorgen;
extern integer8 __tfcbk_MOD_ierrorprint;
extern integer8 __tfcbk_MOD_ierrorth;

extern integer8 __tfcbk_MOD_ierrorexp;
extern integer8 __tfcbk_MOD_ifunbase;
extern integer8 __tfcbk_MOD_initmessage;
extern integer8 __tfcbk_MOD_levelcompile;

extern real8 __tfcbk_MOD_dinfinity;
extern real8 __tfcbk_MOD_dnotanumber;

extern integer4 __tfcbk_MOD_levele;
extern integer4 __tfcbk_MOD_levelp;
extern integer4 __tfcbk_MOD_lgeneration;
extern integer4 __tfcbk_MOD_ltrace;
extern integer4 __tfcbk_MOD_modethrow;
extern integer4 __tfcbk_MOD_iordless;

/* tffvp common block macro */
/*
#define itfcontroot	(tffvp_.__itfcontroot)
#define itfcontext	(tffvp_.__itfcontext)
#define itfcontextpath	(tffvp_.__itfcontextpath)
#define itflocal	(tffvp_.__itflocal)

#define kxeof	(tffvp_.__kxeof)
#define kxfailed	(tffvp_.__kxfailed)
#define iaxhold	(tffvp_.__iaxhold)
#define iaximmediate	(tffvp_.__iaximmediate)

#define iaxline	(tffvp_.__iaxline)
#define kxliteral	(tffvp_.__kxliteral)
#define kxnull	(tffvp_.__kxnull)
#define kxnulll	(tffvp_.__kxnulll)
#define kxnulls	(tffvp_.__kxnulls)

#define iaxout	(tffvp_.__iaxout)
#define iaxpriority	(tffvp_.__iaxpriority)
#define iaxschar	(tffvp_.__iaxschar)

#define iaxslotnull	(tffvp_.__iaxslotnull)
#define iaxslotpart	(tffvp_.__iaxslotpart)
#define iaxslotseqnull	(tffvp_.__iaxslotseqnull)

#define kxvect	(tffvp_.__kxvect)
#define kxvect1	(tffvp_.__kxvect1)
#define iavpw	(tffvp_.__iavpw)

#define kerror	(tffvp_.__kerror)
#define ierrorf	(tffvp_.__ierrorf)
#define ierrorgen	(tffvp_.__ierrorgen)
#define ierrorprint	(tffvp_.__ierrorprint)
#define ierrorth	(tffvp_.__ierrorth)

#define ierrorexp	(tffvp_.__ierrorexp)
#define ifunbase	(tffvp_.__ifunbase)
#define initmessage	(tffvp_.__initmessage)

#define dinfinity	(tffvp_.__dinfinity)
#define dnotanumber	(tffvp_.__dnotanumber)

#define levelcompile	(tffvp_.__levelcompile)
#define levele	(tffvp_.__levele)
#define levelp	(tffvp_.__levelp)
#define lgeneration	(tffvp_.__lgeneration)
#define ltrace	(tffvp_.__ltrace)
#define modethrow	(tffvp_.__modethrow)
#define iordless	(tffvp_.__iordless)
*/

#define itfcontroot	__tfcbk_MOD_itfcontroot
#define itfcontext	__tfcbk_MOD_itfcontext
#define itfcontextpath	__tfcbk_MOD_itfcontextpath
#define itflocal	__tfcbk_MOD_itflocal

#define kxeof	__tfcbk_MOD_kxeof
#define kxfailed	__tfcbk_MOD_kxfailed
#define iaxhold	__tfcbk_MOD_iaxhold
#define iaximmediate	__tfcbk_MOD_iaximmediate

#define iaxline	__tfcbk_MOD_iaxline
#define kxliteral	__tfcbk_MOD_kxliteral
#define kxnull	__tfcbk_MOD_kxnull
#define kxnulll	__tfcbk_MOD_kxnulll
#define kxnulls	__tfcbk_MOD_kxnulls

#define iaxout	__tfcbk_MOD_iaxout
#define iaxpriority	__tfcbk_MOD_iaxpriority
#define iaxschar	__tfcbk_MOD_iaxschar

#define iaxslotnull	__tfcbk_MOD_iaxslotnull
#define iaxslotpart	__tfcbk_MOD_iaxslotpart
#define iaxslotseqnull	__tfcbk_MOD_iaxslotseqnull

#define kxvect	__tfcbk_MOD_kxvect
#define kxvect1	__tfcbk_MOD_kxvect1
#define iavpw	__tfcbk_MOD_iavpw

#define kerror	__tfcbk_MOD_kerror
#define ierrorf	__tfcbk_MOD_ierrorf
#define ierrorgen	__tfcbk_MOD_ierrorgen
#define ierrorprint	__tfcbk_MOD_ierrorprint
#define ierrorth	__tfcbk_MOD_ierrorth

#define ierrorexp	__tfcbk_MOD_ierrorexp
#define ifunbase	__tfcbk_MOD_ifunbase
#define initmessage	__tfcbk_MOD_initmessage

#define dinfinity	__tfcbk_MOD_dinfinity
#define dnotanumber	__tfcbk_MOD_dnotanumber

#define levelcompile	__tfcbk_MOD_levelcompile
#define levele	__tfcbk_MOD_levele
#define levelp	__tfcbk_MOD_levelp
#define lgeneration	__tfcbk_MOD_lgeneration
#define ltrace	__tfcbk_MOD_ltrace
#define modethrow	__tfcbk_MOD_modethrow
#define iordless	__tfcbk_MOD_iordless

#endif /* _TFCBK_H_ */
