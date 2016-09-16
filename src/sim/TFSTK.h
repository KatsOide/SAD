#ifndef _TFSTK_H_
#define _TFSTK_H_

/* sim for inc/TFSTK.inc */
#include <sim/sad_f2c.h>
#include <sim/sad_memory.h>

/* tfstk common block */
typedef struct {
  integer4 __mstk;
  integer4 __isp;
  integer4 __ivstkoffset;
  integer4 __ipurefp;
  integer4 __napuref;
  integer4 __isporg;
  integer8 __ispbase;
} tfstk_t;

extern tfstk_t tfstk_;
extern integer4		__tfstk_MOD_mstk;
extern integer4		__tfstk_MOD_isp;
extern integer4		__tfstk_MOD_ivstkoffset;
extern integer4		__tfstk_MOD_ipurefp;
extern integer4		__tfstk_MOD_napuref;
extern integer4		__tfstk_MOD_isporg;
extern integer8		__tfstk_MOD_ispbase;


/* *stk macro */
#define   vstk(     index)	rlist(     index+ispbase)
#define  ivstk(off, index)	ilist(off, index+ispbase)
#define rtastk(     index)	rlist(     index+ispbase)
#define itastk(off, index)	ilist(off, index+ispbase)
#define jtastk(off, index)	\
	(((integer2*)(lm_table->rlist0 + (index+ispbase)))[(off) - 1])
#define ktastk(     index)      klist(     index+ispbase)

/* tfstk common block macro */
/*
#define mstk		(tfstk_.__mstk)
#define isp		(tfstk_.__isp)
#define ivstkoffset	(tfstk_.__ivstkoffset)
#define ipurefp		(tfstk_.__ipurefp)
#define napuref		(tfstk_.__napuref)
#define isporg		(tfstk_.__isporg)
*/
#define mstk		__tfstk_MOD_mstk
#define isp		__tfstk_MOD_isp
#define ivstkoffset	__tfstk_MOD_ivstkoffset
#define ipurefp		__tfstk_MOD_ipurefp
#define napuref		__tfstk_MOD_napuref
#define isporg		__tfstk_MOD_isporg
#define ispbase		__tfstk_MOD_ispbase

#endif /* _TFSTK_H_ */
