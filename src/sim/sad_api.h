#ifndef _SAD_API_H_
#define _SAD_API_H_

#include <sim/sad_f2c.h>
#include <sim/TFSTK.h>

/* Helper to store object pointer into Real*8 */
real8 pointer2double(void*);
void* double2pointer(real8);

/* Provide asprintf() if system does not have it */
/* Note: asprintf() is not ISO C99/POSIX function */
#ifndef	HAVE_ASPRINTF
#include <stdio.h>
#define	asprintf	sad_asprintf
extern int asprintf(char**, const char*, ...);
#endif	/* HAVE_ASPRINTF */

/* Require runtime null-character termination for using SAD string? */
#define SAD_REQUIRE_STRING_TERMINATION	0

/*
 * Offset between SAD epoch and Unix time epoch in seconds
 * SAD  epoch: 1900/01/01/ 00:00:00 JST
 * Unix epoch: 1970/01/01/ 00:00:00 UTC
 */
#define SAD_EPOCH_OFFSET	2209021200.0

/* Return real8 object as True/False */
#define SAD_BOOLEAN(cond)	((cond) ? 1. : 0.)

/* Typedef for SADScript function pointer */
typedef int (*sad_func_t)(integer4*, integer4*, integer4*, real8*, integer4*);
typedef int (*sad_func_t_8)(integer4*, integer8*, integer4*);

/* Prototype declaration for SADScript function */
#define DECLARE_SAD_FUNC(symbol) \
  int symbol(integer4*, integer4*, integer4*, real8*, integer4*)
#define DECLARE_SAD_FUNC8(symbol) \
  int symbol(integer4*, integer8*, integer4*)

/* Portable bit field operator */
typedef unsigned int	bit_field_t;
#define __BF_WIDTH	(sizeof(bit_field_t) * 8)
#define __BF_LENGTH(n)	(((n) + __BF_WIDTH - 1) / __BF_WIDTH)
#define __BF_MASK(x)	((bit_field_t)1 << ((x) % __BF_WIDTH))

#define BF_ISSET(n, p)	((p[(n) / __BF_WIDTH] & __BF_MASK(n)) != 0)
#define BF_CLR(n, p)	(p[(n) / __BF_WIDTH] &= ~__BF_MASK(n))
#define BF_SET(n, p)	(p[(n) / __BF_WIDTH] |=  __BF_MASK(n))
#define BF_ZERO(n, p)	do { \
    bit_field_t *_p; \
    size_t _n; \
    _p = (p); \
    _n = __BF_LENGTH(n); \
    while(_n > 0) _p[--_n] = 0; \
  } while(0)
#define BF_ALLOC(n, p)	(((p) = malloc(__BF_LENGTH(n) * sizeof(bit_field_t))) != NULL)
#define BF_FREE(p)	(free(p))

#define ktfrealq(k) ((ktrmask & (k)) != ktfnr)
#define ktfnonrealq(k) ((ktrmask & (k)) == ktfnr)
#define ktflistq(k) ((ktfmask & (k)) == ktflist)
#define ktfnonlistq(k) ((ktfmask & (k)) != ktflist)
#define ktfsymbolq(k) ((ktfmask & (k)) == ktfsymbol)
#define ktfstringq(k) ((ktfmask & (k)) == ktfstring)
#define ktfnonstringq(k) ((ktfmask & (k)) != ktfstring)
#define ktfoperq(k) ((ktfmask & (k)) == ktfoper)
#define ktfaddr(k) (ktamask & (k))
#define ktfreallistq(k) ((lnonreallist & ilist(2,(k)-3)) == 0)
#define ktfnonreallistq(k) ((lnonreallist & ilist(2,(k)-3)) != 0)
#define kfromr(x) (*((integer8 *) &x))
#define rfromk(k) (*((real8 *) &k))

static real8 r_true = 1.0;

/* Proxy functions between module and global scope */

/* Fortran based Internal API Prototypes */
extern real8    atof_(const_character, integer4*, ftnlen);
extern integer4 hsrch_(const_character, ftnlen);
extern integer4 hsrchz_(const_character, ftnlen);
extern integer4 igetgl_(const_character, ftnlen);
extern real8    rgetgl1_(const_character, ftnlen);
extern void capita_(character, ftnlen);
extern void tfreadbuf_(integer4*, integer4*, integer4*);
extern void __tfrbuf_MOD_trbinit(integer4*, integer4*);
extern void gettok_(character, integer4*, integer4*, real8*, integer4*, ftnlen);
extern void __tfmem_MOD_tfree(integer8*);
extern integer4 itfdownlevel_(void);
extern logical8 tfruleqk_(integer8*);
extern logical4 tfsamesymbolq_(integer4*, integer4*);
extern logical4 tfsamesymbolqk_(integer8*, integer8*);
extern integer4 itopenbuf_(integer4*, integer4*);
extern integer4 itfsyserr_(integer4*);
extern integer8 ktfcopy_(integer8*);
extern integer8 kfromr_(real8*);
extern real8 rfromk_(integer8*);
extern integer4 itfmessage_(integer4*, const_character, const_character,
			    ftnlen, ftnlen);
extern void tfreseterror_();
extern integer4 itfsymbolc_(const_character, integer4*, integer4*, ftnlen);
extern integer8 ktfsymbolz_(const_character, integer4*, ftnlen);
extern integer8 ktfsymbolf_(const_character, integer4*, logical4*, ftnlen);
extern integer8 ktsalocb_(integer4*, const_character, integer4*, ftnlen);
extern integer4 italoc_(integer4*);

extern integer8 ktaaloc_(integer4*,integer4*);
extern integer8 ktavaloc_(integer4*, integer4*);
extern integer8 ktadaloc_(integer4*, integer4*);
extern integer8 ktfmaloc_(integer8*, integer4*, integer4*,
			  logical4*, logical4*, integer4*);
extern integer4 itfunaloc_(const_character, integer4*, integer4*,
			   integer4*, integer4*, integer4*, ftnlen);
extern integer8 ktfmakelist_(integer4*);
extern void tfsetlist_(integer8*, integer8*, integer4*);
extern void tfgetstr_(character, ftnlen, integer8*, integer4*);
extern void tfevalb_(const_character, integer4, integer8*, integer4*);
extern void tfevalc_(const_character, integer4);
extern void tfdeval_(integer4*, integer8*, integer8*, 
		     integer4*, logical4 *, integer4*, integer4*);
extern void tflocal_(integer8*);
extern void tflocal1_(integer8*);
extern void tfmakerulestk_(integer8*, integer8*);

/* Extended Allocation API Prototypes */
extern integer4 itfgetoptionstk(integer4, const char**);
extern int tfinitstk(tfstk_t*, integer4);

/* Internal API Wrapper Prototypes */
extern void tfreadbuf(integer4, integer4*, integer4*);
extern void trbinit(integer4, integer4*);
extern real8 rgetgl1(const char*);
extern integer4 itfsyserr(integer4);
extern integer8 ktfcopy(integer8);
extern integer4 itfmessage(integer4, const char*, const char*);
extern integer4 itfsymbolc(const char*, integer4);
extern integer8 ktfsymbolz(const char*);
extern integer8 ktfsymbolf(const char*, logical4);
extern integer4 itfsymbol(const char*, int);
extern void tfree(integer8);
extern integer8 ktsalocb(integer4 , const char*);
extern integer8 ktsalocbl(integer4, const char*, integer4);
extern integer4 italoc(integer4);
extern integer8 ktaaloc(integer4, integer4);
extern integer8 ktavaloc(integer4, integer4);
extern integer8 ktadaloc(integer4, integer4);
extern integer8 ktfmaloc(integer8, integer4*, integer4*,
			 logical4, logical4, integer4*);
extern integer8 ktfm2l(real8*, integer4, integer4, integer4, logical4);
#define	itfdownlevel	itfdownlevel_
extern integer8 ktfmakelist(integer4);
extern void tfsetlist(integer8, integer8, integer4);
extern void tfmakerulestk(integer8, integer8);
extern void tfevals(const char*, integer8*, integer4*);
extern void tfevalc(const char*);
extern void tfdeval(integer4, integer8, integer8*, 
		    integer4, logical4,
		    integer4*, integer4*);
extern void tflocal(integer8);
extern void tflocal1(integer8);

/* EPICS Channel Access Value Callback API Prototypes */
extern void tfcavaluecb_(real8*, integer4*, integer4*, real8*,
			 integer4*, integer4*, integer8*);
extern void tfepicsvaluecb_(real8*, integer4*, integer4*, real8*,
			    integer4*, integer4*, integer8*);
extern void tfepicsconstatcb_(real8*, integer4*);

/* Dynamic function registration API */
extern int dlfunaloc(const char*, sad_func_t, int, const int*, const int*, int);
extern int dlfunaloc8(const char*, sad_func_t_8, int, const int*, const int*, int);

/* Expand `~' prefix by shell rule and return new allocated string */
extern char* expand_tilde(const char *);

#endif /* _SAD_API_H_ */
