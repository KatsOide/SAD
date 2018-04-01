#ifndef _MACCBK_H_
#define _MACCBK_H_

/* sim for inc/MACCBK.inc */
#include <sim/sad_f2c.h>

#define	MAXPNAME	32
/*extern integer4 __maccbk_MOD_MAXPNAME;
  #define MAXPNAME __maccbk_MOD_MAXPNAME*/
#define	HTMAX		65536

/* idlist common block */
/*
typedef struct {
  char __pname[HTMAX][MAXPNAME];
  integer4 __idtype[HTMAX];
  integer8 __idval[HTMAX];
} idlist_t;

extern idlist_t idlist_;
*/

/* idlist common block macro 
#define	pname(i)	(idlist_.__pname[(i) - 1])
#define	idtype(i)	(idlist_.__idtype[(i) - 1])
#define	idval(i)	(idlist_.__idval[(i) - 1])
*/
/*
extern char __maccbk_MOD_pname[HTMAX][MAXPNAME];
extern integer4 __maccbk_MOD_idtype[HTMAX];
extern integer8 __maccbk_MOD_idval[HTMAX];
#define	pname(i)	__maccbk_MOD_pname[(i) - 1]
#define	idtype(i)	__maccbk_MOD_idtype[(i) - 1]
#define	idval(i)	__maccbk_MOD_idval[(i) - 1]
*/

#endif /* _MACCBK_H_ */
