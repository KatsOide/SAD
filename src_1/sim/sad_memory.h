#define USE_MMAP_FOR_MALLOC

#ifndef _SAD_MEMORY_H_
#define _SAD_MEMORY_H_

#include <sim/sad_f2c.h>
#include <sys/types.h>
#include <stddef.h>

/* Property table for memory allocation sub-system for SADScript interpreter */
struct LM_TABLE {
  real8 *rlist0;
  ptrdiff_t offset;
  size_t align;
  size_t bits;
};

extern struct LM_TABLE *lm_table;

/* Access macro set for SADScript interpreter memory block */

/* *list macro */
#define rlist(index)	(lm_table->rlist0[(index)])
#define ilist(off, index)	\
	(((integer4*)(lm_table->rlist0 + (index)))[(off) - 1])
#define jlist(off, index)	\
	(((char*)(lm_table->rlist0 + (index)))[(off) - 1])
#define klist(index)	\
  ((integer8*)(lm_table->rlist0+(index)))[0]

/* Heap/Map allocate/free interface for SAD interpreter */
extern void	lminit_(real8*, const integer4*);
extern integer4	lmalloc_(const integer4*, integer4*);
extern void	lmfree_(const integer4*);
extern integer4	mapalloc_(void*, const integer4*, const integer4*, integer4*);
extern integer4	mapfree_(void*);

extern void	lminit8_(real8*, const integer4*);
extern integer8	lmalloc8_(const integer4*, integer4*);
extern void	lmfree8_(const integer8*);
extern integer8	mapalloc8_(void*, const integer4*, const integer4*, integer4*);
extern integer4	mapfree8_(void*);

#endif /* _SAD_MEMORY_H_ */
