#include <sim/sad_memory.h>
#include <sim/sad_f2c.h>
#include <sim/sad_api.h>

#include <sys/types.h>
#include <sys/mman.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/syscall.h>
#include <fcntl.h>
#include <sys/stat.h>

#if defined(UINT64_MAX) && !defined(PRIuSIZE_T)
#if UINT_FAST64_MAX == SIZE_MAX
#define	PRIuSIZE_T	PRIuFAST64
#endif
#endif
#if defined(UINT32_MAX) && !defined(PRIuSIZE_T)
#if UINT_FAST32_MAX == SIZE_MAX
#define	PRIuSIZE_T	PRIuFAST32
#endif
#endif
#if defined(UINT16_MAX) && !defined(PRIuSIZE_T)
#if UINT_FAST16_MAX == SIZE_MAX
#define	PRIuSIZE_T	PRIuFAST16
#endif
#endif
#if defined(UINT8_MAX) && !defined(PRIuSIZE_T)
#if UINT_FAST8_MAX == SIZE_MAX
#define	PRIuSIZE_T	PRIuFAST8
#endif
#endif
#if !defined(PRIuSIZE_T)
#define	PRIuSIZE_T	"u"
#endif

#ifndef	DEBUG_MEMORY
#define	DEBUG_MEMORY		-1
#endif

#ifndef CHECK_LMALLOC_SANITY
#define CHECK_LMALLOC_SANITY	1
#endif

#ifndef CHECK_MAPALLOC_SANITY
#define CHECK_MAPALLOC_SANITY	1
#endif

#ifndef	MAP_ANONYMOUS
#ifdef	MAP_ANON
#define	MAP_ANONYMOUS	MAP_ANON
#else
#define	MAP_ANONYMOUS	0
#endif
#endif

#ifndef	MAP_VARIABLE
#define	MAP_VARIABLE	0
#endif

#ifndef	MAP_FAILED
#define	MAP_FAILED	((void *)-1)
#endif

#ifdef	TRY_SAD_MAP_ADDR32
#if	defined(MAP_32BIT)
#define	SAD_MAP_ADDR32	MAP_32BIT
#elif	defined(MAP_ADDR32)
#define	SAD_MAP_ADDR32	MAP_ADDR32
#else
#define	SAD_MAP_ADDR32	0
#endif
#else	/* !TRY_SAD_MAP_ADDR32 */
#define	SAD_MAP_ADDR32	0
#endif

#if	DEBUG_MEMORY - 0 >= 3
#define DB1(x)	x
#define DB2(x)	x
#define DB3(x)	x
#elif	DEBUG_MEMORY - 0 >= 2
#define DB1(x)	x
#define DB2(x)	x
#define DB3(x)
#elif	DEBUG_MEMORY - 0 >= 1
#define DB1(x)	x
#define DB2(x)
#define DB3(x)
#else
#define DB1(x)
#define DB2(x)
#define DB3(x)
#endif

#ifdef	USE_ELECTRICFENCE
extern int EF_ALIGNMENT;
extern int EF_PROTECT_FREE;
#endif /* USE_ELECTRICFENCE */

/* Memory allocation sub-system for SADScript interpreter */
#ifndef SAD_OFFSET_T
#define SAD_OFFSET_T integer4
struct LM_TABLE *lm_table = NULL;
#endif

/* allocator initializer MUST be called with rlist(0) at first */
/* call lminit(rlist(0), alignment_size) */
void lminit_(real8 *rlist0, const integer4 *align, char *pname0, integer4 *lpname0, 
             integer4 *idtype0, integer8 *idval0) {
  size_t bits;

  if(lm_table != NULL) {
    fprintf(stderr, "lminit: recursive called!\n");
    exit(1);
  }

  /* Check align is 2^N */
  if(*align < 1 || ((*align - 1) & *align) != 0) {
    fprintf(stderr, "lminit: non-2^N alignment is required!\n");
    exit(1);
  } else {
    bits = 0;
    while(((ptrdiff_t)1 << bits) < *align) bits += 1;
    if(*align != ((ptrdiff_t)1 << bits)) {
      fprintf(stderr, "lminit: bit shift sanity check is failed!\n");
      exit(1);
    }
  }

#ifdef	USE_ELECTRICFENCE
  EF_ALIGNMENT = sizeof(real8);
#ifdef	USE_ELECTRICFENCE_FREECHECK
  EF_PROTECT_FREE = 1;
#endif /* USE_ELECTRICFENCE_FREECHECK */
#endif /* USE_ELECTRICFENCE */

  lm_table = malloc(sizeof(struct LM_TABLE));
  if(lm_table == NULL) {
    fprintf(stderr, "lminit: system table allocation is failed!\n");
    exit(1);
  }

  lm_table->rlist0 = rlist0; /* Copy reference pointer */
  lm_table->align  = *align; /* Copy alignment size */
  lm_table->bits   = sizeof(ptrdiff_t) * CHAR_BIT - bits;
  lm_table->offset = ((const char *)rlist0 - (const char *)lm_table);
  lm_table->offset &= (lm_table->align - 1);

  lm_table->pname0 = pname0; /* Copy reference pointer */
  lm_table->lpname0 = lpname0; /* Copy reference pointer */
  lm_table->idtype0 = idtype0; /* Copy reference pointer */
  lm_table->idval0 = idval0; /* Copy reference pointer */

  DB1(fprintf(stderr, "lminit: rlist0=%p align=%"PRIuSIZE_T" bits=%"PRIuSIZE_T" offset=%"PRIuPTR"\n",
	      (void*)lm_table->rlist0, lm_table->align,
	      lm_table->bits, lm_table->offset));
}

/* allocate memory chunk for *list */
/* rlist_offset = lmalloc(size, irtc) */
#ifdef	USE_MMAP_FOR_MALLOC
SAD_OFFSET_T lmalloc_(const integer4 *size0, integer4 *irtc) {
  integer4 unit;

  unit = lm_table->align;
  printf("lmalloc %d\n",*size0);
  return mapalloc_(lm_table->rlist0, size0, &unit, irtc);
}
#else
SAD_OFFSET_T lmalloc_(const integer4 *size0, integer4 *irtc) {
  char *ptr, *ptr0;
  size_t size = (*size0) * lm_table->align + lm_table->offset;
  ptrdiff_t offset;
  SAD_OFFSET_T roffset;

  DB2(fprintf(stderr, "lmalloc%"PRIuSIZE_T": chunk allocation is requested"
	      "[%"PRIuSIZE_T" * %d = %"PRIuSIZE_T"]\n",
	      sizeof(SAD_OFFSET_T), sizeof(real8),
	      *size0, (*size0) * sizeof(real8)));

  *irtc = 1;
  ptr = malloc(size);
  if(ptr == NULL) return 0;

#if CHECK_LMALLOC_SANITY > 0
  if(((ptr - (char *)lm_table) & (lm_table->align - 1)) != 0) {
    free(ptr);
    fprintf(stderr, "lmalloc%"PRIuSIZE_T": broken malloced chunk alignment\n",
	    sizeof(SAD_OFFSET_T));
    return 0;
  }
#endif /* CHECK_LMALLOC_SANITY */

  ptr0 = ptr;
  ptr += lm_table->offset;
  offset = (real8 *)ptr - lm_table->rlist0;

  /* Try negative offset hack */
  if(offset < 0) {
    offset += ((ptrdiff_t)1 << lm_table->bits);
    if(lm_table->rlist0 + offset != (real8 *)ptr) { /* Point different addr. */
      fprintf(stderr, "lmalloc%"PRIuSIZE_T": negative offset hack is failed!\n",
	      sizeof(SAD_OFFSET_T));
      offset -= ((ptrdiff_t)1 << lm_table->bits);
    }
  }

  roffset = (SAD_OFFSET_T)offset;

  if(roffset >= 0 && (ptrdiff_t)roffset == offset) {
    DB2(fprintf(stderr, "lmalloc%"PRIuSIZE_T": allocate chunk offset=0x%tx (heap=%p)\n",
		sizeof(SAD_OFFSET_T), offset, ptr0));
    *irtc = 0;
    return roffset;
  }

  fprintf(stderr,
	  "lmalloc%"PRIuSIZE_T": allocated chunk offset is out of range\n"
	  "          base=%p heap=%p offset=0x%tx\n",
	  sizeof(SAD_OFFSET_T), (void*)lm_table->rlist0, ptr0, offset);
  free(ptr0);

  return 0;
}
#endif

#ifdef	USE_MMAP_FOR_MALLOC
void lmfree_(const SAD_OFFSET_T *offset) {
  char *ptr = (char *)(lm_table->rlist0 + *offset);

  mapfree_(ptr);
}
#else
void lmfree_(const SAD_OFFSET_T *offset) {
  char *ptr;

  ptr = (char *)(lm_table->rlist0 + *offset);
  ptr -= lm_table->offset;

  DB2(fprintf(stderr, "lmfree%"PRIuSIZE_T": release chunk offset=0x%tx (heap=%p)\n",
	      sizeof(SAD_OFFSET_T), (ptrdiff_t)*offset, ptr));
  free(ptr);
}
#endif

/* memory map sub-system for SAD interpreter */
/*
 * mapped memory structure
 * |size_t size|align salt[1, unit]|usize * unit octet|
 *               -keep salt value- |<-base[offset]
 */

/* mapalloc(base, size, unit, irtc)
 *     call gate for Unix mmap system call
 * Arguments:
 *     base    base address for returned offset and memory alignment
 *     size    map size by `unit'
 *     unit    unit size[octet](MUST be 2^N)
 *     irtc    0:         return with error code if mmap failed
 *             none-zero: abort if mmap failed
 * Results:
 *     irtc    0:        Succeed
 *             Non-zero: Error number
 *     mapalloc          offset from base counted by `unit' octet
 */
SAD_OFFSET_T mapalloc_(void *base, const integer4 *usize,
		  const integer4 *unit, integer4 *irtc) {
  int errcode;
  char align, *map;
  size_t size, bits, *map0;
  ptrdiff_t offset;
  SAD_OFFSET_T roffset;

  DB3(fprintf(stderr, "mapalloc%"PRIuSIZE_T": memory map allocation is requested"
	      "[base=%p size=(%d * %d = %d)]\n",
	      sizeof(SAD_OFFSET_T), base, *unit, *usize, (*unit) * (*usize)));

  /* Check unit is 2^N */
  if(*unit < 1 || ((*unit - 1) & *unit) != 0) {
    fprintf(stderr, "mapalloc%" PRIuSIZE_T ": non-2^N unit is required!\n",
	    sizeof(SAD_OFFSET_T));
    exit(1);
  } else {
    bits = 0;
    while(((ptrdiff_t)1 << bits) < *unit) bits += 1;
    if(*unit != ((ptrdiff_t)1 << bits)) {
      fprintf(stderr, "mapalloc%"PRIuSIZE_T": bit shift sanity check is failed!\n",
	      sizeof(SAD_OFFSET_T));
      exit(1);
    }
    bits = sizeof(ptrdiff_t) * CHAR_BIT - bits;
  }

  size = (*usize) * (*unit);
  size += sizeof(size_t); /* size store */
  size += ((*unit > sizeof(align)) ? *unit : sizeof(align));


  map0 = MAP_FAILED;
  while(true) {
#if	defined(USE_ADDR_HINT_IN_MMAP) && !defined(USE_MMAP_FOR_MALLOC)
    map0 = mmap(base, size, PROT_READ | PROT_WRITE, SAD_MAP_ADDR32 |
		MAP_ANONYMOUS | MAP_SHARED | MAP_VARIABLE, -1, 0);
    if(map0 != MAP_FAILED) break;
    DB1(fprintf(stderr, "mapalloc%"PRIuSIZE_T": failed to apply address hint to map"
		"[base=%p size=%"PRIuSIZE_T"]\n",
		sizeof(SAD_OFFSET_T), base, size));
#endif
    map0 = mmap(NULL, size, PROT_READ | PROT_WRITE, SAD_MAP_ADDR32 |
		MAP_ANONYMOUS | MAP_SHARED | MAP_VARIABLE, -1, 0);

    break;	/* Don't delete this break line */
  }

  if(map0 == MAP_FAILED) {
    errcode = errno;
    fprintf(stderr, "mapalloc%"PRIuSIZE_T": mmap failed because of %s\n",
	    sizeof(SAD_OFFSET_T), strerror(errno));
    if(*irtc == 0) {
      *irtc = errcode;
      return(0);
    } else {
      exit(1);
    }
  }

  map0[0] = size; /* fill size header */
  map = (char*)&(map0[1]);
  offset = ((const char*)map - (const char*)base);
  
  align = offset % *unit;
  if(align < 0) {
    align += *unit;
#if CHECK_MAPALLOC_SANITY > 0
    if(align < 0) {
      fprintf(stderr, "mapalloc%"PRIuSIZE_T": align sanity check failed!\n",
	      sizeof(SAD_OFFSET_T));
      exit(1);
    }
#endif /* CHECK_MAPALLOC_SANITY */
  }

  align = *unit - align;

  map += align; /* align with base pointer by *unit octect */
  offset = (const char *)map - (const char *)base;
  offset /= *unit;
  map[-1] = align; /* keep align shift */

  /* negative offset hack */
  if(offset < 0) {
    offset += ((ptrdiff_t)1 << bits);
    /* Please implement offset hack sanity check! */
  }

  roffset = (SAD_OFFSET_T)offset;

  if(roffset >= 0 && (ptrdiff_t)roffset == offset) {
    DB2(fprintf(stderr, "mapalloc%"PRIuSIZE_T": "
		"allocate %"PRIuSIZE_T"bytes memory map at %p"
		"[offset=0x%tx align=%d]\n",
		sizeof(SAD_OFFSET_T), size, (void*)map0, offset, align));
    *irtc = 0;
    return roffset;
  }

  fprintf(stderr,
	  "mapalloc%"PRIuSIZE_T": allocated chunk offset is out of range!\n"
	  "          base=%p map=%p unit=%d offset=0x%tx\n",
	  sizeof(SAD_OFFSET_T), base, (void*)map0, *unit, offset);
  munmap(map0, size);
  if(*irtc != 0) exit(1);
  *irtc = 1;
  return 0;
}

/* mapfree(base[offset])
 *    call gate for Unix munmap
 * Argument:
 *   base[offset]  top address of allocated memory map by mapalloc
 */
integer4 mapfree_(void *ptr) {
  char *map;
  size_t *map0, size;
  int success;

  DB3(fprintf(stderr, "mapfree%" PRIuSIZE_T ": memory map free is requested[addr=%p]\n",
	      sizeof(SAD_OFFSET_T), ptr));

  map   = (char *)ptr;
  map  -= map[-1];		/* Unwind alignment */
  map0  = (size_t *)map;
  map0 -= 1;			/* Unwind sizeof(size_t) */
  size  = map0[0];

  success = munmap(map0, size);
  if(success == 0) {
    DB2(fprintf(stderr, "mapfree%"PRIuSIZE_T": release %"PRIuSIZE_T"bytes memory map from %p\n",
		sizeof(SAD_OFFSET_T), size, (void*)map0));
  } else {
    fprintf(stderr, "mapfree%" PRIuSIZE_T ": munmap failed because of %s\n",
		sizeof(SAD_OFFSET_T), strerror(errno));
  }
  return success;
}

SAD_OFFSET_T mapallocfixed_(void *base, const integer8 *usize,
		  const integer4 *unit, integer4 *irtc) {
  size_t size, *map0;
  /*  integer8 map; */

  size = (*usize) * (*unit);
/*  if(munmap(base, size)){
    printf("mapallocshared-unmap %d,%d\n",base,size);
    fprintf(stderr,"%s\n",strerror(errno));
    *irtc=-1;} */
  map0 = MAP_FAILED;
  map0 = mmap(base, size, PROT_READ | PROT_WRITE, 
              MAP_ANONYMOUS | MAP_PRIVATE | MAP_FIXED, -1, 0);
  /*  map=syscall(SYS_mmap,base, size, PROT_READ | PROT_WRITE, 
      MAP_ANONYMOUS | MAP_SHARED | MAP_FIXED, -1, 0); 
  if(map < 0){ */
  if(map0 == MAP_FAILED) {
    /*    printf("mapallocfixed %llu,%llu\n",base,size);
          fprintf(stderr,"%s\n",strerror(errno)); */
    *irtc=-1;}
  else{
    *irtc=0;};
  return((integer8) map0);

}

void mapallocshared_(void *base, const integer8 *usize,
		  const integer4 *unit, integer4 *irtc) {
  char align;
  size_t size, *map0;
  /*  integer8 map; */

  size = (*usize) * (*unit);
  /*  if(munmap(base, size)){
    printf("mapallocshared-unmap %d,%d\n",base,size);
    fprintf(stderr,"%s\n",strerror(errno));
    *irtc=-1;} */
  map0 = MAP_FAILED;
  map0 = mmap(base, size, PROT_READ | PROT_WRITE, 
              MAP_ANONYMOUS | MAP_SHARED | MAP_FIXED, -1, 0);
  /*  map=syscall(SYS_mmap,base, size, PROT_READ | PROT_WRITE, 
      MAP_ANONYMOUS | MAP_SHARED | MAP_FIXED, -1, 0);
  if(map < 0){ */
  if(map0 == MAP_FAILED){
    printf("mapallocshared %llu,%llu\n",base,size);
    fprintf(stderr,"%s\n",strerror(errno));
    *irtc=-1;}
  else{
    *irtc=0;};

}

#ifndef caddr_t
typedef char *caddr_t;
#endif

integer8 mapallocfile_(char *fname, int *fd, integer8 *size, integer4 *irtc) {
  char *src, *rpath, *expt;
  struct stat statbuf;

  expt=expand_tilde(fname);
  rpath = realpath(expt,NULL);
  free(expt);

  if((*fd = open(rpath, O_RDONLY)) < 0)
    {free(rpath);
     *irtc=1;
      return 0;}

  if ( fstat(*fd, &statbuf) < 0 )
    {free(rpath);
      *irtc=2;
      return 0;}

  if ((src = mmap(0, *size = statbuf.st_size, PROT_READ, MAP_SHARED, *fd, 0))
      == (caddr_t) -1)
    {free(rpath);
      *irtc=3;
      return 0;}

  free(rpath);
  *irtc=0;
  return (integer8) src;

}

integer8 mapallocfd_(int *fd, integer8 *size, integer4 *irtc) {
  char *src;
  struct stat statbuf;

  if ( fstat(*fd, &statbuf) < 0 )
    {
      *irtc=2;
      return 0;}

  if ((src = mmap(0, *size = statbuf.st_size, PROT_READ, MAP_SHARED, *fd, 0))
      == (caddr_t) -1)
    {
      *irtc=3;
      return 0;}

  *irtc=0;
  return (integer8) src;

}

integer4 unmapfile_(char *map,integer8 *size){
  return munmap(map,*size);
}

integer8 lenfile_(int* fd){
  struct stat statf;
  if ( fstat(*fd, &statf) < 0 )
    {
      return -1;}
return statf.st_size;}

integer8 mapresizefile_(char *map, int* fd, integer8* size0, integer8 *size){
  char *src;
  munmap(map,*size0);
  if((src = mmap(0,*size,PROT_READ,MAP_SHARED,*fd,0)) ==(caddr_t) -1)
    return -1;
  return (integer8) src;
}

/* End of File */
