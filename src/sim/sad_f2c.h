#ifndef _SAD_F2C_H_
#define _SAD_F2C_H_

#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>

typedef float   real4;
typedef double  real8;

typedef int8_t  integer1;
typedef int16_t integer2;
typedef int32_t integer4;
typedef int64_t integer8;

typedef int8_t  logical1;
typedef int16_t logical2;
typedef int32_t logical4;
typedef int64_t logical8;

typedef unsigned int ftnlen;
typedef       char*       character;
typedef const char* const_character;

/* Default length type */
typedef real8   real;

#ifndef INTEGER_SIZE
#define INTEGER_SIZE 4
#endif
#if INTEGER_SIZE == 2
typedef integer2 integer;
typedef logical2 logical;
#elif INTEGER_SIZE == 4
typedef integer4 integer;
typedef logical4 logical;
#elif INTEGER_SIZE == 8
typedef integer8 integer;
typedef logical8 logical;
#else
#error "Invalid default Integer size"
#endif

/* Fortran %LOC SIM Prototypes */
extern integer4 decloc_(void *);

/* Fortran intrinsic for C */
extern integer len_trim(const_character, ftnlen);

/* Fortran Library API Prototypes */
extern integer4 lnblnk_(const_character, ftnlen);
extern integer4 getfd_(integer4 *);
extern integer fork_worker_(void);
extern void tfsavesharedmap_(void);
extern integer __tfshare_MOD_itffork(void);

#endif /* _SAD_F2C_H_ */
