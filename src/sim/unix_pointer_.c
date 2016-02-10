#include <sim/sad_f2c.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/* pointer2double(argv)
 * Return: argv address divide sizeof(void*)
 */
real8 pointer2double(void *var) {
  ptrdiff_t address = (ptrdiff_t)var;
  real8 index;

  index = (double)(address / sizeof(void*));
  if(address % sizeof(void*) != 0 || (((ptrdiff_t)index) * sizeof(void*)) != address) {
    fprintf(stderr, "pointer2double: out of range[%p]\n", var);
    abort();
    exit(1);
  }

  return index;
}

void* double2pointer(real8 var) {
  return (void*)(((ptrdiff_t)var) * sizeof(void*));
}

/* decloc(argv)
 * Result: argv [ address / sizeof(void*) ] in interger4: like a %LOC() Intrinsic
 */
integer4 decloc_(void *var) {
  ptrdiff_t address = (ptrdiff_t)var;
  integer4 decloc;

  decloc = (integer8)address / sizeof(void*);
  if(sizeof(ptrdiff_t) != sizeof(integer4)){
    if((integer8)decloc * sizeof(void*) != (integer8) address) {
      fprintf(stderr, "decloc: out of range  %x %x\n",(integer8)decloc * sizeof(void*),(integer8) address);
      exit(1);
  }}
  return decloc;
}

/* setfnp(table, function pointer) */
void setfnp_(integer4 *table, void (*func)(integer4 *)) {
  ptrdiff_t address = (ptrdiff_t)func;
  integer4 decloc;

  decloc = (integer8)address / sizeof(void*);
  if(sizeof(ptrdiff_t) != sizeof(integer4)){
    if((integer8)decloc * sizeof(void*) != (integer8) address) {
      fprintf(stderr, "setfnp: out of range  %x %x\n",(integer8)decloc * sizeof(void*),(integer8) address);
      exit(1);
  }}

  *table = decloc;
}

/* funcall1(func, argv) */
void funcall1_(integer4 *func, integer4 *arg1) {
  ptrdiff_t address = ((ptrdiff_t)*func)*sizeof(void*);

  /* Cast integer *func as 1 argument function: void (*)(integer4 *) */
  ((void (*)(ptrdiff_t *)) address)((ptrdiff_t *) arg1);
}

