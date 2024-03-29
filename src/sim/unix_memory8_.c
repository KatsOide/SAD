#include <sim/sad_memory.h>

#define SAD_OFFSET_T integer8
#define lminit_   lminit8_
#define lmalloc_  lmalloc8_
#define lmfree_   lmfree8_
#define mapalloc_ mapalloc8_
#define mapfree_  mapfree8_
#define mapallocshared_ mapallocshared8_
#define mapallocfixed_ mapallocfixed8_
#define mapallocfile_ mapallocfile8_
#define mapallocfd_ mapallocfd8_
#define unmapfile_ unmapfile8_
#define lenfile_ lenfile8_
#define mapresizefile_ mapresizefile8_

#include "sim/unix_memory_.c"

integer4 getpagesize_(void) {
  return sysconf(_SC_PAGESIZE);
}

