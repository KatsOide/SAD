#include <sys/types.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <fcntl.h>
#include <sys/stat.h>

int64_t maprwfile_(char *fname, int *fd, int64_t *size, int32_t *irtc) {
  caddr_t src;
  struct stat statbuf;

  if((*fd = open(fname, O_RDONLY)) < 0)
    {*irtc=1;
      return 0;}

  if ( fstat(*fd, &statbuf) < 0 )
    {*irtc=2;
      return 0;}

  if ((src = mmap(0, *size = statbuf.st_size, PROT_READ|PROT_WRITE, MAP_SHARED, *fd, 0))
      == (caddr_t) -1)
    { *irtc=3;
      return 0;}

  *irtc=0;
  return (int64_t) src;

}
