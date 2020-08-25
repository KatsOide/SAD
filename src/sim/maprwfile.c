#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <semaphore.h>

int64_t maprwfile_(char *fname, int *fd, int64_t *size, int32_t *irtc) {
  void* src;
  struct stat statbuf;

  if((*fd = open(fname, O_RDWR)) < 0)
    {*irtc=1;
      fprintf(stderr,"%s\n",strerror(errno));
      return 0;}

  /*  fprintf(stderr,"size: %d\n",*size); */

  if(*size > 0)ftruncate(*fd,*size);

  if ( fstat(*fd, &statbuf) < 0 )
    {*irtc=2;
      return 0;}

  if ((src = mmap(0, *size = statbuf.st_size, PROT_READ|PROT_WRITE, MAP_SHARED, *fd, 0))
      == (void*) -1)
    { *irtc=3;
      return 0;}

  *irtc=0;
  return (int64_t) src;

}

int unmap_(void* *addr, int32_t *len, int *fd) {
  int irtc;
  irtc=munmap(*addr,*len);
  close(*fd);
  return(irtc);
}
