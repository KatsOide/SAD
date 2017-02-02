#include <sim/sad_f2c.h>
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>

int fseek_(integer4 *lun, integer4 *offset, integer4 *whence) {
  fprintf(stderr,
	  "fseek is invoked with LUN=%d offset=%d whence=%d...aborted!\n",
	  *lun, *offset, *whence);
  exit(EX_IOERR);
}
