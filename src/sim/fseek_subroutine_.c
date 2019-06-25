#include <sim/sad_f2c.h>

extern integer4 fseek_subroutine_(integer4*, integer4*, integer4*);

integer4 fseek_(integer4 *lun, integer4 *offset, integer4 *whence) {
  return fseek_subroutine_(lun, offset, whence);
}

/* End of File */
