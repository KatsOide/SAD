#include <sim/sad_tcltk.h>
#include <sim/sad_api.h>
#include <stdbool.h>

sadtk_hooks_t sadtk_hooks = {NULL, NULL, NULL, false};

/* SAD Fortran interfaces */
void tfsetpendio_(void) {
  sadtk_hooks.PendIO = true;
}

void tfresetpendio_(void) {
  sadtk_hooks.PendIO = false;
}

void tftclupdate_(integer4 *mode0) {
  static integer8 itfInterruptMask = 0;
  int mode, mask;

  if(itfInterruptMask == 0) {
    itfInterruptMask = ktfsymbolz("System`FFS$InterruptMask") - 4;
  }

  mode = *mode0 & TCL_UPDATE_TASK_MASK;
  mask = rlist(itfInterruptMask);

  if(mode != 0) {
    mode &= ~mask;		/* Mask by FFS$InterruptMask */
    if(sadtk_hooks.PendIO)	/* Mask except TASK_IDLE if PendIO is true */
      mode &= TCL_UPDATE_TASK_IDLE;

    if(mode == 0) return;	/* Skip update if all event is masked */
  }

  if(sadtk_hooks.TclUpdate)
    sadtk_hooks.TclUpdate(mode);
}

/* Bridge functions */
void sadTk_CreateFileHandler(int fd, int mask,
			     sadTk_FileProc *proc, void *data) {
  if(sadtk_hooks.CreateFileHandler != NULL)
    sadtk_hooks.CreateFileHandler(fd, mask, proc, data);
}

void sadTk_DeleteFileHandler(int fd) {
  if(sadtk_hooks.DeleteFileHandler != NULL)
    sadtk_hooks.DeleteFileHandler(fd);
}

/* End of File */
