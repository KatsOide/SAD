#ifndef _SAD_TCLTK_H_
#define _SAD_TCLTK_H_

#include <sim/sad_f2c.h>
#include <stdbool.h>

/* TclUpdate task flags */
#define TCL_UPDATE_TASK_BITS	4
#define TCL_UPDATE_TASK_MASK	((1L << TCL_UPDATE_TASK_BITS) - 1)
#define TCL_UPDATE_TASK_FILE	(1L << 3)
#define TCL_UPDATE_TASK_TIMER	(1L << 2)
#define TCL_UPDATE_TASK_WINDOW	(1L << 1)
#define TCL_UPDATE_TASK_IDLE	(1L << 0)

/* sadTk_CreateFileHandler masks */
#define SAD_TK_READABLE		(1L << 1)
#define SAD_TK_WRITABLE		(1L << 2)
#define SAD_TK_EXCEPTION	(1L << 3)

typedef void (sadTk_FileProc)(void*, int);

typedef struct {
  void (*TclUpdate)(int);
  void (*CreateFileHandler)(int, int, sadTk_FileProc*, void*);
  void (*DeleteFileHandler)(int);
  bool PendIO;
} sadtk_hooks_t;

extern sadtk_hooks_t sadtk_hooks;

/* SAD Fortran interfaces */
extern void tfsetpendio_(void);
extern void tfresetpendio_(void);
extern void tftclupdate_(integer4*);

/* Bridge functions */
extern void sadTk_CreateFileHandler(int, int, sadTk_FileProc*, void*);
extern void sadTk_DeleteFileHandler(int);

#endif /* _SAD_TCLTK_H_ */
