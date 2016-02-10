#ifndef	_SAD_XLIB_H_
#define	_SAD_XLIB_H_

#include <stdbool.h>

typedef struct {
  bool (*SetNewDisplay)(void *);
} sadxlib_hooks_t;

extern sadxlib_hooks_t sadxlib_hooks;

/* Bridge functions */
extern bool sadXlib_SetNewDisplay(void*);

#endif /* _SAD_XLIB_H_ */
