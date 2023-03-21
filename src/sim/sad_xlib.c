#include <sim/sad_xlib.h>
#include <sim/sad_f2c.h>
#include <stdbool.h>

sadxlib_hooks_t sadxlib_hooks = {NULL};

/* Bridge functions */
bool sadXlib_SetNewDisplay(void *new) {
  if(sadxlib_hooks.SetNewDisplay == NULL)
    return false;

  return sadxlib_hooks.SetNewDisplay(new);
}

/* End of File */
