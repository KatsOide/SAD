/* MD part for dynamic object loader implemented by dlopen(3) */

#include <sim/dynl.h>

#include <sysexits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <dlfcn.h>

#ifndef RTLD_SELF
#define RTLD_SELF NULL
#endif

static bool dynl_initialized = false;

int dynl_error_report_level = 0;

bool dynl_provide(void) {
  return true;
}

bool dynl_init(void) {
  if(!dynl_initialized) {
    dynl_initialized = true;
    return true;
  }

  return true;
}

void *dynl_link(const char *fname, int mode) {
  void *handle;
  int md_mode = 0;

  switch(mode & DYNL_MODEMASK) {
  case DYNL_LAZY:
    md_mode = RTLD_LAZY;
    break;

  case DYNL_NOW: 
    md_mode = RTLD_NOW;
   break;

  default:
    md_mode = 0;
    break;
  }

  if(mode & DYNL_GLOBAL)
    md_mode |= RTLD_GLOBAL;

  handle = dlopen(fname, md_mode);

  if (handle == NULL) fprintf(stderr, "dynl_link: %s\n", dlerror());

  return handle;
}

int dynl_unlink(void *handle) {
  int status = dlclose(handle);

  if (status != 0) fprintf(stderr, "dynl_unlink: %s\n", dlerror());

  return status;
}

/* --- CAUTION ---
 * NULL/DYNL_SELF/DYNL_NEXT is not correctly worked from shared/loaded object
 */

/* Expand dlsym wrapper function from template file */
#define FUNC		dynl_symbol
#define FUNCNAME	"dynl_symbol"
#define FUNC_T_		void *
#define FUNC_T		void *
#define DLSYM_T		void *
#define DLSYM		dlsym
#include <sim/dynl-dlsym.c.in>

#if defined(__FreeBSD__)
/* Use dlfunc() interface to get function pointer from (handle, symbol) */

#undef	FUNC
#undef	FUNCNAME
#undef	FUNC_T_
#undef	FUNC_T
#undef	DLSYM_T
#undef	DLSYM
#define FUNC		dynl_function
#define FUNCNAME	"dynl_function"
#define FUNC_T_		dynl_func_t
#define FUNC_T		dynl_func_t
#define DLSYM_T		dlfunc_t
#define DLSYM		dlfunc
#include <sim/dynl-dlsym.c.in>

#else
/* Use pointer-cast hack to function pointer from dlsym() */

#undef	FUNC
#undef	FUNCNAME
#undef	FUNC_T_
#define FUNC		__dynl_function
#define FUNCNAME	"dynl_function"
#define FUNC_T_		static void *
#include <sim/dynl-dlsym.c.in>
dynl_func_t dynl_function(void *handle, const char *symbol) {
  dynl_func_t func;

  if(sizeof(dynl_func_t) != sizeof(void*)) {
    fprintf(stderr,
	    "The function pointer representation disagrees with"
	    " the object pointer!\n");
    exit(EX_CONFIG);
  }

  /*
   * This code is a workaround for cast problem between object pointer
   * and function pointer. This workaround assumes that
   * the representation of the function pointer is equivalent with
   * the object pointer.
   * But, this assumption is not guaranteed by ISO C standard.
   */
  *(void **)&func = __dynl_function(handle, symbol);

  return func;
}
#endif

/* End of File */
