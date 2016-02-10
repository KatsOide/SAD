#ifndef _DYNL_H_
#define _DYNL_H_

#include <stdbool.h>

extern int dynl_error_report_level;	/* Error report level */

/* MD-API link flags */
#define DYNL_LAZY	1
#define DYNL_NOW	2
#define DYNL_MODEMASK	0x3
#define DYNL_GLOBAL	0x100
#define DYNL_LOCAL	0

/* MD-API special handle symbol */
#define DYNL_NEXT	((void *) -1)
#define DYNL_DEFAULT	((void *) -2)
#define DYNL_SELF	((void *) -3)
#define DYNL_NOHANDLE	((void *) -9)

/* Define function type for dynl_function() */
struct __dynl_func_arg {
  int __dynl_func_dummy;
};

typedef void (*dynl_func_t)(struct __dynl_func_arg);

#define DYNL_FUNC_NULL	((dynl_func_t)0)

/* MD-API function prototype */
/* dynl_provide() - Query dynamic loader backend support
 * Return value:
 *  true:	dynl MD routines is supported
 *  false:	dynl MD routines in not supported
 */
bool  dynl_provide(void);

/* dynl_init() - Initialize dynamic loader backend
 * Return value:
 *  true:	success
 *  false:	failed to initialize dynamic loader backend
 */
bool  dynl_init(void);

/* dynl_link(fname, mode) - Load shared object ``fname'' with mode ``mode''
 * Retuen value:
 *  NULL(failed to load object)
 *  otherwise(access handle pointer)
 */
void *dynl_link(const char *, int);

/* dynl_unlink(handle) - Unload shared object pointed by ``handle''
 * Retuen value:
 *  0(success)
 *  -1(failed to unload object)
 */
int   dynl_unlink(void *);

/* dynl_symbol(handle, symbol) - Resolve ``symbol'' address in ``handle'' object
 * Special handles to control symbol resolving scope:
 *  NULL(the executable or shared object from which the call is being made)
 *  DYNL_SELF(the shared object issuing dynl_symbol())
 *  DYNL_NEXT(the shared objects loaded after the one issuing dynl_symbol())
 *  DYNL_DEFAULT(search by default algorithm)
 * Return value:
 *  NULL(failed to resolve symbol address)
 *  otherwise(symbol address)
 */
void *dynl_symbol(void *, const char *);

/* dynl_function(handle, symbol) - dynl_symbol() for function pointer */
dynl_func_t dynl_function(void *, const char *);

/* dl_* family MI service functions */
int dl_link(const char*, int);
int dl_unlink(int);
int   dl_reg_handle(void*);
void* dl_ref_handle(int);

#endif /* _DYNL_H_ */
