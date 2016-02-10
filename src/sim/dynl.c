#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sim/dynl.h>

#include <stdbool.h>
#include <stdio.h>

/* Handler table for dynamic loader */
#ifndef MAX_DL_HANDLES
#define MAX_DL_HANDLES	128
#endif
static void *dl_handle_tbl[MAX_DL_HANDLES];
static bool  dl_handle_tbl_initialized = false;

/* Dynamic loader MI helper function */
static void dl_handle_tbl_init(void) {
  if(!dl_handle_tbl_initialized) {
    for(int i = 0; i < MAX_DL_HANDLES; i++)
      dl_handle_tbl[i] = NULL;
    dl_handle_tbl_initialized = true;
  }
}

/* Dynamic loader MI API */
int dl_reg_handle(void *handle) {
  int i;

  for(i = 0; i < MAX_DL_HANDLES; i++)
    if(dl_handle_tbl[i] == NULL) {
      dl_handle_tbl[i] = handle;
      return i;
    } else if(dl_handle_tbl[i] == handle) return -1;

  return -2;
}

void* dl_ref_handle(int id) {
  switch(id) {
  case 0:
    return NULL;

  case -1:
    return DYNL_NEXT;

  case -2:
    return DYNL_DEFAULT;

  case -3:
    return DYNL_SELF;

  default:
    if(id >= 1 && id <= MAX_DL_HANDLES)
      return dl_handle_tbl[id - 1];
    break;
  }

  return DYNL_NOHANDLE;
}

int dl_link(const char *fname, int mode) {
  int id;
  void *handle;

  handle = dynl_link(fname, mode);
  if(handle == NULL) return -1;

  id = dl_reg_handle(handle);
  if(id >= 0) return id + 1;

  switch(id) {
  case -1:
    fprintf(stderr, "dl_link: %s is already loaded\n", fname);
    break;

  case -2:
    fprintf(stderr, "dl_link: Object handle table exceeded\n");
    dynl_unlink(handle);
    break;

  default:
    fprintf(stderr, "dl_link: Catch unkown ID code[%d]\n", id);
    dynl_unlink(handle);
    break;
  }
  return -1;
}

int dl_unlink(int id) {
  void *handle;
  int status;

  if(id >= 1 && id <= MAX_DL_HANDLES) {
    handle = dl_handle_tbl[id - 1];
    if(handle != NULL) {
      status = dynl_unlink(handle);
      if(status != 0) return -1;
      dl_handle_tbl[id - 1] = NULL;
      return 0;
    }
  }
  return -2;
}

/* Initialization hook */
int init_framework_DynlMI(void) {
  dl_handle_tbl_init();		// Initialize handle table
  dynl_init();			// Initialize MD part

  return 0;
}

/* End of File */
