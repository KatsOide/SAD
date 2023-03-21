/* MD part for dynamic object loader(dummy) */

#include <sim/dynl.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

static bool dynl_initialized = false;

int dynl_error_report_level = 0;

bool dynl_provide(void) {
  return false;
}

bool dynl_init(void) {
  if(!dynl_initialized) {
    dynl_initialized = true;
    return true;
  }

  return true;
}

void *dynl_link(const char *fname, int mode) {
  if(dynl_error_report_level > 0)
    fprintf(stderr, "dynl_link: Not supported\n");

  return NULL;
}

int dynl_unlink(void *handle) {
  if(dynl_error_report_level > 0)
    fprintf(stderr, "dynl_unlink: Not supported\n");

  return -1;
}

void *dynl_symbol(void *handle, const char *symbol) {
  if(dynl_error_report_level > 0)
    fprintf(stderr, "dynl_symbol: Not supported\n");

  return NULL;
}

dynl_func_t dynl_function(void *handle, const char *symbol) {
  if(dynl_error_report_level > 0)
    fprintf(stderr, "dynl_symbol: Not supported\n");

  return NULL;
}

/* End of File */
