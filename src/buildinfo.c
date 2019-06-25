#define _BUILDINFO_C_
#include <buildinfo.h>
#include <feature.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sysexits.h>

/* Internal API */
static bool buildinfo_check_db(void) {
  bool success = true;

  for(const buildinfo_t* p = buildinfo_db; p->keyword != NULL; p++) {
    if(strlen(p->keyword) < 1) {
      success = false;
      fprintf(stderr, "build-info database contains "
	      "null keyword entry!\n");
    }
    switch(p->type) {
    case _INT:
      break;

    case _STR:
      if(p->string == NULL) {
	success = false;
	fprintf(stderr,	"build-info database contains "
		"invalid string type entry[%s]!\n", p->keyword);
      }
      break;

    default:
      success = false;
      fprintf(stderr, "build-info database contains "
	      "unknown type(%d) entry[%s]!\n", p->type, p->keyword);
      break;
    }
  }

  return success;
}

/* Public API */
size_t buildinfo_scan(buildinfo_scan_callback_t callback) {
  size_t n = 0;
  for(const buildinfo_t* p = buildinfo_db; p->keyword != NULL; p++) {
    callback(p); n += 1;
  }
  return n;
}

const buildinfo_t* buildinfo_search(const char* keyword) {
  for(const buildinfo_t* p = buildinfo_db; p->keyword != NULL; p++) {
    if(strcmp(p->keyword, keyword) == 0) {
      return p;
    }
  }
  return NULL;
}

bool buildinfo_string_entryQ(const buildinfo_t* entry) {
  return entry && entry->type == _STR;
}

bool buildinfo_integer_entryQ(const buildinfo_t* entry) {
  return entry && entry->type == _INT;
}

const char* buildinfo_extract_keyword(const buildinfo_t* entry) {
  return entry ? entry->keyword : NULL;
}

const char* buildinfo_extract_string(const buildinfo_t* entry) {
  return buildinfo_string_entryQ(entry) ? entry->string : NULL;
}

int buildinfo_extract_integer(const buildinfo_t* entry) {
  return buildinfo_integer_entryQ(entry) ? entry->integer : 0;
}

bool buildinfo_init(void) {
  static bool initialized = false;

  if(!initialized) {
    if(buildinfo_check_db()) {
      initialized = true;
    }
  }
  return initialized;
}

int init_framework_BuildInfo(void) {
  if(!buildinfo_init()) { exit(EX_CONFIG); }

  feature_provide("BuildInfo/Framework",
		  FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0));

  return 0;
}
/* End of File */
