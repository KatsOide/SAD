#ifndef _BUILDINFO_H_
#define _BUILDINFO_H_

#include <sys/types.h>
#include <stdbool.h>
#include <stdint.h>

#ifndef _BUILDINFO_C_
typedef struct {
  int dummy_argument;
} buildinfo_t;
#else /* _BUILDINFO_C_ */
typedef enum {
  _NOP,
  _INT,
  _STR,
} buildinfo_type_t;

typedef struct {
  const char* keyword;
  buildinfo_type_t type;
  int32_t integer;
  const char* string;
} buildinfo_t;
#endif /* _BUILDINFO_C_ */

/* Database structure */
extern const buildinfo_t buildinfo_db[];

/* Callback proto-type for Scan API */
typedef void (*buildinfo_scan_callback_t)(const buildinfo_t*);

/* Scan API */
size_t buildinfo_scan(buildinfo_scan_callback_t callback);

/* Search API */
const buildinfo_t* buildinfo_search(const char* keyword);

/* Type test API */
bool buildinfo_string_entryQ(const buildinfo_t* entry);
bool buildinfo_integer_entryQ(const buildinfo_t* entry);

/* Value extract API */
const char* buildinfo_extract_keyword(const buildinfo_t* entry);
const char* buildinfo_extract_string(const buildinfo_t* entry);
int buildinfo_extract_integer(const buildinfo_t* entry);

/* Initialize test */
bool buildinfo_init(void);

#endif /* _BUILDINFO_H_ */
