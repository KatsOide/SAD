#ifndef _FEATURE_H_
#define _FEATURE_H_

#include <stdbool.h>
#include <stdint.h>

#define FEATURE_VERSION(major, minor)	((major) * 1000UL + (minor))
#define FEATURE_MAJOR(version)		((version) / 1000UL)
#define FEATURE_MINOR(version)		((version) % 1000UL)

extern bool feature_provide(const char*, uint32_t, uint32_t);
extern bool feature_require(const char*, uint32_t);
extern uint32_t feature_version(const char*);
extern int feature_scan(void (*)(const char*));

#endif /* _FEATURE_H_ */
