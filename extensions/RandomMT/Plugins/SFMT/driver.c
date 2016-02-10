/*
 * Mersenne Twister(SFMT) Plugin
 */

#include <random_driver.h>
#include <feature.h>

#include <math.h>

/* Macro for mtdriver */
#define DRIVER_VER_MAJOR	1UL
#define DRIVER_VER_MINOR	0UL
#define PLUGIN_HEADER_MAGIC	0x0000deadUL

/* Internal State Macro */
#define N_MTI		N32
#define mt_index	idx
#define mt_state	psfmt32

/* MT API Macro */
typedef uint32_t	init_array_t;
#define MT_SEED_DEFAULT	1234UL
#define	init_by_seed	init_gen_rand
#define	init_by_array	init_by_array
#define	_genrand_uint32	gen_rand32
#ifndef USE_INLINE_GENRAND
#define	genrand_uint32	gen_rand32
#else
#define	genrand_uint32()	((mt_index >= N_MTI) ? (gen_rand_all(), mt_index = 1, mt_state[0]) : mt_state[mt_index++])
#endif

/* Include Mersenne Twister reference code */
#include "SFMT.h"
#include "SFMT.c"

/* Include Mersenne Twister Plugin Common Part */
#include "mtdriver.c"

#define _cat2(a, b)		a##b
#define _cat3(a, b, c)		a##b##c

/* Plugin header */
#define _plugin_SFMT(mexp, bit)	_cat3(_plugin_SFMT, mexp, _##bit)
#define init_framework_RandomPlugin_SFMT(mexp)	_cat2(init_framework_RandomPlugin_SFMT, mexp)

static const random_plugin_t _plugin_SFMT(MEXP, 53) = {
  RANDOM_ABI_VERSION,
  FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
  "SFMT"STR_MEXP"/53bit", IDSTR,
  INTERNAL_STATE_SIZE,
  dump_internal_state,
  restore_internal_state,
  is_supported,
  generate53,
  _genrand_uint32,
  NULL
};

static const random_plugin_t _plugin_SFMT(MEXP, 32) = {
  RANDOM_ABI_VERSION,
  FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
  "SFMT"STR_MEXP"/32bit", IDSTR,
  INTERNAL_STATE_SIZE,
  dump_internal_state,
  restore_internal_state,
  is_supported,
  generate32,
  _genrand_uint32,
  NULL
};

/* Plugin auto-registration at Library@Require[] */
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
  int init_framework_RandomPlugin_SFMT(MEXP)(void) {
#endif /* WITH_EXTENSION_MODULE */
  uint32_t init_vector[] = {0x1234, 0x5678, 0x9abc, 0xdef0};

  if(!feature_require("Random/Framework", FEATURE_VERSION(1, 6)))
    return -1;

  if(!feature_provide("Random/SFMT"STR_MEXP,
		      FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
		      FEATURE_VERSION(DRIVER_VER_MAJOR, 0)))
    return -1;

  /* Initialize internal state */
  restore_internal_state(init_vector,
			 sizeof(init_vector) / sizeof(uint32_t));

  /* Register SFMT plugins */
  random_register(&_plugin_SFMT(MEXP, 53));
  random_register(&_plugin_SFMT(MEXP, 32));

  return 0;
}

/* End of File */
