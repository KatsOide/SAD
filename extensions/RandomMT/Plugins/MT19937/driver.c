/*
 * Mersenne Twister(MT19937) Plugin
 */

#include <random_driver.h>
#include <feature.h>

#include <math.h>

/* Macro for mtdriver */
#define DRIVER_VER_MAJOR	1UL
#define DRIVER_VER_MINOR	0UL
#define PLUGIN_HEADER_MAGIC	0x0000beafUL

/* Internal State Macro */
#define N_MTI		N
#define mt_index	mti
#define mt_state	mt

/* MT API Macro */
typedef unsigned long	init_array_t;
#define MT_SEED_DEFAULT	5489UL
#define	init_by_seed	init_genrand
#define	init_by_array	init_by_array
#define	genrand_uint32	genrand_int32

/* Include Mersenne Twister reference code */
#define main core_main
#include "core.c"

/* Include Mersenne Twister Plugin Common Part */
#include "mtdriver.c"

/* proto-type matching for plugin API */
static uint32_t _genrand_uint32(void) {
  return genrand_uint32();
}

/* Plugin header */
static const random_plugin_t _plugin_MT19937_53 = {
  RANDOM_ABI_VERSION,
  FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
  "MT19937/53bit", NULL,
  INTERNAL_STATE_SIZE,
  dump_internal_state,
  restore_internal_state,
  is_supported,
  generate53,
  _genrand_uint32,
  NULL
};

static const random_plugin_t _plugin_MT19937_32 = {
  RANDOM_ABI_VERSION,
  FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
  "MT19937/32bit", NULL,
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
int init_framework_RandomPlugin_MT19937(void) {
#endif /* WITH_EXTENSION_MODULE */
  uint32_t init_vector[] = {291, 564, 837, 1110};

  if(!feature_require("Random/Framework", FEATURE_VERSION(1, 6)))
    return -1;

  if(!feature_provide("Random/MT19937",
		      FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
		      FEATURE_VERSION(DRIVER_VER_MAJOR, 0)))
    return -1;

  /* Initialize internal state */
  restore_internal_state(init_vector,
			 sizeof(init_vector) / sizeof(uint32_t));

  /* Register MT19937 plugins */
  random_register(&_plugin_MT19937_53);
  random_register(&_plugin_MT19937_32);

  return 0;
}

/* End of File */
