/*
 * SAD Default Random Number Generator Plugin
 */

#include <random_driver.h>
#include <feature.h>

#include <stdlib.h>
#include <math.h>

#define DRIVER_VER_MAJOR	1UL
#define DRIVER_VER_MINOR	1UL

#define True	(0 == 0)
#define False	(0 != 0)

#define A	48828125
#define B	0

/* Size of internal state in uint32_t unit */
#define INTERNAL_STATE_SIZE	4

/* Internal state for random number generator */
struct internal_state_t {
  int32_t seed;
  int32_t filled;
  double next;
};

static union {
  struct internal_state_t state;
  uint32_t __ui[INTERNAL_STATE_SIZE];
} internal_state = {{0, False, 0.0}};

/* Random number generator state init/dump/restore */
static size_t dump_internal_state(uint32_t *buffer, size_t length) {
  int i;

  if(length >= 4) {
    for(i = 0; i < 4; i++)
      buffer[i] = internal_state.__ui[i];
    return 4;
  }

  return 0;
}

static int restore_internal_state(const uint32_t *buffer, size_t length) {
  int i;

  if(length >= 4) {
    for(i = 0; i < 4; i++)
      internal_state.__ui[i] = buffer[i];
    if(internal_state.state.seed < 0)
      internal_state.state.seed = -(abs(internal_state.state.seed / 2) * 2 + 1);
    else
      internal_state.state.seed = abs(internal_state.state.seed / 2) * 2 + 1;
    return 0;
  }

  if (length > 0) {
    internal_state.__ui[0] = buffer[0];
    if(internal_state.state.seed < 0)
      internal_state.state.seed = -(abs(internal_state.state.seed / 2) * 2 + 1);
    else
      internal_state.state.seed = abs(internal_state.state.seed / 2) * 2 + 1;
  } else {
    internal_state.state.seed = 17;
  }

  internal_state.state.filled = False;
  internal_state.state.next   = 0.0;

  return 0;
}

/* Random number generator */
static int is_supported(int mode) {
  switch(mode) {
  case RANDOM_PARABOLA:		/* Parabola */
  case RANDOM_GAUSS:		/* Gaussian */
  case RANDOM_GAUSS_NOCUT:	/* Gaussian */
  case RANDOM_UNIFORM:		/* (0, 1) */
  case RANDOM_UNIFORM0:		/* [0, 1) */
  case RANDOM_UNIFORM1:		/* (0, 1] */
  case RANDOM_UNIFORM01:	/* [0, 1] */
    return 0;
    break;

  default:
    return -1;
    break;
  }
  return -1;
}

static int generate32(int mode, size_t length, double *array) {
  struct internal_state_t *state = &internal_state.state;

  int i;
  double x, y;

  switch(mode) {
  case RANDOM_PARABOLA:		/* Parabola */
    for(i = 0; i < length; i++) {
#ifdef USE_PARABOLA_TRANSFORM
      state->seed = A * state->seed + B;
      x = ((double)state->seed + 0.5)*(1.0/2147483648.0);	/* (-1, 1) */
      x = 2.0 * sin(asin(x) / 3.0);
#else
      do {
	state->seed = A * state->seed + B;
	x = state->seed*(1.0/2147483648.0);		/* [-1, 1) */
	state->seed = A * state->seed + B;
	y = state->seed*(1.0/4294967296.0) + 0.5;	/* [ 0, 1) */
      } while(x * x > y);
#endif
      array[i] = x;
    }
    break;

  case RANDOM_GAUSS:		/* Gaussian */
    for(i = 0; i < length; i++) {
      if(!state->filled) {
	state->seed = A * state->seed + B;
	x = state->seed*(M_PI/2147483648.0); /* phi [-Pi, Pi) */
	state->seed = A * state->seed + B;
	y = sqrt(-2.0 * log(state->seed*(1.0/4294967296.0) + 0.5)); /* rho */

	array[i] = y * cos(x); i += 1;
	if(fabs(array[i-1]) > random_gauss_cut) i -= 1;
	if(i < length) {
	  array[i] = y * sin(x);
	  if(fabs(array[i]) > random_gauss_cut) i -= 1;
	} else {
	  state->next = y * sin(x); state->filled = True;
	}
      } else {
	array[i] = state->next; state->filled = False; state->next = 0.0;
	if(fabs(array[i]) > random_gauss_cut) i -= 1;
      }
    }
    break;

  case RANDOM_GAUSS_NOCUT:	/* Gaussian */
    for(i = 0; i < length; i++) {
      if(!state->filled) {
	state->seed = A * state->seed + B;
	x = state->seed*(M_PI/2147483648.0); /* phi [-Pi, Pi) */
	state->seed = A * state->seed + B;
	y = sqrt(-2.0 * log(state->seed*(1.0/4294967296.0) + 0.5)); /* rho */

	array[i] = y * cos(x); i += 1;
	if(i < length) {
	  array[i] = y * sin(x);
	} else {
	  state->next = y * sin(x); state->filled = True;
	}
      } else {
	array[i] = state->next; state->filled = False; state->next = 0.0;
      }
    }
    break;

  case RANDOM_UNIFORM:		/* (0, 1) */
    for(i = 0; i < length; i++) {
      state->seed = A * state->seed + B;
      array[i] = ((double)state->seed + 2147483648.5)*(1.0/4294967296.0);
    }
    break;

  case RANDOM_UNIFORM0:		/* [0, 1) */
    for(i = 0; i < length; i++) {
      state->seed = A * state->seed + B;
      array[i] = ((double)state->seed + 2147483648.0)*(1.0/4294967296.0);
    }
    break;

  case RANDOM_UNIFORM1:		/* (0, 1] */
    for(i = 0; i < length; i++) {
      state->seed = A * state->seed + B;
      array[i] = ((double)state->seed + 2147483649.0)*(1.0/4294967296.0);
    }
    break;

  case RANDOM_UNIFORM01:	/* [0, 1] */
    for(i = 0; i < length; i++) {
      state->seed = A * state->seed + B;
      array[i] = ((double)state->seed + 2147483648.0)*(1.0/4294967295.0);
    }
    break;

  default:
    return -1;
    break;
  }

  return 0;
}

static uint32_t genrand_uint32(void) {
  struct internal_state_t *state = &internal_state.state;

  state->seed = A * state->seed + B;
  return state->seed;
}

/* Plugin header */
static const random_plugin_t _plugin_SAD = {
  RANDOM_ABI_VERSION,
  FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
  "SAD", NULL,
  INTERNAL_STATE_SIZE,
  dump_internal_state,
  restore_internal_state,
  is_supported,
  generate32,
  genrand_uint32,
  NULL
};

/* Plugin auto-registration at Library@Require[] */
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int init_framework_RandomPlugin_SAD(void) {
#endif /* WITH_EXTENSION_MODULE */
  uint32_t init_vector[] = {17};

  if(!feature_require("Random/Framework", FEATURE_VERSION(1, 6)))
    return -1;

  if(!feature_provide("Random/SAD",
		      FEATURE_VERSION(DRIVER_VER_MAJOR, DRIVER_VER_MINOR),
		      FEATURE_VERSION(DRIVER_VER_MAJOR, 0)))
    return -1;

  /* Initialize internal state */
  restore_internal_state(init_vector,
			 sizeof(init_vector) / sizeof(uint32_t));

  /* Register plugins */
  random_register(&_plugin_SAD);

  return 0;
}

/* End of File */
