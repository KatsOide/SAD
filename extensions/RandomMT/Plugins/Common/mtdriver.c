/*
 * Mersenne Twister Plugin Common Part
 */

#include <stdlib.h>

#define True	(0 == 0)
#define False	(0 != 0)

/* Size of internal state in uint32_t unit
 * Header Magic:		1
 * Mersenne Twister:		N_MTI + 2
 * Box-Muller Gauss Generator:	3
 */
#define MERSENNE_TWISTER_STATE_SIZE	(N_MTI + 2)
#define BOX_MULLER_STATE_SIZE		(3)
#define INTERNAL_STATE_SIZE	(1 \
				 + MERSENNE_TWISTER_STATE_SIZE \
				 + BOX_MULLER_STATE_SIZE \
				 + 0)

#define STATE_HEADER_MAGIC	((  DRIVER_VER_MAJOR << 24) \
				 | (DRIVER_VER_MINOR << 16) \
				 | (PLUGIN_HEADER_MAGIC & 0x0000ffffUL))

/* Internal state for Box-Muller generator */
struct box_muller_gauss_t {
  double next;
  int32_t filled;
};

static union {
  struct box_muller_gauss_t state;
  uint32_t __ui[BOX_MULLER_STATE_SIZE];
} box_muller_gauss = {{0.0, False}};

/* Mersenne Twister generator state dump/restore */
static uint32_t *dump_mersenne_twister(uint32_t *buffer) {
  int i;

  /* dump state[mt_index, N_MTI, mt_state[0,1,...,N_MTI-1]] */
  *buffer = mt_index;	buffer++; /* dump state register */
  *buffer = N_MTI;	buffer++; /* dump state vector length */
  for(i = 0; i < N_MTI; i++, buffer++)
    *buffer = mt_state[i];

  return buffer;
}

static const uint32_t *restore_mersenne_twister(const uint32_t *buffer) {
  int i;

  if(buffer[1] != N_MTI) return NULL;

  mt_index = *buffer; /* restore state register */
  buffer += 2;
  for(i = 0; i < N_MTI; i++, buffer++)
    mt_state[i] = *buffer;

  return buffer;
}

/* Box-Mullder gauss generator state init/dump/restore */
static void init_box_muller_gauss(void) {
  struct box_muller_gauss_t *state = &box_muller_gauss.state;

  state->filled = False; state->next = 0.0;
}

static uint32_t *dump_box_muller_gauss(uint32_t *buffer) {
  uint32_t *src = box_muller_gauss.__ui;
  int i;

  for(i = 0; i < BOX_MULLER_STATE_SIZE; i++, src++, buffer++)
    *buffer = *src;

  return buffer;
}

static const uint32_t *restore_box_muller_gauss(const uint32_t *buffer) {
  uint32_t *dst = box_muller_gauss.__ui;
  int i;

  for(i = 0; i < BOX_MULLER_STATE_SIZE; i++, dst++, buffer++)
    *dst = *buffer;

  return buffer;
}

/* Random number generator state init/dump/restore */
static int init_internal_state(const uint32_t *buffer, size_t length) {
  init_array_t *array;
  int i;

  if(length > 0)
    array = malloc(length * sizeof(init_array_t));

  if(length > 0 && array != NULL) {
    for(i = 0; i < length; i++)
      array[i] = buffer[i];
    init_by_array(array, length);
    free(array);
  } else
    init_by_seed(MT_SEED_DEFAULT);

  init_box_muller_gauss();

  return 0;
}

static size_t dump_internal_state(uint32_t *buffer, size_t length) {
  uint32_t *buffer0 = buffer;

  if(length - (buffer - buffer0) >= 1) {
    *buffer = STATE_HEADER_MAGIC; buffer++;
  } else return 0;

  /* Mersenne Twister */
  if(length - (buffer - buffer0) >= MERSENNE_TWISTER_STATE_SIZE) {
    buffer  = dump_mersenne_twister(buffer);
  } else return 0;

  /* Box-Muller Gauss */
  if(length - (buffer - buffer0) >= BOX_MULLER_STATE_SIZE) {
    buffer = dump_box_muller_gauss(buffer);
  } else return 0;

  return (buffer - buffer0);
}

static int restore_internal_state(const uint32_t *buffer, size_t length) {
  const uint32_t *buffer0 = buffer;

  if(length - (buffer - buffer0) >= 1) {
    if(*buffer != STATE_HEADER_MAGIC)
      return init_internal_state(buffer0, length);
    buffer += 1;
  } else return init_internal_state(buffer0, length);

  /* Mersenne Twister */
  if(length - (buffer - buffer0) >= MERSENNE_TWISTER_STATE_SIZE) {
    buffer  = restore_mersenne_twister(buffer);
    if(buffer == NULL)
      return init_internal_state(buffer0, length);
  } else return init_internal_state(buffer0, length);

  /* Box-Muller Gauss */
  if(length - (buffer - buffer0) >= BOX_MULLER_STATE_SIZE) {
    buffer  = restore_box_muller_gauss(buffer);
    if(buffer == NULL)
      init_box_muller_gauss();    
  } else init_box_muller_gauss();    

  return 0;
}

/* Random number generator state init/dump/restore */
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

static int generate53(int mode, size_t length, double *array) {
  struct box_muller_gauss_t *state = &box_muller_gauss.state;

  int i;
  unsigned long a, b, c;
  double x, y, r;

  switch(mode) {
  case RANDOM_PARABOLA:		/* Parabola */
    for(i = 0; i < length; i++) {
#ifdef USE_PARABOLA_TRANSFORM
      do {
	a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
      } while(a == 0UL && b == 0UL);
      x = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0;	/* (-1, 1) */
      x = 2.0 * sin(asin(x) / 3.0);
#else
      do {
	a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	x = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0;	/* [-1, 1) */
	a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	y = (a*67108864.0+b)*(1.0/9007199254740992.0);		/* [ 0, 1) */
      } while(x * x > y);
#endif
      array[i] = x;
    }
    break;

  case RANDOM_GAUSS:		/* Gaussian */
    for(i = 0; i < length; i++) {
      if(!state->filled) {
	do { /* [-1, 1) x [-1, 1) */
	  a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	  x = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0; 
	  a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	  y = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0;
	  r = x * x + y * y;
	} while(r >= 1 || r == 0);
	r = sqrt(-2.0 * log(r) / r);

	array[i] = x * r; i += 1;
	if(fabs(array[i-1]) > random_gauss_cut) i -= 1;
	if(i < length) {
	  array[i] = y * r;
	  if(fabs(array[i]) > random_gauss_cut) i -= 1;
	} else {
	  state->next = y * r; state->filled = True;
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
	do { /* [-1, 1) x [-1, 1) */
	  a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	  x = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0; 
	  a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
	  y = (a*67108864.0+b)*(1.0/4503599627370496.0) - 1.0;
	  r = x * x + y * y;
	} while(r >= 1 || r == 0);
	r = sqrt(-2.0 * log(r) / r);

	array[i] = x * r; i += 1;
	if(i < length) {
	  array[i] = y * r;
	} else {
	  state->next = y * r; state->filled = True;
	}
      } else {
	array[i] = state->next; state->filled = False; state->next = 0.0;
      }
    }
    break;

  case RANDOM_UNIFORM:		/* (0, 1) */
    for(i = 0; i < length; i++) {
      do {
	a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
      } while(a == 0UL && b == 0UL);
      array[i] = (a*67108864.0+b)*(1.0/9007199254740992.0);
    }
    break;

  case RANDOM_UNIFORM0:		/* [0, 1) */
    for(i = 0; i < length; i++) {
      a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
      array[i] = (a*67108864.0+b)*(1.0/9007199254740992.0);
    }
    break;

  case RANDOM_UNIFORM1:		/* (0, 1] */
    for(i = 0; i < length; i++) {
      a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
      if(a == 0UL && b == 0UL)
	array[i] = 1.0;
      else
	array[i] = (a*67108864.0+b)*(1.0/9007199254740992.0);
    }
    break;

  case RANDOM_UNIFORM01:	/* [0, 1] */
    for(i = 0; i < length; i++) {
#ifdef USE_TRUE_53BIT_RES
      do { /* Modular by 2^53+1 */
	a = genrand_uint32(); b = genrand_uint32();
      } while(a > 0xffe00000UL || (a == 0xffe00000UL && b >= 0x07ffUL));
      c = a >> 21;
      a = (a & 0x001fffffUL) << 6 | (b >> 26);
      b = b & 0x03ffffffUL;
      if(c > b) {
	if(a > 0UL) {
	  a -= 1UL;
	  b += 0x04000000UL;
	} else if(c == b + 1) {
	  a += 0x08000000UL;
	  b += 0x00000001UL;
	} else {
	  a += 0x07ffffffUL;
	  b += 0x04000001UL;
	}
      }
      b -= c;
      if(a > 0x07ffffffUL)
	array[i] = 1.0;
      else
	array[i] = (a*67108864.0+b)*(1.0/9007199254740992.0);
#else
      a = genrand_uint32() >> 5; b = genrand_uint32() >> 6;
      array[i] = (a*67108864.0+b)*(1.0/9007199254740991.0);
#endif
    }
    break;

  default:
    return -1;
    break;
  }

  return 0;
}

static int generate32(int mode, size_t length, double *array) {
  struct box_muller_gauss_t *state = &box_muller_gauss.state;

  int i;
  double x, y, r;

  switch(mode) {
  case RANDOM_PARABOLA:		/* Parabola */
    for(i = 0; i < length; i++) {
#ifdef USE_PARABOLA_TRANSFORM
      x = ((double)genrand_uint32() + 0.5)*(1.0/2147483648.0) - 1.0;	/* (-1, 1) */
      x = 2.0 * sin(asin(x) / 3.0);
#else
      do {
	x = genrand_uint32()*(1.0/2147483648.0) - 1.0;	/* [-1, 1) */
	y = genrand_uint32()*(1.0/4294967296.0);	/* [ 0, 1) */
      } while(x * x > y);
#endif
      array[i] = x;
    }
    break;

  case RANDOM_GAUSS:		/* Gaussian */
    for(i = 0; i < length; i++) {
      if(!state->filled) {
	do { /* [-1, 1) x [-1, 1) */
	  x = genrand_uint32()*(1.0/2147483648.0) - 1.0;
	  y = genrand_uint32()*(1.0/2147483648.0) - 1.0;
	  r = x * x + y * y;
	} while(r >= 1 || r == 0);
	r = sqrt(-2.0 * log(r) / r);

	array[i] = x * r; i += 1;
	if(fabs(array[i-1]) > random_gauss_cut) i -= 1;
	if(i < length) {
	  array[i] = y * r;
	  if(fabs(array[i]) > random_gauss_cut) i -= 1;
	} else {
	  state->next = y * r; state->filled = True;
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
	do { /* [-1, 1) x [-1, 1) */
	  x = genrand_uint32()*(1.0/2147483648.0) - 1.0;
	  y = genrand_uint32()*(1.0/2147483648.0) - 1.0;
	  r = x * x + y * y;
	} while(r >= 1 || r == 0);
	r = sqrt(-2.0 * log(r) / r);

	array[i] = x * r; i += 1;
	if(i < length) {
	  array[i] = y * r;
	} else {
	  state->next = y * r; state->filled = True;
	}
      } else {
	array[i] = state->next; state->filled = False; state->next = 0.0;
      }
    }
    break;

  case RANDOM_UNIFORM:		/* (0, 1) */
    for(i = 0; i < length; i++)
      array[i] = ((double)genrand_uint32() + 0.5)*(1.0/4294967296.0);
    break;

  case RANDOM_UNIFORM0:		/* [0, 1) */
    for(i = 0; i < length; i++)
      array[i] = ((double)genrand_uint32()      )*(1.0/4294967296.0);
    break;

  case RANDOM_UNIFORM1:		/* (0, 1] */
    for(i = 0; i < length; i++)
      array[i] = ((double)genrand_uint32() + 1.0)*(1.0/4294967296.0);
    break;

  case RANDOM_UNIFORM01:	/* [0, 1] */
    for(i = 0; i < length; i++)
      array[i] = ((double)genrand_uint32()      )*(1.0/4294967295.0);
    break;

  default:
    return -1;
    break;
  }

  return 0;
}

/* End of File */
