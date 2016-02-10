#ifndef RANDOM_DRIVER_H
#define RANDOM_DRIVER_H

#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>

/* Random Plugin Framework ABI */
#define RANDOM_ABI_VERSION	0x00000105UL	/* Version 1.5 */
#define RANDOM_ABI_MASK		0x0000ffffUL
#define RANDOM_ABI_MAJOR(ver)	((ver & 0x0000ff00UL) >> 8)
#define RANDOM_ABI_MINOR(ver)	((ver & 0x000000ffUL))

/* Mode macro for random number generator */
#define RANDOM_PARABOLA		(-3)	/* Parabolic distribution: (-1,+1) */
#define RANDOM_GAUSS		(-2)	/* Gauss distribution: Sigma=1 */
#define RANDOM_GAUSS_NOCUT	(-1)	/* Gauss distribution: Sigma=1 */
#define RANDOM_UNIFORM		(0)	/* Uniform distribution: (0,1) */
#define RANDOM_UNIFORM0		(1)	/* Uniform distribution: [0,1) */
#define RANDOM_UNIFORM1		(2)	/* Uniform distribution: (0,1] */
#define RANDOM_UNIFORM01	(3)	/* Uniform distribution: [0,1] */

#define RANDOM_GCUT_MIN		1e-8

typedef struct {
  uint32_t abi;			/* plugin ABI version */
  uint32_t version;		/* plugin version */
  char *name;			/* plugin name string */
  char *sequence;		/* plugin sequence id string[optional] */
  size_t state_size;		/* internal state size */
  size_t  (*dump)(uint32_t *buffer, const size_t size);
  int  (*restore)(const uint32_t *buffer, const size_t size);
  int (*is_supported)(int mode);
  int (*generate)(int mode, size_t size, double *array);
  uint32_t (*generate32)(void);
  uint64_t (*generate64)(void);
} random_plugin_t;

extern double random_gauss_cut;

extern double	   random_get_gcut(void);
extern double	   random_set_gcut(double);

extern int	   random_register(const random_plugin_t *plugin);
extern int	   random_select(const char *name, const char *sequence);
extern const char *random_name(void);
extern const char *random_sequence(void);

extern void   random_scan_plugin_list(void (*func)(random_plugin_t *plugin, void *param), void *param);
extern size_t random_size(void);
extern size_t random_dump(uint32_t *buffer, size_t length);
extern int    random_restore(const uint32_t *buffer, size_t length);

extern int    random_check_mode(int mode);
extern int    random_generate(int mode, size_t length, double *array);
extern uint32_t random_generate32(void);
extern uint32_t random_generate64(void);

#endif
/* End of File */
