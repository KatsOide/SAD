#ifndef LEAPSECS_H
#define LEAPSECS_H

#include <libtai/tai.h>

extern int leapsecs_init(void);
extern int leapsecs_read(void);

extern void leapsecs_add(struct tai *, int);
extern int leapsecs_sub(struct tai *);

#endif
