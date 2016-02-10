#include "tai.h"

void tai_add(t,u,v)
struct tai *t;
struct tai *u;
struct tai *v;
{
  t->x = u->x + v->x;
}
