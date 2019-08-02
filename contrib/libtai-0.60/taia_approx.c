#include "taia.h"

double taia_approx(t)
struct taia *t;
{
  return tai_approx(&t->sec) + taia_frac(t);
}
