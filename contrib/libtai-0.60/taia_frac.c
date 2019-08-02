#include "taia.h"

double taia_frac(t)
struct taia *t;
{
  return (t->atto * 0.000000001 + t->nano) * 0.000000001;
}
