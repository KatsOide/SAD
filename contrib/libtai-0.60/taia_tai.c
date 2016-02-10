#include "taia.h"

void taia_tai(ta,t)
struct taia *ta;
struct tai *t;
{
  *t = ta->sec;
}
