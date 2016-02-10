#include "caldate.h"

void caldate_normalize(cd)
struct caldate *cd;
{
  caldate_frommjd(cd,caldate_mjd(cd),(int *) 0,(int *) 0);
}
