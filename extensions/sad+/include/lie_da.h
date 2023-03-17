#ifndef LIE_DA_H
#define LIE_DA_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_da.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <dacpp.h>
#include <c_da.h>
#include <matrix.h>
#include <lin_map.h>
#include <map_double.h>
#include <map_da.h>
map_da lie_exp(const da&,const map_da&);
map_da lie_exp(const da&);
da lie_exp(const da&,const da&);
map_da fac_map_type1(const da&);
map_da fac_map_type2(const da&);
map_da fac_map_type3(const da&);
da fac_drg_type1(const map_da&);
da fac_drg_type2(const map_da&);
da fac_drg_type3(const map_da&);
da fac_drg_type1(const map_da&,map_da&);
da fac_drg_type2(const map_da&,map_da&);
da fac_drg_type3(const map_da&,map_da&);
da can_perturbation(const map_da&);
da can_perturbation(const map_da&,da&);
c_da Normal_expression(const da&);
c_da Real_expression(const c_da&);
void nonlinear_dist(const map_da&);

#endif
