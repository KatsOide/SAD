#ifndef QUAD_H
#define QUAD_H
#include <element.h>

class Quad : public Element
{
   int N_div;
   double K1,ds0,rk,sqrtk,f1,f2;
   double a11,a12,b11,b12,a11h,a12h,b11h,b12h;
   int ent_edge,exit_edge;

public:

   Quad(double* K);
   Quad(char*);

   void print(void) ;
   void set_ent_edge(void) { ent_edge=1;}
   void set_exit_edge(void) { exit_edge=1;}
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

#endif
