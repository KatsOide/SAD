#ifndef BEND_H
#define BEND_H
#include <element.h>

class Bend : public Element
{
   double phi_0,e1,e2,phi,dphi,omega,rho;
   double dphix,dphiy;
   double sin_e1,sin_e2,cos_e1,cos_e2,tan_e1,tan_e2;
   double sin_omega,cos_omega,sin_p2,cos_p2,sin_pp;

public:
   Bend(double* K);
   Bend(char*);
   void print(void) { Element::print(); 
      cout<< ", phi_0=" << phi_0 << " ,dphi=" << dphi << 
	", e1=" << e1 << ", e2=" << e2 << ")\n";}

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

#endif
