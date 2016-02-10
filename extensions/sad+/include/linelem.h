#ifndef LINELEM_H
#define LINELEM_H
#include <element.h>

#define nbslice 4

class LBend : public Element
{
   double phi_0,e1,e2,phi,dphi,omega,rho;
   double dphix,dphiy;
   double sin_e1,sin_e2,cos_e1,cos_e2,tan_e1,tan_e2;
   double a11,a12,a33,a34,a16,a51,a52,a56;
   double a11h,a12h,a33h,a34h,a16h,a51h,a52h,a56h;
   double G,ds0,K1,fx,fy;

public:
   LBend(double* K);
   LBend(char*);
   void print(void) { Element::print(); 
      cout<< ", phi_0=" << phi_0 << " ,dphi=" << dphi << 
	", e1=" << e1 << ", e2=" << e2 << ")\n";}

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

#endif
