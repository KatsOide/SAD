#ifndef IMPEDANCE_H
#define IMPEDANCE_H
#include <element.h>
#include <phys_const.h>

class Impedance : public Element
{
   int niz;
   double Zw,dZ,kfac;
   double *Z,*rho0,*rho1x,*rho1y,*W0P,*W1x,*W1y,*F0,*F1x,*F1y,*G0,*G1x,*G1y;
   double K0,K1x,K1y; 
public:

   Impedance(double* K);
   Impedance(char*);
   ~Impedance(void);

   void beam_set(Beam& BEAM) {
      kfac=e*BEAM.N_particle()/(BEAM.Energy()*1.e9); 
      cout << "kfac is set to " << kfac << '\n';
   }
   void print(void) ;
   void printrho(int);
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};


class Chromaticity : public Element
{
   double Kxp,Kyp;
public:

   Chromaticity(double* K);
   Chromaticity(char*);
   //   ~Chromaticity(void);

   void print(void) ;
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

#endif
