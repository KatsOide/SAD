#include <element.h>

class Wiggler : public Element
{
   int N_pole,N_div;
   double Bx,By,dphase;

   double Kz,Kx,Ky,Qx,Qy,lambda,ds;
   double F_0,G_0;

public:
   Wiggler(void)
	  : Element() {;}

   Wiggler(char*);
   Wiggler(double*);
   void print(void)
   { Element::print(); 
      cout << ", N_pole=" << N_pole << " , N_div=" << N_div
	 << ", By=" << By << ", Bx=" << Bx << ", dphase=" << dphase << ")\n";
      cout << "Kz=" << Kz << " Kx=" << Kx << ", Ky=" << Ky << '\n';
   }

   void beam_set(Beam& BEAM) {
      F_0=By/(C*BEAM.Energy()*1.e-9); G_0=Bx/(C*BEAM.Energy()*1.e-9);}
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   
};
