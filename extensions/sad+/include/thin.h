#ifndef THIN_H
#define THIN_H
#include <element.h>
 
class Thin : public Element
{
   int n_MP;
   double K_n;

public:
   Thin(double* K)
	  : Element(K) 
   {
     n_MP=(int)K[1];
     K_n=K[2];
     //cout <<" Thin element\n";  
   }
   Thin(char*s) : Element(s)
   {
      n_MP=(int)get_parm(s,"n_MP");
      K_n=get_parm(s,"K_n");
   }
   void print(void) { Element::print(); 
      cout << ", n_MP=" << n_MP << ", K_n=" << K_n << ")\n";}

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   
};

#endif
