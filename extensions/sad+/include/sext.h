#ifndef SEXT_H
#define SEXT_H
#include <element.h>
 
class Sext : public Element
{
   double K_2;

public:

   Sext(char* s) : Element(s)
   {
      K_2=get_parm(s,"K2");
   }
   Sext(double* K) : Element(K)
   {
      K_2=K[1];
   }
   void print(void) { Element::print(); 
      cout << ", K2=" << K_2 << ")\n";}

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   
};

class Octu : public Element
{
   double K_3;

public:

   Octu(char* s) : Element(s)
   {
      K_3=get_parm(s,"K3");
   }
   Octu(double* K) : Element(K)
   {
      K_3=K[1];
   }
   void print(void) { Element::print(); 
      cout << ", K3=" << K_3 << ")\n";}

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   
};

#endif
