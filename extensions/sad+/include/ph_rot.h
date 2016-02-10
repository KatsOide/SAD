#ifndef PH_ROT_H
#define PH_ROT_H
#include <element.h>

class Ph_rot : public Element
{
   double nux,nuy,nuz;
   double mux,muy,muz;
   double twiss[18];
   matrix* Linmap;
protected:
   
public:
   Ph_rot(double *K);
   Ph_rot(char*);
   ~Ph_rot() {delete Linmap;}
   void print(void) { Element::print();
      cout << *Linmap << "\n"; }

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);

};


#endif
