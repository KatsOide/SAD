#ifndef SOLENOID_H
#define SOLENOID_H
#include <element.h>

class Solenoid : public Element
{
   double Ks,rk;

public:
   Solenoid(char*);
   Solenoid(double*);
   void print(void);

   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

#endif
