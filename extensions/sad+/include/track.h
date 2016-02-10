#ifndef TRACK_H
#define TRACK_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_double.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <iostream>
using std::istream;
using std::ostream;

#include <dacpp.h>
#include <map_double.h>
#include <matrix.h>
//#include <element.h>

//#include <map_da.h>

class pBeam		//:private da
{
public:
   int N_particle;
   int np;
   int ndump;
   double rn_particle;
   double* x;
   double *y,*z,*px,*py,*pz;
   int* sv;
public:
   pBeam(char* str,int N=1000,double rn=0.);
   pBeam(int N=1000,double rn=0.);
   pBeam(const pBeam&);
   ~pBeam(void){ delete [] sv; delete [] x;}
   void Initialize(const char*, const matrix&);
//   void Initialize(const char*, EMIT&);
//   void Initialize(const char*, const Beam&);
//   void GaussianDistribution(const Beam&);
   void SetNdump(int n) {ndump=n;}

   friend ostream& operator<<(ostream&,pBeam&);
   friend istream& operator>>(istream&,pBeam&);
   pBeam& operator=(const pBeam&);
   pBeam& operator=(double);
   friend void is_survive(pBeam&);
   friend void BeamSizeMonitor(const pBeam&,double*);
   friend void BeamSizeMonitor(const pBeam&,double*,matrix&);

   map_double particle(int);

   double* operator[](int i) {
     return x+i*N_particle;
   }
// .m(0)=x
   double* m(int i) {
     return x+i*N_particle;
   }

   int N(void) { return N_particle;}

};


#endif

