#ifndef MAP_DOUBLE_H
#define MAP_DOUBLE_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_double.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <iostream>
using std::ostream;

#include <dacpp.h>

class map_da;

class map_double    //:private da
{
   double *m;
   int sv;
public:
   map_double(void) {
      m=new double[N_cv];
   }
   map_double(const map_double&);
   ~map_double(void){ delete [] m;}
   friend ostream& operator<<(ostream&,map_double&);

   void dBase(void);
   void dBase(double);
   void dBase(const double*);
   void NBase(void);
   void print(const char*);
   int fsv(void) {return sv;}
   void setsv(int xsv) {sv=xsv;}
   friend void is_symplectic(const map_double&);
   friend int is_survive(map_double&);
   map_double& operator=(const map_double&);
   map_double& operator=(double);
   map_double& operator=(double*);
   map_double operator+=(const map_double&);
   map_double operator-=(const map_double&);
   friend map_double operator+(const map_double&,const map_double&);
   friend map_double operator-(const map_double&,const map_double&);
   friend map_double operator*(double,const map_double&);
   friend map_double operator*(const map_double&,double);
   friend double concatenate(const da&,const map_double&);
   friend map_double concatenate(const map_da&,const map_double&);
   double& operator[](int i) {return m[i];}

};


#endif
