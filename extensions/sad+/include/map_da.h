#ifndef MAP_DA_H
#define MAP_DA_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_da.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <iostream>
using std::istream;
using std::ostream;

#include <dacpp.h>
#include <c_da.h>
#include <matrix.h>
#include <lin_map.h>
#include <map_double.h>
#include <track.h>

class map_c_da;

class map_da    //:private da
{
   da *m;
public:
   map_da(void) {
      m=new da[N_cv];
      //cout << " map constructed " << m << '\n'; 
   }
   map_da(const map_da&);
   map_da(const matrix&);

   ~map_da(void){ 
      //cout << "map destracted " << m << '\n'; 
      delete [] m;}

   void dBase(void);
   void dBase(double);
   void dBase(const double*);
   void NBase(void);
   void print(const char*);
   friend ostream& operator<<(ostream&,const map_da&);
   friend istream& operator>>(istream&,map_da&);
   friend void is_symplectic(const map_da&);
   map_da& operator=(const map_da&);
   map_da& operator=(const matrix&);
   map_da& operator=(double);
   map_da operator+=(const map_da&);
   map_da operator-=(const map_da&);
   da& operator[](int i) { return m[i];}
   friend matrix lin_da(const map_da&);

   friend da concatenate(const da&,const map_da&);
   friend da concatenate(const da&,const map_da&,int);
   friend map_da concatenate(const map_da&,const map_da&);
   friend map_da concatenate(const map_da&,const map_da&,int);
   friend map_double concatenate(const map_da&,const map_double&);
   friend pBeam concatenate(const map_da&,const pBeam&);

   friend map_c_da concatenate(const map_da&,const map_c_da&);
   friend map_c_da concatenate(const map_da&,const map_c_da&,int);

   map_double mapping(const map_double& x) {
      return concatenate(*this,x);}
   pBeam mapping(const pBeam& x) {
      return concatenate(*this,x);}
   
   friend map_da operator*(const map_da& x,const map_da& y) 
      { return(map_da(concatenate(x,y))); }
   friend map_da oldmul(const map_da&,const map_da&);
   friend map_da oldmul(const map_da&,const map_da&,int);
   friend map_da operator*(matrix&,const map_da&);
   friend map_da operator*(const map_da&,const matrix&);
   friend map_da operator*(const lin_map& x,const map_da& y)
     {return ((matrix&) x*y);}
   friend map_da operator*(const map_da& x,const lin_map& y)
     {return (x*(matrix&) y);}
   
   friend map_da Sym_trans(const map_da&,const matrix&);
   //Lie map
   friend da poi_itg(const map_da&);
   
   friend map_da lie_exp(const da& f,const map_da& x);
   friend map_da lie_exp(const da& f);

};


#endif
