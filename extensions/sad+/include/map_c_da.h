#ifndef MAP_C_DA_H
#define MAP_C_DA_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_c_da.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <dacpp.h>
#include <c_da.h>
class map_c_da    //:private da
{
   c_da *m;
public:
   map_c_da(void) {
      m=new c_da[N_cv];
      //      cout << " map_c_da alocated " << m << '\n'; 
   }
   map_c_da(const map_c_da&);
   ~map_c_da(void){ delete [] m;}

   void dBase(void);
   void dBase(double);
   void dBase(const double*);
   void NBase(void);
   void NiBase(void);
   void print(const char*);
   friend ostream& operator<<(ostream&,const map_c_da&);
   friend void is_symplectic(const map_c_da&);
   map_c_da& operator=(const map_c_da&);
   map_c_da& operator=(double);
   map_c_da operator+=(const map_c_da&);
   map_c_da operator-=(const map_c_da&);
   c_da& operator[](int i) { return m[i];}

   friend c_da concatenate(const da&,const map_c_da&);
   friend c_da concatenate(const da&,const map_c_da&,int);
   friend c_da concatenate(const c_da&,const map_c_da&);
   friend c_da concatenate(const c_da&,const map_c_da&,int);
   friend map_c_da concatenate(const map_da&,const map_c_da&);
   friend map_c_da concatenate(const map_da&,const map_c_da&,int);
   friend map_c_da operator*(const map_da& x,const map_c_da& y) 
      { return(map_c_da(concatenate(x,y))); }
   friend map_c_da concatenate(const map_c_da&,const map_c_da&);
   friend map_c_da concatenate(const map_c_da&,const map_c_da&,int);
   friend map_c_da operator*(const map_c_da& x,const map_c_da& y) 
      { return(map_c_da(concatenate(x,y))); }
   
   //Lie map
   
   friend map_c_da lie_exp(const da&,const map_c_da&);
   friend da lie_exp(const da&,const da&);
   friend map_c_da fac_map(const c_da&);
   friend map_c_da faci_map(const c_da&);
   friend da fac_drg(const map_c_da&);
   friend da fac_drgi(const map_c_da&);
   friend da fac_lie(const map_c_da&);
   friend da fac_drg(const map_c_da&,map_c_da&);
   friend da fac_lie(const map_c_da&,map_c_da&);
   friend da poi_itg(const map_c_da&);
   
};


#endif
