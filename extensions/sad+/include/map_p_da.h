#ifndef MAP_P_DA_H
#define MAP_P_DA_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               map_p_da.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <cstdio>

#include <dacpp.h>
#include <p_da.h>
class map_p_da    //:private da
{
   p_da *m;
public:
   map_p_da(void) {
      m=new p_da[N_cv];
      std::printf(" map_p_da alocated %p \n",(void*)m); 
   }
   map_p_da(const map_p_da&);
   ~map_p_da(void){ delete [] m;}

   void dBase(void);
   void dBase(double);
   void dBase(const double*);
   void NBase(void);
   void print(const char*);
   friend void is_symplectic(const map_p_da&);
   map_p_da& operator=(const map_p_da&);
   map_p_da& operator=(double);
   map_p_da operator+=(const map_p_da&);
   map_p_da operator-=(const map_p_da&);
   p_da& operator[](int i) { return m[i];}
   friend ostream& operator<<(ostream&,const map_p_da&);

/*   friend da concatenate(const da&,const map_p_da&);
   friend da concatenate(const da&,const map_p_da&,int);
   friend map_p_da concatenate(const map_p_da&,const map_p_da&);
   friend map_p_da concatenate(const map_p_da&,const map_p_da&,int);
*/
//   friend map_p_da operator*(const map_p_da& x,const map_p_da& y) 
//      { return(map_p_da(concatenate(x,y))); }
//   friend map_p_da oldmul(const map_p_da&,const map_p_da&);
//   friend map_p_da oldmul(const map_p_da&,const map_p_da&,int);
   
   //Lie map
   
/*   friend map_p_da lie_exp(const da&,const map_p_da&);
   friend da lie_exp(const da&,const da&);
   friend map_p_da fac_map_p_da(const da&);
   friend map_p_da faci_map_p_da(const da&);
   friend da fac_drg(const map_p_da&);
   friend da fac_drgi(const map_p_da&);
   friend da fac_lie(const map_p_da&);
   friend da fac_drg(const map_p_da&,map_p_da&);
   friend da fac_lie(const map_p_da&,map_p_da&);
   friend da poi_itg(const map_p_da&);
*/
};


#endif
