#ifndef P_DA_H
#define P_DA_H
#define Nresist 6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               p_da.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
using std::ostream;

#include <dacpp.h>
 
extern int N_ord,N_var,L_vec,N_cv;
extern int Np_ord,Np_var,Lp_vec;


class p_da
{
   da *pv;
public:
   p_da(void) {pv=new da[Lp_vec];}
   p_da(const p_da&);
   ~p_da(void){ delete [] pv;}

   p_da& operator=(const p_da& x) {
      for(int i=0;i<Lp_vec;i++) pv[i]=x.pv[i];
      return *this;
   }
   p_da& operator=(const da& x) {
      for(int i=1;i<Lp_vec;i++) pv[i]=0.;
      pv[0]=x;
      return *this;
   }
   p_da& operator=(const double f) {
      for(int i=1;i<Lp_vec;i++) pv[i]=0.;
      pv[0]=f;
      return *this;
   }
   
   friend ostream& operator<<(ostream&,const p_da&);
   da co(int i) { return pv[i]; }
   void ci(int i,da& x) { pv[i]=x;}
   void ci(int i,double f) { pv[i]=f;}
   void Dcin(int i,int j,double f) { pv[i].Dcin(j,f);}
   double lin_da(int i,int j) { return pv[i].lin_da(j);}
   friend void DUnit(p_da& x,int i,double f) { DUnit(x.pv[0],i,f); }
   friend void NUnit(p_da& x,int i,double f) { DUnit(x.pv[0],i,f); }

   p_da& operator+(void) {return(*this);}
   p_da  operator-(void);
   p_da& operator+=(const p_da&);
   p_da& operator-=(const p_da&);
   p_da& operator+=(double);
   p_da& operator-=(double);
   p_da& operator=(const char*);

   friend p_da operator+(const p_da&,const p_da&);
   friend p_da operator+(const p_da&,double);
   friend p_da operator+(double,const p_da&);
   friend p_da operator+(const p_da&,const da&);
   friend p_da operator+(const da&,const p_da&);
   friend p_da operator-(const p_da&,const p_da&);
   friend p_da operator-(const p_da&,double);
   friend p_da operator-(const da&,const p_da&);
   friend p_da operator-(const p_da&,const da&);
   friend p_da operator-(double,const p_da&);

   friend p_da operator*(const p_da&,const p_da&);
   friend p_da operator*(double,const p_da&);
   friend p_da operator*(const p_da&,double);
   friend p_da operator*(const da&,const p_da&);
   friend p_da operator*(const p_da&,const da&);
   friend p_da operator/(double,const p_da&);
   friend p_da operator/(const p_da&,double);
   friend p_da operator/(const p_da&,const da&);
   friend p_da operator/(const p_da&,const p_da&);

   friend p_da muln(const p_da&,const p_da&, int);
   friend p_da da_mul(const p_da&,const p_da&, int, int, int);

   void varmsk(const da&,int);

   friend p_da dif(const p_da&,int);
   friend p_da p_dif(const p_da&,int);
   friend p_da itg(const p_da&,int);
   friend p_da line_itg(const p_da&);
   friend p_da poi(const p_da&,const p_da&);
   //friend da pow(const p_da&,int,int);
   friend p_da sin(const p_da&);
   friend p_da cos(const p_da&);
   friend p_da sinh(const p_da&);	
   friend p_da cosh(const p_da&);
   friend p_da exp(const p_da&);
   friend p_da log(const p_da&);
   friend p_da sqrt(const p_da&);
   friend p_da asin(const p_da&);
   friend p_da atan(const p_da&);
   friend p_da pow(const p_da&,double);
};


//
//
//  ----------------------------------------------------------------------
//     inline functions
//  ----------------------------------------------------------------------
//

inline p_da p_da::operator-(void)
{
   p_da z;
   for(int i=0;i<Lp_vec;i++) z.pv[i]=-pv[i];
   return z;
}

inline p_da& p_da::operator+=(const p_da& x)
{
   for(int i=0;i<Lp_vec;i++) pv[i]+=x.pv[i];
   return *this;
}

inline p_da& p_da::operator-=(const p_da& x)
{
   for(int i=0;i<L_vec;i++) pv[i]-=x.pv[i];
   return *this;
}

inline p_da& p_da::operator+=(double x)
{
   pv[0]+=x;
   return *this;
}

inline p_da& p_da::operator-=(double x)
{
   pv[0]-=x;
   return *this;
}


inline p_da operator+(const p_da& y,const p_da& x)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=y.pv[i]+x.pv[i];
  return z;
}
inline p_da operator+(const p_da& y,double f)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=y.pv[i];
  z.pv[0]=y.pv[0]+f;
  return z;
}
inline p_da operator+(double f,const p_da& x)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=x.pv[i];
  z.pv[0]=x.pv[0]+f;
  return z;
}
inline p_da operator+(const p_da& x,const da& y)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=x.pv[i];
  z.pv[0]=x.pv[0]+y;
  return z;
}
inline p_da operator+(const da& x,const p_da& y)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=y.pv[i];
  z.pv[0]=y.pv[0]+x;
  return z;
}
inline p_da operator-(const p_da& y,const p_da& x)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=y.pv[i]-x.pv[i];
  return z;
}
inline p_da operator-(const p_da& x,const da& y)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=x.pv[i];
  z.pv[0]=x.pv[0]-y;
  return z;
}
inline p_da operator-(const da& x,const p_da& y)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=y.pv[i];
  z.pv[0]=x-y.pv[0];
  return z;
}
inline p_da operator-(const p_da& y,double f)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=y.pv[i];
  z.pv[0]=y.pv[0]-f;
  return z;
}
inline p_da operator-(double f,const p_da& x)
{
  p_da z;
  for(int i=1;i<Lp_vec;i++) z.pv[i]=-x.pv[i];
  z.pv[0]=f-x.pv[0];
  return z;
}

inline p_da operator*(const p_da& x,double f)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=f*x.pv[i];
  return z;
}
inline p_da operator*(double f,const p_da& x)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=f*x.pv[i];
  return z;
}
inline p_da operator*(const p_da& x,const da& f)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=f*x.pv[i];
  return z;
}
inline p_da operator*(const da& f,const p_da& x)
{
  p_da z;
  for(int i=0;i<Lp_vec;i++) z.pv[i]=f*x.pv[i];
  return z;
}


inline p_da operator/(const p_da& x,double f)
{
   p_da z;
   double g=1./f;
   for(int i=0;i<Lp_vec;i++) z.pv[i]=x.pv[i]*g;
   return z;
}
inline p_da operator/(const p_da& x,const da& f)
{
   p_da z;
   da g=1./f;
   for(int i=0;i<Lp_vec;i++) z.pv[i]=x.pv[i]*g;
   return z;
}

#endif

