#ifndef DACPP_H
#define DACPP_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               dacpp.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define Nresist 6

#include <iostream>
using std::istream;
using std::ostream;

int ccal(int,int,int*);
void ivcalp(int,int,int,int*);
void ivcal(int,int,int,int*);
int icomb(int,int);
void da_init(char*,int,int);
void da_mul(double*,double*,double*);

extern int N_ord,N_var,L_vec,N_cv;

void InitializeDifferentialAlgebra
(int No=N_ord,int Nv=N_var,int Npo=0,int Npv=0);


class da
{
   double* v;
public:
   da(void);
   da(const da&);
   da(double* z) {v=z;}

   friend ostream& operator<<(ostream&,da&);
   friend istream& operator>>(istream&,da&);
   void daprint(const char*);
   void daprint(int);
   void NFPrint(ostream&);
   friend void xprint(const char* s,da& x) { x.daprint(s); }
   friend void pprint(da&);
   double co(int i) { return (v[i]);}
   friend double vget(const da& x,int i) { return (x.v[i]);}
   void ci(int i,double f) { v[i]=f; }
   double lin_da(int);
   
   da& operator=(const da&);
   da& operator=(double);
   da& operator=(const char*);
   da& operator+(void) {return(*this);}
   da  operator-(void);

   da& operator+=(const da&);
   da& operator-=(const da&);
   da& operator+=(double x) {v[0]+=x; return *this;}
   da& operator-=(double x) {v[0]-=x; return *this;}
   da& operator*=(const da& x) { da_mul(v,x.v,v); return *this; }
   da& operator/=(const da&);
   friend da operator+(const da&,const da&);
   friend da operator+(const da&,double);
   friend da operator+(double,const da&);
   friend da operator-(const da&,const da&);
   friend da operator-(const da&,double);
   friend da operator-(double,const da&);

   friend da operator*(const da& y,const da& x) { 
      da z; da_mul(y.v,x.v,z.v); return z;}
   friend da operator*(double,const da&);
   friend da operator*(const da&,double);
   friend da operator/(double,const da&);
   friend da operator/(const da&,double);
   friend da operator/(const da&,const da&);
   
   int operator==(double);
   int operator!=(double);

   friend void da_mul(double*,double*,double*);
   void mul(const da& x,const da& y) { da_mul(x.v,y.v,v); }
   void mul(const da&,double);
   void mul(double,const da&);
   friend da mul(const da&,const da&,int,int,int);
   friend da muln(const da&,const da&,int);

   void valset(int i,double f) { *(v+i)=f;}
   void valset(const double*);
   void msk(const da&,const int);
   friend da msk(const da&,int);
   void varmsk(const da&,const int);
   double coef(const char*);
   friend int min_ord(const da&); 
   friend da diag_part(const da&);

   void dBase(int, double);
   void Dcin(int,double);
   void dBase(double*);
   void test(void);
   friend void DUnit(da& x,int i,double f) { x.dBase(i,f); }
   friend void NUnit(da& x,int i,double f) { x.dBase(i,f); }
	
   ~da(void){ delete [] v;}

   friend da dpow(const da&,int,int);
   friend da dif(const da&,const int);
   friend da itg(const da&,const int);
   //friend da pow(const da&,const int,const int);
   friend da sin(const da&);
   friend da cos(const da&);
   friend da sinh(const da&);	
   friend da cosh(const da&);
   friend da exp(const da&);
   friend da log(const da&);
   friend da sqrt(const da&);
   friend da asin(const da&);
   friend da atan(const da&);
   friend da pow(const da&,double);
   
   friend da poi(const da&,const da&);
   friend da line_itg(const da&);
};


//
//
//  ----------------------------------------------------------------------
//     inline functions
//  ----------------------------------------------------------------------
//


inline void da::test(void)
{
   for(int i=0;i<L_vec;i++) *(v+i)=double(i+1);
}

inline void da::valset(const double* f)
{
   for(int i=0;i<L_vec;i++) *(v+i)=*(f+i);
}


inline da& da::operator+=(const da& x)
{
   for(int i=0;i<L_vec;i++) v[i]+=x.v[i];
   return *this;
}

inline da& da::operator-=(const da& x)
{
   for(int i=0;i<L_vec;i++) v[i]-=x.v[i];
   return *this;
}

inline da& da::operator=(const da& x)
{
  for(int i=0;i<L_vec;i++) *(v+i)=*(x.v+i);
  return *this;
}
inline da& da::operator=(double f)
{
  for(int i=1;i<L_vec;i++) *(v+i)=0.;
  *v=f;
  return *this;
}
inline da operator+(const da& y,const da& x)
{
  da z;
  for(int i=0;i<L_vec;i++) z.v[i]=(*(y.v+i))+x.v[i];
  return z;
}
inline da operator+(const da& y,double f)
{
  da z;
  for(int i=1;i<L_vec;i++) z.v[i]=(*(y.v+i));
  z.v[0]=y.v[0]+f;
  return z;
}
inline da operator+(double f,const da& x)
{
  da z;
  for(int i=1;i<L_vec;i++) z.v[i]=(*(x.v+i));
  z.v[0]=x.v[0]+f;
  return z;
}
inline da operator-(const da& y,const da& x)
{
  da z;
  for(int i=0;i<L_vec;i++) z.v[i]=(*(y.v+i))-x.v[i];
  return z;
}
inline da operator-(const da& y,double f)
{
  da z;
  for(int i=1;i<L_vec;i++) z.v[i]=(*(y.v+i));
  z.v[0]=y.v[0]-f;
  return z;
}
inline da operator-(double f,const da& x)
{
  da z;
  for(int i=1;i<L_vec;i++) z.v[i]=-(*(x.v+i));
  z.v[0]=f-x.v[0];
  return z;
}

inline da operator*(const da& x,double f)
{
  da z;
  for(int i=0;i<L_vec;i++) z.v[i]=f*x.v[i];
  return z;
}
inline da operator*(double f,const da& x)
{
  da z;
  for(int i=0;i<L_vec;i++) z.v[i]=f*x.v[i];
  return z;
}



inline da operator/(const da& x,double f)
{
   da z;
   double g=1./f;
   for(int i=0;i<L_vec;i++) z.v[i]=x.v[i]*g;
   return z;
}
inline da operator/(const da& x,const da& y)
{
   da z;
   z=(1./y)*x;
   return z;
}

inline da da::operator-(void)
{
   da z;
   for(int i=0;i<L_vec;i++) z.v[i]=-v[i];
   return z;
}

inline int da::operator==(double f)
{
   for(int i=0;i<L_vec;i++) {
      if(v[i]!=f) return(1);
   }
   return(0);
}

inline int da::operator!=(double f)
{
   for(int i=0;i<L_vec;i++) {
      if(v[i]!=f) return(0);
   }
   return(1);
}   


inline void da::mul(const da& x,double f)
{
  for(int i=0;i<L_vec;i++) v[i]=f*x.v[i];
}

inline void da::mul(double f,const da& x)
{
  for(int i=0;i<L_vec;i++) v[i]=f*x.v[i];
}



inline da msk(const da& x,int i)
{
   da z;
   z.msk(x,i);
   return z;
}
#endif
