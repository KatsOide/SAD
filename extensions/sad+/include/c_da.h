#ifndef C_DA_H
#define C_DA_H
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;

#include <Complex.h>
#include <dacpp.h>


class c_da    //:private a
{
   da *z;

public:
   c_da(void) {z=new da[2];}
   c_da(const c_da& x) {
      z=new da[2];
      z[0]=x.z[0]; 
      z[1]=x.z[1];
      
   }
   c_da(const da& x,const da& y) { z=new da[2]; z[0]=x; z[1]=y;}
   ~c_da(void){ delete [] z; }

   friend void DUnit(c_da& x,int i,double f) { 
      DUnit(x.z[0],i,f); DUnit(x.z[1],i,0.); }
   friend void DUnit(c_da& x,int i,Complex f) { 
      DUnit(x.z[0],i,real(f)); DUnit(x.z[1],i,imag(f)); }
   friend void NUnit(c_da& x,int i,double f,int j,double g) { 
      DUnit(x.z[0],i,f); DUnit(x.z[1],j,g); }

   c_da& operator=(const c_da& x) {
      z[0]=x.z[0]; z[1]=x.z[1];
      return(*this);
   }
   c_da& operator=(double f) {
      z[0]=f;  z[1]=0.;
      return(*this);
   }
   c_da& operator=(Complex x) {
      z[0]=real(x); z[1]=imag(x);
      return(*this);
   }
   c_da& operator+(void) {return(*this);}
   c_da operator-(void) {return (c_da(-z[0], -z[1]));}
   c_da& operator+=(const c_da& x) {
      z[0]+=x.z[0];    z[1]+=x.z[1];
      return(*this);
   }
   c_da& operator-=(const c_da& x){
      z[0]-=x.z[0];    z[1]-=x.z[1];
      return(*this);
   }
   c_da operator*=(const c_da& y) {
      return(c_da(z[0]*y.z[0]-z[1]*y.z[1], z[0]*y.z[1]+z[1]*y.z[0]));
   }
      
   
   friend c_da operator+(const c_da& x,const c_da& y) {
      return(c_da(x.z[0]+y.z[0], x.z[1]+y.z[1]));
   }
   friend c_da operator-(const c_da& x,const c_da& y) { 
      return(c_da(x.z[0]-y.z[0], x.z[1]-y.z[1]));
   }
   
   friend c_da operator*(const c_da& x,const c_da& y) {
      return(c_da(x.z[0]*y.z[0]-x.z[1]*y.z[1], x.z[0]*y.z[1]+x.z[1]*y.z[0]));
   }
   friend c_da operator*(const da& f,const c_da& x) {
      return(c_da(f*x.z[0], f*x.z[1]));
   }
   friend c_da operator*(const c_da& x,const da& f) {
      return(c_da(f*x.z[0], f*x.z[1]));
   }
   friend c_da operator*(double f,const c_da& x) {
      return(c_da(f*x.z[0], f*x.z[1]));
   }
   friend c_da operator*(const c_da& x,double f) {
      return(c_da(f*x.z[0], f*x.z[1]));
   }
   friend c_da operator*(Complex f,const c_da& x) {
    return(c_da(real(f)*x.z[0]-imag(f)*x.z[1], real(f)*x.z[1]+imag(f)*x.z[0]));
   }
   friend c_da operator*(const c_da& x,Complex f) {
    return(c_da(real(f)*x.z[0]-imag(f)*x.z[1], real(f)*x.z[1]+imag(f)*x.z[0]));
   }
   friend c_da mul(const c_da& x,const c_da& y,int No,int Mx=0, int My=0) {
      c_da w;
      w.z[0]=mul(x.z[0],y.z[0],No,Mx,My)-mul(x.z[1],y.z[1],No,Mx,My);
      w.z[1]=mul(x.z[1],y.z[0],No,Mx,My)+mul(x.z[0],y.z[1],No,Mx,My);
      return w;
   }
   friend c_da muln(const c_da& x,const c_da& y,int No) {
      c_da w;
      w.z[0]=muln(x.z[0],y.z[0],No)-muln(x.z[1],y.z[1],No);
      w.z[1]=muln(x.z[1],y.z[0],No)+muln(x.z[0],y.z[1],No);
      return w;
   }
   
   friend c_da operator/(const c_da& z,const da& x) {
      return( c_da(z.z[0]/x , z.z[1]/x));
   }
   friend c_da operator/(double f,const c_da& z) {
      return( f*conj(z)/norm(z) );
   }
   friend c_da operator/(const c_da& z,double f) {
      return (z*(1./f));
   }
   
   friend da real(const c_da& x) { return da(x.z[0]);}
   friend da imag(const c_da& x) { return da(x.z[1]);}
   friend c_da conj(const c_da& x) { 
      c_da z1;
      z1.z[0]=x.z[0];
      z1.z[1]=x.z[1]*(-1.);
      return z1;
   }              //return(c_da(z.z[0], -z.z[1]));}
   friend da norm(const c_da& x) { return(x.z[0]*x.z[0]+x.z[1]*x.z[1]); }
   
   friend c_da pow(const c_da&, int);
   friend c_da exp(const c_da&);
   
   friend c_da erf(const c_da&);
   friend c_da cerf(const c_da&);
   friend c_da cerf_odd(const c_da&);
   
   void daprint(const char* s) {
      cout << "Complex da print out : " << s << endl;
      z[0].daprint("Real part");
      z[1].daprint("Imaginary part"); }
   friend c_da diag_part(const c_da& x) { 
      return c_da(diag_part(x.z[0]),diag_part(x.z[1]));}
   friend void xprint(const char* s,c_da& x) { x.daprint(s); }
   friend ostream& operator<<(ostream&,c_da&);
   void ci(int i,Complex f) { 
      z[0].ci(i,real(f)); z[1].ci(i,imag(f)); }
   void ci(int i,double f) { 
      z[0].ci(i,f); z[1].ci(i,0.); }
   friend Complex vget(const c_da& x,int i) {
      return(Complex(vget(x.z[0],i),vget(x.z[1],i)));}
   friend c_da p_gen_f(const double* rmu,const c_da& x);
};
#endif
