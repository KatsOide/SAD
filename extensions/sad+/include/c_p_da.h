#ifndef C_P_DA_H
#define C_P_DA_H
#include <iostream>
using std::ostream;
using std::cout;

#include <Complex.h>
#include <dacpp.h>
#include <c_da.h>
#include <p_da.h>


class c_p_da    //:private a
{
   p_da *z;

public:
   c_p_da(void) {z=new p_da[2];}
   c_p_da(const c_p_da& x) {
      z=new p_da[2];
      z[0]=x.z[0]; 
      z[1]=x.z[1];
      
   }
   c_p_da(const p_da& x,const p_da& y) { z=new p_da[2]; z[0]=x; z[1]=y;}
   ~c_p_da(void){ delete [] z; }

   friend void DUnit(c_p_da& x,int i,double f) { 
      DUnit(x.z[0],i,f); DUnit(x.z[1],i,0.); }
   friend void DUnit(c_p_da& x,int i,Complex f) { 
      DUnit(x.z[0],i,real(f)); DUnit(x.z[1],i,imag(f)); }
   friend void NUnit(c_p_da& x,int i,double f,int j,double g) { 
      DUnit(x.z[0],i,f); DUnit(x.z[1],j,g); }

   c_p_da& operator=(const c_p_da& x) {
      z[0]=x.z[0]; z[1]=x.z[1];
      return(*this);
   }
   c_p_da& operator=(double f) {
      z[0]=f;  z[1]=0.;
      return(*this);
   }
   c_p_da& operator=(Complex x) {
      z[0]=real(x); z[1]=imag(x);
      return(*this);
   }
   c_p_da& operator+(void) {return(*this);}
   c_p_da operator-(void) {return (c_p_da(-z[0], -z[1]));}
   c_p_da& operator+=(const c_p_da& x) {
      z[0]+=x.z[0];    z[1]+=x.z[1];
      return(*this);
   }
   c_p_da& operator-=(const c_p_da& x){
      z[0]-=x.z[0];    z[1]-=x.z[1];
      return(*this);
   }
   c_p_da operator*=(const c_p_da& y) {
      return(c_p_da(z[0]*y.z[0]-z[1]*y.z[1], z[0]*y.z[1]+z[1]*y.z[0]));
   }
      
   friend p_da real(const c_p_da& x) { return p_da(x.z[0]);}
   friend p_da imag(const c_p_da& x) { return p_da(x.z[1]);}
   friend c_p_da conj(const c_p_da& x) { 
      c_p_da z1;
      z1.z[0]=x.z[0];
      z1.z[1]=x.z[1]*(-1.);
      return z1;
   }              //return(c_da(z.z[0], -z.z[1]));}
   friend p_da norm(const c_p_da& x) { return(x.z[0]*x.z[0]+x.z[1]*x.z[1]); }
   
   friend c_p_da operator+(const c_p_da& x,const c_p_da& y) {
      return(c_p_da(x.z[0]+y.z[0], x.z[1]+y.z[1]));
   }
   friend c_p_da operator-(const c_p_da& x,const c_p_da& y) { 
      return(c_p_da(x.z[0]-y.z[0], x.z[1]-y.z[1]));
   }
   
   friend c_p_da operator*(const c_p_da& x,const c_p_da& y) {
      return(c_p_da(x.z[0]*y.z[0]-x.z[1]*y.z[1], x.z[0]*y.z[1]+x.z[1]*y.z[0]));
   }
   friend c_p_da operator*(const p_da& f,const c_p_da& x) {
      return(c_p_da(f*x.z[0], f*x.z[1]));
   }
   friend c_p_da operator*(const c_p_da& x,const p_da& f) {
      return(c_p_da(f*x.z[0], f*x.z[1]));
   }
   friend c_p_da operator*(double f,const c_p_da& x) {
      return(c_p_da(f*x.z[0], f*x.z[1]));
   }
   friend c_p_da operator*(const c_p_da& x,double f) {
      return(c_p_da(f*x.z[0], f*x.z[1]));
   }
   friend c_p_da operator*(Complex f,const c_p_da& x) {
    return(c_p_da(real(f)*x.z[0]-imag(f)*x.z[1], real(f)*x.z[1]+imag(f)*x.z[0]));
   }
   friend c_p_da mul(const c_p_da& x,const c_p_da& y,int No,int Mx,int My) {
      c_p_da w;
      w.z[0]=da_mul(x.z[0],y.z[0],No,Mx,My)-da_mul(x.z[1],y.z[1],No,Mx,My);
      w.z[1]=da_mul(x.z[1],y.z[0],No,Mx,My)+da_mul(x.z[0],y.z[1],No,Mx,My);
      return w;
   }
   friend c_p_da muln(const c_p_da& x,const c_p_da& y,int No) {
      c_p_da w;
      w.z[0]=muln(x.z[0],y.z[0],No)-muln(x.z[1],y.z[1],No);
      w.z[1]=muln(x.z[1],y.z[0],No)+muln(x.z[0],y.z[1],No);
      return w;
   }
   
   friend c_p_da operator/(const c_p_da& z,const da& x) {
      return( c_p_da(z.z[0]/x , z.z[1]/x));
   }
   friend c_p_da operator/(const c_p_da& z,const p_da& x) {
      return( c_p_da(z.z[0]/x , z.z[1]/x));
   }
   friend c_p_da operator/(double f,const c_p_da& z) {
      return( f*conj(z)/norm(z) );
   }
   friend c_p_da operator/(const c_p_da& z,double x) {
      return( c_p_da(z.z[0]/x , z.z[1]/x));
   }
   
   
   friend c_p_da pow(const c_p_da&, int);
   friend c_p_da erf(const c_p_da&);
   friend c_p_da cerf(const c_p_da&);
   
   friend ostream& operator<<(ostream& s,c_p_da& x) {
      cout << "\n Complex da print out \n\n Real part \n" 
	   << x.z[0] << "\n Imaginary part\n" << x.z[1];
      return s;}
   void ci(int i,Complex f) { 
      z[0].ci(i,real(f)); z[1].ci(i,imag(f)); }
   c_da co(int i) {
      return(c_da(z[0].co(i),z[1].co(i)));}
   friend c_p_da p_gen_f(const double* rmu,const c_p_da& x);
};
#endif
