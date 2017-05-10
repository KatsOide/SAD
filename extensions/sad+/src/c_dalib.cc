#include <iostream>
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>

#include <Complex.h>
#include <phys_const.h>
#include <dacpp.h>
#include <c_da.h>

#define N_dim_max 20       // affect only erf and cerf table

Complex erf(Complex);

ostream& operator<<(ostream& s,c_da& x) 
{
  s << "RealPart =" 
    << x.z[0] << "\nImagPart=" << x.z[1];
  return s;
}

c_da pow(const c_da& x,int n)
{
   c_da z=x;
   if(n==0) z=1.;
   for(int i=1;i<n;i++) {z=z*x;}
   return z;
}

c_da exp(const c_da& x)
{
   da a=exp(real(x));
   return(c_da(a*cos(imag(x)), a*sin(imag(x))));
}

c_da erf(const c_da& x)
{
   Complex z0,a[N_dim_max];
   c_da xf,z;
   int i;
   
   if(N_ord>=N_dim_max) {
    cout<< " Nord is larger than max value (cerf) "
	<< N_ord << '>' << N_dim_max;
    exit(1);
   }
   xf=0.;
   
   z0=vget(x,0);
   z=x;
   z.ci(0,0.);

   if(z0!=0.) {
   Complex e1;
   e1=2.*exp(-z0*z0)/sqr_pi;
   a[0]=erf(z0);
   //a[i]=H(i-1)
   a[1]=1.;
   a[2]=2.*z0;
   for(i=3;i<=N_ord;i++) {
      a[i]=2.*z0*a[i-1]-2.*(i-2)*a[i-2];
   }
   // a[n]=erf(n)
   for(i=1;i<=N_ord;i++) {
      if(i%2) a[i]=a[i]*e1; else a[i]=-a[i]*e1;
   }
   if(N_ord==0) {xf=a[0];}
   else {
     Complex f=a[N_ord]/a[N_ord-1]/(double)N_ord;
     xf=z*f;
     xf.ci(0,1.);
     for(i=1;i<N_ord;i++) {
        f=a[N_ord-i]/a[N_ord-i-1]/(double)(N_ord-i);
        xf=mul(xf,z,i+1,0,1);
        //xf=xf*z;
        xf=xf*f;
        xf.ci(0,1.);
     }
     xf=xf*a[0];
   }   
   
   
   }
   else {
   z=x*x;

   int No2=(N_ord-1)/2;
  //   z  : Expanding variable

   double f=-(double)(2*No2-1)/((double)((2*No2+1)*No2));
   xf=f*z;
   xf.ci(0,1.);
   for(i=1;i<No2;i++) {
      //      xf=mul(xf,z,2*(i+1),0,2);
      xf=xf*z;
      f=-(double)(2*(No2-i)-1)/((double)((2*(No2-i)+1)*(No2-i)));
      xf=xf*f;
      xf.ci(0,1.);
   }
   xf=xf*x*(2./sqr_pi);
   }
   return xf;
}

c_da cerf(const c_da&x)
{
   c_da xf,z;
   xf=0.;
   
   z=Complex(0.,1.)*x;
   double a[N_dim_max];
   int i;
   if(N_ord>=N_dim_max) {
    cout<< " Nord is larger than max value (cerf) "
	<< N_ord << '>' << N_dim_max;
    exit(1);
   }
   if(vget(x,0)!=0.) {
      cout << " x.v[0] must be zero in this version \n";
   }
   a[0]=1.;
   a[1]=sqr_pi*0.5;
   for(i=2;i<=N_ord;i+=2) {a[i]=a[i-2]*(double)i/2.;}
   for(i=3;i<=N_ord;i+=2) {a[i]=a[i-2]*(double)i/2.;}
   
  if(N_ord==0) {xf=1.;}
  else if(N_ord==1) { xf=z/sqr_pi*2.;}
  else {
    double f=a[N_ord-1]/a[N_ord];
    xf=z*f;
    xf.ci(0,1.);
    for(i=1;i<N_ord;i++) {
       f=a[N_ord-i-1]/a[N_ord-i];
       xf=mul(xf,z,i+1,0,1);
       //xf=xf*z;
       xf=xf*f;
       xf.ci(0,1.);
    }
  }   
  return xf;
}

c_da cerf_odd(const c_da&x)
{
   c_da xf,z,z2;
   xf=0.;
   
   z=Complex(0.,1.)*x;
   double a[N_dim_max];
   int i;
   if(N_ord>=N_dim_max) {
    cout<< " Nord is larger than max value (cerf) "
	<< N_ord << '>' << N_dim_max;
    exit(1);
   }
   if(vget(x,0)!=0.) {
      cout << " x.v[0] must be zero in this version \n";
   }
   a[0]=0.;
   a[1]=sqr_pi*0.5;
   for(i=2;i<=N_ord;i+=2) {a[i]=0.;}
   for(i=3;i<=N_ord;i+=2) {a[i]=a[i-2]*(double)i/2.;}
   
  if(N_ord==0) {xf=0.;}
  else if(N_ord==1) { xf=z/sqr_pi*2.;}
  else {
     int No;
     int nn=N_ord%2;
     if(nn) No=N_ord; else No=N_ord-1;
     z2=z*z;
    double f=a[No-2]/a[No];
    xf=z2*f;
    xf.ci(0,1.);
    for(i=2;i<No-2;i+=2) {
       f=a[No-i-2]/a[No-i];
       if(nn) xf=mul(xf,z2,i+2,0,2); else xf=mul(xf,z2,i+3,0,2);
       //xf=xf*z2;
       xf=xf*f;
       xf.ci(0,1.);
    }
    xf=xf*z/a[1];
  }   
  return xf;
}

