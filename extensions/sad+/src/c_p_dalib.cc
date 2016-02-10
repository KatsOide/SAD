#include <iostream>
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>

#include <dacpp.h>
#include <c_da.h>
#include <c_p_da.h>
#include <phys_const.h>

#define N_dim_max 50       // affect only cerf table

c_p_da pow(const c_p_da& x,int n)
{
   c_p_da z=x;
   if(n==0) z=1.;
   for(int i=1;i<n;i++) {z=z*x;}
   return z;
}

c_p_da exp(const c_p_da& x)
{
   p_da a=exp(real(x));
   return(c_p_da(a*cos(imag(x)), a*sin(imag(x))));
}

c_p_da erf(const c_p_da& x)
{
   c_p_da xf,z;
   xf=0.;

   z=x*x;

   int No2=(N_ord-1)/2;
  //   z  : Expanding variable

  xf=-(double)(2*No2-1)*z/((double)((2*No2+1)*No2));
  xf.ci(0,1.);
  for(int i=1;i<No2;i++) {
     xf=mul(xf,z,i,0,1);
     xf=-(double)(2*(No2-i)-1)*z/((double)((2*(No2-i)+1)*(No2-i)));
     xf.ci(0,1.);
  }
  return xf*x;
}

c_p_da cerf(const c_p_da&x)
{
   c_p_da xf,z;
   xf=0.;
   
   z=Complex(0.,1.)*x;
   double a[N_dim_max];
   int i;
   if(N_ord>=N_dim_max) {
    cout<< " Nord is larger than max value (cerf) "
	<< N_ord << '>' << N_dim_max;
    exit(1);
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
       xf=xf*f;
       xf.ci(0,1.);
    }
  }   
  return xf;
}
