#include <Complex.h>
c_da p_gen_f(const double* rmu,const c_da& x)
{
  c_da y; 
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A,j;
  Complex xi,x2,x3;

  y=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	   C1M=*(Id1+i1A);
	   ivcal(C1M,Nord,Nvr2,iv1);
	   double a=0.;
	   for(j=0;j<Nvr2;j++) 
	      a+=(iv1[j]-iv2[j])*rmu[j];
	   if(fabs(a)>1.e-5) {
//	      xi*=1./(1.-exp(Complex(0.,a)));
//	      xi*=Complex((1.-cos(a)),sin(a))*0.5/(1.-cos(a));
//	      xi*=Complex(0.5,0.5*sin(a)/(1.-cos(a)));
//	      xi*=Complex(0.5,0.5*(1.+cos(a))/sin(a));
	      double cos_a=cos(a);
	      double sin_a=sin(a);
	      if(cos_a<0) {
		 xi*=Complex(0.5,0.5*sin_a/(1.-cos_a));
	      } else {
		 xi*=Complex(0.5,0.5*(1.+cos_a)/sin_a);
	      }
	      y.ci(i,xi);
	   }
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return y;
}

