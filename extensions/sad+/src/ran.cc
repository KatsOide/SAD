#define IA 16807
#define IM 2147483647
#define AM (1./IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.-EPS)
#define pi 3.14159265358979324

#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
using std::log;

static long idum=17;
static double pi2;
static double f;

float ran1()
{
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   float temp;
   
   if(idum<=0 || !iy) {
      if(-idum < 1) idum=1;
      else idum=-idum;
      for(j=NTAB+7;j>=0;j--) {
	 k=idum/IQ;
	 idum=IA*(idum-k*IQ)-IR*k;
	 if(idum < 0) idum+=IM;
	 if(j<NTAB) iv[j]=idum;
      }
      iy=iv[0];
   }
   k=idum/IQ;
   idum=IA*(idum-k*IQ)-IR*k;
   if(idum < 0) idum+=IM;
   j=int(iy/NDIV);
   iy=iv[j];
   iv[j]=idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}

void SetRandum(int iseed)
{
  idum=iseed;
  pi2=2.*pi;
  f=pi2/4294967296.0;
}

float rgauss()
{
  float ran1();
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;

  if(iset==0) {
    do {
      v1=2.*ran1()-1.;
      v2=2.*ran1()-1.;
      rsq=v1*v1+v2*v2;
    } while(rsq>=1. || rsq==0.);
    fac=sqrt(-2.*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}


double tran()
{
  double ran;
  idum=idum*48828125;
  ran=idum/4294967296.0+.5;
  return ran;
}


double tgauss()
{
  double ran;
  static int flag=0;
  static double r,phi;

  if(idum==0) {
    flag=0;
    return 0;
  }
  if(flag) {
    ran=r*sin(phi);}
  else {
    idum=idum*48828125;
    phi=idum*f;
    idum=idum*48828125;
    r=sqrt(-2.0*log(idum/4294967296.0+.5));
    ran=r*cos(phi);
  }
  flag=!flag;
  return ran;
}


