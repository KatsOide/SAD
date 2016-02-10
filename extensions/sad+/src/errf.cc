#include <cmath>
using std::fabs;
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::exp;
using std::log;

#include <Complex.h>
#define MWFLT	1
#define MREAL	3
//*     IMPORTANT: MWNAM must be an integer multiple of MWFLT.
#define MCWRD	4
#define MCNAM	16
#define MWNAM	4		//MWNAM = MCNAM / MCWRD)
#define CC	1.12837916709551
#define XLIM	5.33
#define YLIM	4.29

// ERRF *************************************************
Complex cerf(Complex x)
{
/*----------------------------------------------------------------------*
 * Purpose:                                                             *
 *   Modification of WWERF, double precision Complex error function,    * 
 *   written at CERN by K. Koelbig.                                     *
 * Input:                                                               *
 *   XX, YY    (real)    Argument to CERF.                              *
 * Output:                                                              *
 *   WX, WY    (real)    Function result.                               *
 * conversion to C++ by K.Ohmi  (94.3.10)                               *
 *----------------------------------------------------------------------*/

   double RX[33], RY[33];
   double X,Y;
   double Q,H,XL,XH,YH,TX,TY,TN,SX,SY,SAUX,WX,WY;
   int NC,NU,n;

   X = fabs(real(x));
   Y = fabs(imag(x));

   if(Y < YLIM && X < XLIM) {
      Q  = (1. - Y / YLIM) * sqrt(1. - X*X/(XLIM*XLIM));
      H  = 1. / (3.2 * Q);
      NC = 7 + (int)(23.0*Q);
      XL = pow(H,(double)(1 - NC));
      XH = Y + 0.5/H;
      YH = X;
      NU = 10 + (int)(21.0*Q);
      RX[NU] = 0.;
      RY[NU] = 0.;

      for(n=NU;n>0;n--) {	//DO 10 N = NU, 1, -1
          TX = XH + n * RX[n];
          TY = YH - n * RY[n];
          TN = TX*TX + TY*TY;
          RX[n-1] = 0.5 * TX / TN;
          RY[n-1] = 0.5 * TY / TN;
      }

      SX = 0.;
      SY = 0.;

      for(n=NC;n>0;n--) {	//DO 20 N = NC, 1, -1
	 SAUX = SX + XL;
	 SX = RX[n-1] * SAUX - RY[n-1] * SY;
	 SY = RX[n-1] * SY + RY[n-1] * SAUX;
	 XL = H * XL;
      }

      WX = CC * SX;
      WY = CC * SY;
   }
   else {
      XH = Y;
      YH = X;
      RX[0] = 0.;
      RY[0] = 0.;

      for(n=9;n>0;n--) {	//DO 30 N = 9, 1, -1
         TX = XH + n * RX[0];
         TY = YH - n * RY[0];
         TN = TX*TX + TY*TY;
         RX[0] = 0.5 * TX / TN;
         RY[0] = 0.5 * TY / TN;
      }

      WX = CC * RX[0];
      WY = CC * RY[0];
   }

   if(Y == 0.) WX = exp(-X*X);
   if(imag(x) < 0.) {
       WX =   2. * exp(Y*Y-X*X) * cos(2.*X*Y) - WX;
       WY = - 2. * exp(Y*Y-X*X) * sin(2.*X*Y) - WY;
       if(real(x) > 0.) WY = -WY;
   }
   else {
       if(real(x) < 0.) WY = -WY;
   }
   return Complex(WX,WY);
}

double erf_fl(double x)
{
   double t,z,ans;
   z=fabs(x);
   t=1./(1.+0.5*z);
   ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
       t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
       t*(-0.82215223+t*0.17087277)))))))));
   return x>=0.0 ? 1.-ans : ans-1.;
}

Complex erf(Complex x)
{
#if 0
   Complex z;
   z=1.-exp(-x*x)*cerf(Complex(0.,1.)*x);
#else
   Complex z(1.);
   z-=exp(-x*x)*cerf(Complex(0.,1.)*x);
#endif
   return z;
}
