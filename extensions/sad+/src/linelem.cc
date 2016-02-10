#include <cmath>
using std::fabs;
using std::sqrt;
using std::sin;
using std::cos;
using std::sinh;
using std::cosh;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <element.h>
#include <linelem.h>

#define nbslice 4

LBend::LBend(char* s) : Element(s)
{
  double dsh,sqr_K,sqr_G;

  phi_0=get_parm(s,"phi0");
  dphi=get_parm(s,"dphi");
  K1=get_parm(s,"K1");
  phi=phi_0+dphi;
 
  e1=get_parm(s,"e1");
  e2=get_parm(s,"e2");

  omega=phi_0-e1-e2;
  if(phi!=0.) rho=length/phi; else rho=0.;

  sin_e1=sin(e1);
  cos_e1=cos(e1);
  sin_e2=sin(e2);
  cos_e2=cos(e2);
  tan_e1=sin_e1/cos_e1;
  tan_e2=sin_e2/cos_e2;

  ds0=length/nbslice;
  dsh=ds0*0.5;

  if(phi!=0) G=(1./rho/rho+K1); else G=0.;
  if(G==0.) {
    a11=1.;
    a12=ds0;
    a11h=1.;
    a12h=dsh;
  }
  else if(G>0) {
    sqr_G=sqrt(G);
    a11=cos(sqr_G*ds0);
    a12=sin(sqr_G*ds0)/sqr_G;
    a11h=cos(sqr_G*dsh);
    a12h=sin(sqr_G*dsh)/sqr_G;
  }
  else {
    sqr_G=sqrt(-G);
    a11=cosh(sqr_G*ds0);
    a12=sinh(sqr_G*ds0)/sqr_G;
    a11h=cosh(sqr_G*dsh);
    a12h=sinh(sqr_G*dsh)/sqr_G;
  }

  if(K1==0.) {
    a33=1.;
    a34=ds0;
    a33h=1.;
    a34h=dsh;
  }
  else if(K1>0.) {
    sqr_K=sqrt(K1);
    a33=cosh(sqr_K*ds0);
    a34=sinh(sqr_K*ds0)/sqr_K;
    a33h=cosh(sqr_K*dsh);
    a34h=sinh(sqr_K*dsh)/sqr_K;
  }
  else {
    sqr_K=sqrt(-K1);
    a33=cos(sqr_K*ds0);
    a34=sin(sqr_K*ds0)/sqr_K;
    a33h=cos(sqr_K*dsh);
    a34h=sin(sqr_K*dsh)/sqr_K;
  }

  a16=(1.-a11)/rho/G;
  a51=-a12/rho;
  a52=-a16;
  a56=-(ds0-a12)/rho/rho/G;

  a16h=(1.-a11h)/rho/G;
  a51h=-a12h/rho;
  a52h=-a16h;
  a56h=-(dsh-a12h)/rho/rho/G;

   // cout<< "initialize Bend\n";
}

LBend::LBend(double* K) : Element(K)
{
  double dsh,sqr_K,sqr_G;

  phi_0=K[1];
  dphi=K[9];
  phi=phi_0+dphi;
  K1=K[6];
  if(K1!=0 && length>0) K1=K1/length;

  if(fabs(K[4])<1.e-3) K[5]=K[4];
  if(fabs(K[4]-pi/2.)<1.e-3) K[5]=K[4]-pi/2.;
  if(fabs(K[4]+pi/2.)<1.e-3) K[5]=K[4]+pi/2.;

  dphix=phi_0*sin(0.5*K[5])*sin(0.5*K[5]);
  dphiy=0.5*phi_0*sin(K[5]);


  e1=K[2]*phi_0;
  e2=K[3]*phi_0;

  omega=phi_0-e1-e2;
  if(phi!=0.) rho=length/phi; else rho=0.;
//  rho=length/phi;
//  rho_0=length/phi_0;

  sin_e1=sin(e1);
  cos_e1=cos(e1);
  sin_e2=sin(e2);
  cos_e2=cos(e2);
  tan_e1=sin_e1/cos_e1;
  tan_e2=sin_e2/cos_e2;

  ds0=length/nbslice;
  dsh=ds0*0.5;

  if(phi!=0) G=(1./rho/rho+K1); else G=0.;
  if(G==0.) {
    a11=1.;
    a12=ds0;
    a11h=1.;
    a12h=dsh;
  }
  else if(G>0) {
    sqr_G=sqrt(G);
    a11=cos(sqr_G*ds0);
    a12=sin(sqr_G*ds0)/sqr_G;
    a11h=cos(sqr_G*dsh);
    a12h=sin(sqr_G*dsh)/sqr_G;
  }
  else {
    sqr_G=sqrt(-G);
    a11=cosh(sqr_G*ds0);
    a12=sinh(sqr_G*ds0)/sqr_G;
    a11h=cosh(sqr_G*dsh);
    a12h=sinh(sqr_G*dsh)/sqr_G;
  }

  if(K1==0.) {
    a33=1.;
    a34=ds0;
    a33h=1.;
    a34h=dsh;
  }
  else if(K1>0.) {
    sqr_K=sqrt(K1);
    a33=cosh(sqr_K*ds0);
    a34=sinh(sqr_K*ds0)/sqr_K;
    a33h=cosh(sqr_K*dsh);
    a34h=sinh(sqr_K*dsh)/sqr_K;
  }
  else {
    sqr_K=sqrt(-K1);
    a33=cos(sqr_K*ds0);
    a34=sin(sqr_K*ds0)/sqr_K;
    a33h=cos(sqr_K*dsh);
    a34h=sin(sqr_K*dsh)/sqr_K;
  }

  if(G!=0. && rho>0.) {
    a16=(1.-a11)/rho/G;
    a51=-a12/rho;
    a52=-a16;
    a56=-(ds0-a12)/rho/rho/G;

    a16h=(1.-a11h)/rho/G;
    a51h=-a12h/rho;
    a52h=-a16h;
    a56h=-(dsh-a12h)/rho/rho/G;
  }
  else {
    a16=0.;a51=0.;a52=0.;a56=0.;
    a16h=0.;a51h=0.;a52h=0.;a56h=0.;
  }

  if(K[14]==3. && rho!=0.) {
    fx=K[13]*K[13]/24./rho;
    fy=K[13]/6./rho/rho;
  }
  else {fx=0.; fy=0.;}
      // cout<< "initialize Bend\n";
}


// Linearized Bending magnet

void LBend::Mapping(map_double& x)
{
  double x0;
  double Krs;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }

  x[0]-=fx*x[5];
  x[2]+=fx*x[3];
  x[4]+=fy*x[1];

  x[3]+=tan_e1/rho*x[0];
  x[4]-=tan_e1/rho*x[1];

  x[0]+=x[1]*x[1]/(2.*rho);
  x[4]-=x[3]*x[1]/rho;

//  Body

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

  Krs=K1/2./rho*ds0;

  for(int i=0;i<nbslice-1;i++) {
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//
    x[2]+=a51*x[0]+a52*x[3]+a56*x[5];
    x0=a11*x[0]+a12*x[3]+a16*x[5];
    x[3]=-a12*G*x[0]+a11*x[3]+a12/rho*x[5];
    x[0]=x0;
    x0=a33*x[1]+a34*x[4];
    x[4]=a34*K1*x[1]+a33*x[4];
    x[1]=x0;
  }      
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

// Body end
  x[0]-=x[1]*x[1]/(2.*rho);
  x[4]+=x[3]*x[1]/rho;

  x[3]+=tan_e2/rho*x[0];
  x[4]-=tan_e2/rho*x[1];

  x[0]+=fx*x[5];
  x[2]-=fx*x[3];
  x[4]+=fy*x[1];

}

void LBend::Mapping(map_da& x)
{
  da x0;
  double Krs;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }

  x[0]-=fx*x[5];
  x[2]+=fx*x[3];
  x[4]+=fy*x[1];

  x[3]+=tan_e1/rho*x[0];
  x[4]-=tan_e1/rho*x[1];

  x[0]+=x[1]*x[1]/(2.*rho);
  x[4]-=x[3]*x[1]/rho;

//  Body

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

  Krs=K1/2./rho*ds0;

  for(int i=0;i<nbslice-1;i++) {
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//
    x[2]+=a51*x[0]+a52*x[3]+a56*x[5];
    x0=a11*x[0]+a12*x[3]+a16*x[5];
    x[3]=-a12*G*x[0]+a11*x[3]+a12/rho*x[5];
    x[0]=x0;
    x0=a33*x[1]+a34*x[4];
    x[4]=a34*K1*x[1]+a33*x[4];
    x[1]=x0;
  }      
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

// Body end

  x[0]-=x[1]*x[1]/(2.*rho);
  x[4]+=x[3]*x[1]/rho;

  x[3]+=tan_e2/rho*x[0];
  x[4]-=tan_e2/rho*x[1];

  x[0]+=fx*x[5];
  x[2]-=fx*x[3];
  x[4]+=fy*x[1];

}

void LBend::Mapping(map_p_da& x)
{
  p_da x0;
  double Krs;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }

  x[0]-=fx*x[5];
  x[2]+=fx*x[3];
  x[4]+=fy*x[1];

  x[3]+=tan_e1/rho*x[0];
  x[4]-=tan_e1/rho*x[1];

  x[0]+=x[1]*x[1]/(2.*rho);
  x[4]-=x[3]*x[1]/rho;

//  Body

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

  Krs=K1/2./rho*ds0;

  for(int i=0;i<nbslice-1;i++) {
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//
    x[2]+=a51*x[0]+a52*x[3]+a56*x[5];
    x0=a11*x[0]+a12*x[3]+a16*x[5];
    x[3]=-a12*G*x[0]+a11*x[3]+a12/rho*x[5];
    x[0]=x0;
    x0=a33*x[1]+a34*x[4];
    x[4]=a34*K1*x[1]+a33*x[4];
    x[1]=x0;
  }      
// Nonlinear kick
    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];

    x[2]-=0.5*ds0*(x[3]*x[3]+x[4]*x[4]);
    x[0]=(x[0]-x[3]*x[5]*ds0)/(1.-ds0/rho*x[3]);
    x[3]-=(x[3]*x[3]+x[4]*x[4])*ds0/rho;
    x[1]+=(x[0]/rho-x[5])*x[4]*ds0;

    x[3]-=Krs*0.5*(3*x[0]*x[0]-x[1]*x[1]);
    x[4]+=Krs*x[0]*x[1];
//

  x[2]+=a51h*x[0]+a52h*x[3]+a56h*x[5];
  x0=a11h*x[0]+a12h*x[3]+a16h*x[5];
  x[3]=-a12h*G*x[0]+a11h*x[3]+a12h/rho*x[5];
  x[0]=x0;
  x0=a33h*x[1]+a34h*x[4];
  x[4]=a34h*K1*x[1]+a33h*x[4];
  x[1]=x0;

// Body end

  x[0]-=x[1]*x[1]/(2.*rho);
  x[4]+=x[3]*x[1]/rho;

  x[3]+=tan_e2/rho*x[0];
  x[4]-=tan_e2/rho*x[1];

  x[0]+=fx*x[5];
  x[2]-=fx*x[3];
  x[4]+=fy*x[1];

}

void LBend::Mapping(pBeam& x)
{
  int i;
  double x0;
  double Krs;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) {
    for(i=0;i<x.np;i++) {
      x.px[i]+=phi_0*x.pz[i]-dphi; x.z[i]-=phi_0*x.x[i];
    }
    return;
  }

  for(i=0;i<x.np;i++) {
    x.x[i]-=fx*x.pz[i];
    x.z[i]+=fx*x.px[i];
    x.py[i]+=fy*x.y[i];

    x.px[i]+=tan_e1/rho*x.x[i];
    x.py[i]-=tan_e1/rho*x.y[i];

    x.x[i]+=x.y[i]*x.y[i]/(2.*rho);
    x.py[i]-=x.px[i]*x.y[i]/rho;

//  Body

    x.z[i]+=a51h*x.x[i]+a52h*x.px[i]+a56h*x.pz[i];
    x0=a11h*x.x[i]+a12h*x.px[i]+a16h*x.pz[i];
    x.px[i]=-a12h*G*x.x[i]+a11h*x.px[i]+a12h/rho*x.pz[i];
    x.x[i]=x0;
    x0=a33h*x.y[i]+a34h*x.py[i];
    x.py[i]=a34h*K1*x.y[i]+a33h*x.py[i];
    x.y[i]=x0;

    Krs=K1/2./rho*ds0;

    for(i=0;i<nbslice-1;i++) {
// Nonlinear kick
      x.px[i]-=Krs*0.5*(3*x.x[i]*x.x[i]-x.y[i]*x.y[i]);
      x.py[i]+=Krs*x.x[i]*x.y[i];

      x.z[i]-=0.5*ds0*(x.px[i]*x.px[i]+x.py[i]*x.py[i]);
      x.x[i]=(x.x[i]-x.px[i]*x.pz[i]*ds0)/(1.-ds0/rho*x.px[i]);
      x.px[i]-=(x.px[i]*x.px[i]+x.py[i]*x.py[i])*ds0/rho;
      x.y[i]+=(x.x[i]/rho-x.pz[i])*x.py[i]*ds0;

      x.px[i]-=Krs*0.5*(3*x.x[i]*x.x[i]-x.y[i]*x.y[i]);
      x.py[i]+=Krs*x.x[i]*x.y[i];
//
      x.z[i]+=a51*x.x[i]+a52*x.px[i]+a56*x.pz[i];
      x0=a11*x.x[i]+a12*x.px[i]+a16*x.pz[i];
      x.px[i]=-a12*G*x.x[i]+a11*x.px[i]+a12/rho*x.pz[i];
      x.x[i]=x0;
      x0=a33*x.y[i]+a34*x.py[i];
      x.py[i]=a34*K1*x.y[i]+a33*x.py[i];
      x.y[i]=x0;
  }
// Nonlinear kick
    x.px[i]-=Krs*0.5*(3*x.x[i]*x.x[i]-x.y[i]*x.y[i]);
    x.py[i]+=Krs*x.x[i]*x.y[i];

    x.z[i]-=0.5*ds0*(x.px[i]*x.px[i]+x.py[i]*x.py[i]);
    x.x[i]=(x.x[i]-x.px[i]*x.pz[i]*ds0)/(1.-ds0/rho*x.px[i]);
    x.px[i]-=(x.px[i]*x.px[i]+x.py[i]*x.py[i])*ds0/rho;
    x.y[i]+=(x.x[i]/rho-x.pz[i])*x.py[i]*ds0;

    x.px[i]-=Krs*0.5*(3*x.x[i]*x.x[i]-x.y[i]*x.y[i]);
    x.py[i]+=Krs*x.x[i]*x.y[i];
//
    x.z[i]+=a51h*x.x[i]+a52h*x.px[i]+a56h*x.pz[i];
    x0=a11h*x.x[i]+a12h*x.px[i]+a16h*x.pz[i];
    x.px[i]=-a12h*G*x.x[i]+a11h*x.px[i]+a12h/rho*x.pz[i];
    x.x[i]=x0;
    x0=a33h*x.y[i]+a34h*x.py[i];
    x.py[i]=a34h*K1*x.y[i]+a33h*x.py[i];
    x.y[i]=x0;

    x.x[i]-=x.y[i]*x.y[i]/(2.*rho);
    x.py[i]+=x.px[i]*x.y[i]/rho;

    x.px[i]+=tan_e2/rho*x.x[i];
    x.py[i]-=tan_e2/rho*x.y[i];

    x.x[i]+=fx*x.pz[i];
    x.z[i]-=fx*x.px[i];
    x.py[i]+=fy*x.y[i];
  }
}

