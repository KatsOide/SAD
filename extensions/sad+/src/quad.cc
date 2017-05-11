#include <cmath>
using std::fabs;
using std::sqrt;
using std::sin;
using std::cos;
using std::sinh;
using std::cosh;
using std::exp;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <quad.h>


// Quadrupole magnet

Quad::Quad(char* s) : Element(s) 
{
   N_div=(int)get_parm(s,"Ndiv");
   K1=get_parm(s,"K1");
   ent_edge=(int)get_parm(s,"ent_edge");
   exit_edge=(int)get_parm(s,"exit_edge");

   ds0=length/(double)N_div;
   rk=fabs(K1/length);
   sqrtk=sqrt(rk);
   if(K1!=0.) {
   a11=cos(sqrtk*ds0);
   a12=sin(sqrtk*ds0)/sqrtk;
   b11=cosh(sqrtk*ds0);
   b12=sinh(sqrtk*ds0)/sqrtk;
   a11h=cos(sqrtk*ds0*.5);
   a12h=sin(sqrtk*ds0*.5)/sqrtk;
   b11h=cosh(sqrtk*ds0*.5);
   b12h=sinh(sqrtk*ds0*.5)/sqrtk;}
   // cout << "initialize Quad\n";
}

Quad::Quad(double* K) : Element(K) 
{
  N_div=3;
  K1=K[1];
  ent_edge=0;
  exit_edge=0;

  ds0=length/(double)N_div;
  if(length>0) {
    rk=fabs(K1/length);
    f1=-K1/length*K[8]*fabs(K[8])/24.;
    f2=K1/length*K[9];
  }
  else {
    rk=0.;
    f1=0.;
    f2=0.;
  }
  sqrtk=sqrt(rk);
  if(K1!=0.) {
    a11=cos(sqrtk*ds0);
    a12=sin(sqrtk*ds0)/sqrtk;
    b11=cosh(sqrtk*ds0);
    b12=sinh(sqrtk*ds0)/sqrtk;
    a11h=cos(sqrtk*ds0*.5);
    a12h=sin(sqrtk*ds0*.5)/sqrtk;
    b11h=cosh(sqrtk*ds0*.5);
    b12h=sinh(sqrtk*ds0*.5)/sqrtk;
  }
  else {
    a11=1.;    a12=0.;    b11=1.;    b12=0.;
    a11h=1.;    a12h=0.;    b11h=1.;    b12h=0.;
  }
  // cout << "initialize Quad\n";
}


void Quad::print(void)
{
  Element::print(); 
  cout<< ", K1="<< K1 << ", Ndiv=" << N_div << ", ent_edge=" << ent_edge
      << ", exit edge=" << exit_edge << ")\n";
}

void Quad::Mapping(map_double& x)
{
  double x_f,y_f,E,pxy2,dl,pz;
  double a21,b21,a21h,b21h,ef1;

  if(K1==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) {
    x[3]-=K1*x[0];
    x[4]+=K1*x[1];
    return;
  }
  if(ent_edge) { // entrance edge effect
    double a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;
    E=1./(1.+x[5]);
    Kdelta=0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;

    ef1=exp(f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;
  }

  for(int i=0;i<=N_div;i++){
    if(N_div==0 || (i>0 && i<N_div)) { 
      a21=-a12*rk;
      b21=b12*rk; 

      if(K1>0.) {
	x_f=a11*x[0]+a12*x[3];
	x[3]=a21*x[0]+a11*x[3];
	x[0]=x_f;
	y_f=b11*x[1]+b12*x[4];
	x[4]=b21*x[1]+b11*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11*x[0]+b12*x[3];
	x[3]=b21*x[0]+b11*x[3];
	x[0]=x_f;
	y_f=a11*x[1]+a12*x[4];
	x[4]=a21*x[1]+a11*x[4];
	x[1]=y_f;
      }
    }
    else {
      a21h=-a12h*rk; 
      b21h=b12h*rk; 

      if(K1>0.) {
	x_f=a11h*x[0]+a12h*x[3];
	x[3]=a21h*x[0]+a11h*x[3];
	x[0]=x_f;
	y_f=b11h*x[1]+b12h*x[4];
	x[4]=b21h*x[1]+b11h*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11h*x[0]+b12h*x[3];
	x[3]=b21h*x[0]+b11h*x[3];
	x[0]=x_f;
	y_f=a11h*x[1]+a12h*x[4];
	x[4]=a21h*x[1]+a11h*x[4];
	x[1]=y_f;
      }
    }

   // nonlinear kick
    if(i!=N_div){
      E=1.+x[5];
      pxy2=x[3]*x[3]+x[4]*x[4];
      pz=sqrt(E*E-pxy2);
      dl=(pxy2-x[5]*(2.+x[5]))/(pz*(1.+pz))*ds0;
      x[0]+=x[3]*dl;
      x[1]+=x[4]*dl;
      x[2]-=pxy2/(pz*(pz+E))*ds0;
    }
  }

  if(exit_edge) { // exit edge effect
    double a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;

    ef1=exp(-f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;

    E=1./(1.+x[5]);
    Kdelta=-0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;
  }
}
	 

// Quadrupole magnet

void Quad::Mapping(map_da& x)
{
  da x_f,y_f,E,pxy2,dl,pz;
  double a21,b21,a21h,b21h,ef1;

  if(K1==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) {
    x[3]-=K1*x[0];
    x[4]+=K1*x[1];
    return;
  }
   
  if(ent_edge) { // entrance edge effect
    da a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;
    E=1./(1.+x[5]);
    Kdelta=0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;

    ef1=exp(f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;

  }
  //cout << "ent edge kick\n" << x0;
  //is_symplectic(x0);
  //   cout << x; cout.flush();
  for(int i=0;i<=N_div;i++){
    if(N_div==0 || (i>0 && i<N_div)) { 
      a21=-a12*rk;
      b21=b12*rk; 

      if(K1>0.) {
	x_f=a11*x[0]+a12*x[3];
	x[3]=a21*x[0]+a11*x[3];
	x[0]=x_f;
	y_f=b11*x[1]+b12*x[4];
	x[4]=b21*x[1]+b11*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11*x[0]+b12*x[3];
	x[3]=b21*x[0]+b11*x[3];
	x[0]=x_f;
	y_f=a11*x[1]+a12*x[4];
	x[4]=a21*x[1]+a11*x[4];
	x[1]=y_f;
      }
    }
    else {
      a21h=-a12h*rk; 
      b21h=b12h*rk; 

      if(K1>0.) {
	x_f=a11h*x[0]+a12h*x[3];
	x[3]=a21h*x[0]+a11h*x[3];
	x[0]=x_f;
	y_f=b11h*x[1]+b12h*x[4];
	x[4]=b21h*x[1]+b11h*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11h*x[0]+b12h*x[3];
	x[3]=b21h*x[0]+b11h*x[3];
	x[0]=x_f;
	y_f=a11h*x[1]+a12h*x[4];
	x[4]=a21h*x[1]+a11h*x[4];
	x[1]=y_f;
      }
    }

   // nonlinear kick
    if(i!=N_div){
      E=1.+x[5];
      pxy2=x[3]*x[3]+x[4]*x[4];
      pz=sqrt(E*E-pxy2);
      dl=(pxy2-x[5]*(2.+x[5]))/(pz*(1.+pz))*ds0;
      x[0]+=x[3]*dl;
      x[1]+=x[4]*dl;
      x[2]-=pxy2/(pz*(pz+E))*ds0;
    }
  }
 
  if(exit_edge) { // exit edge effect
    da a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;

    ef1=exp(-f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;

    E=1./(1.+x[5]);
    Kdelta=-0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;
  }
}
	 
	 

// Quadrupole magnet

void Quad::Mapping(map_p_da& x)
{
  p_da x_f,y_f,E,pxy2,dl,pz;
  double a21,b21,a21h,b21h,ef1;
 
  if(K1==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) {
    x[3]-=K1*x[0];
    x[4]+=K1*x[1];
    return;
  }

  if(ent_edge) { // entrance edge effect
    p_da a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;
    E=1./(1.+x[5]);
    Kdelta=0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;

    ef1=exp(f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;

  }

  for(int i=0;i<=N_div;i++){
    if(N_div==0 || (i>0 && i<N_div)) { 
      a21=-a12*rk;
      b21=b12*rk; 

      if(K1>0.) {
	x_f=a11*x[0]+a12*x[3];
	x[3]=a21*x[0]+a11*x[3];
	x[0]=x_f;
	y_f=b11*x[1]+b12*x[4];
	x[4]=b21*x[1]+b11*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11*x[0]+b12*x[3];
	x[3]=b21*x[0]+b11*x[3];
	x[0]=x_f;
	y_f=a11*x[1]+a12*x[4];
	x[4]=a21*x[1]+a11*x[4];
	x[1]=y_f;
      }
    }
    else {
      a21h=-a12h*rk; 
      b21h=b12h*rk; 

      if(K1>0.) {
	x_f=a11h*x[0]+a12h*x[3];
	x[3]=a21h*x[0]+a11h*x[3];
	x[0]=x_f;
	y_f=b11h*x[1]+b12h*x[4];
	x[4]=b21h*x[1]+b11h*x[4];
	x[1]=y_f;
      }
      else{
	x_f=b11h*x[0]+b12h*x[3];
	x[3]=b21h*x[0]+b11h*x[3];
	x[0]=x_f;
	y_f=a11h*x[1]+a12h*x[4];
	x[4]=a21h*x[1]+a11h*x[4];
	x[1]=y_f;
      }
    }

   // nonlinear kick
    if(i!=N_div){
      E=1.+x[5];
      pxy2=x[3]*x[3]+x[4]*x[4];
      pz=sqrt(E*E-pxy2);
      dl=(pxy2-x[5]*(2.+x[5]))/(pz*(1.+pz))*ds0;
      x[0]+=x[3]*dl;
      x[1]+=x[4]*dl;
      x[2]-=pxy2/(pz*(pz+E))*ds0;
    }
  }
  if(exit_edge) { // exit edge effect
    p_da a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;

    ef1=exp(-f1);
    x[0]=ef1*x[0]+f2*x[3];
    x[1]=x[1]/ef1-f2*x[4];
    x[3]=x[3]/ef1;
    x[4]=x[4]*ef1;

    E=1./(1.+x[5]);
    Kdelta=-0.25*K1/length*E;
    a=Kdelta*x[0]*x[0];
    b=Kdelta*x[1]*x[1];
    ab=Kdelta*(x[0]-x[1])*(x[0]+x[1]);
    t=ab*ab/6.;
    xx=1.+a+b+ab*(5.*a-b)/6.;
    yy=1.-a-b+ab*(a-5.*b)/6.;
    h=2.*Kdelta*x[0]*x[1]*(1.-ab/3.);
    f=1./(xx*yy+h*h);
    tmpx=(yy*x[3]+h*x[4])*f;
    x[4]=(xx*x[4]-h*x[3])*f;
    x[3]=tmpx;
    tmpx=x[0]*(a/3.+b+t);
    tmpy=x[1]*(a+b/3.-t);
    x[2]+=E*((tmpy-t*x[1])*x[4]-(tmpx+t*x[0])*x[3]);
    x[0]+=tmpx;
    x[1]-=tmpy;
  }
}
	 

void Quad::Mapping(pBeam& x)
{
  int i;
  double x_f,y_f,E,pxy2,dl,pz;
  double a21,b21,a21h,b21h,ef1;

  if(K1==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) {
    for(i=0;i<x.np;i++) {
      x.px[i]-=K1*x.x[i];
      x.py[i]+=K1*x.y[i];
    }
    return;
  }

  if(ent_edge) { // entrance edge effect
    for(i=0;i<x.np;i++) {
      double a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;

      ef1=exp(f1);
      x.x[i]=ef1*x.x[i]+f2*x.px[i];
      x.y[i]=x.y[i]/ef1-f2*x.py[i];
      x.px[i]=x.px[i]/ef1;
      x.py[i]=x.py[i]*ef1;

      E=1./(1.+x.pz[i]);
      Kdelta=0.25*K1/length*E;
      a=Kdelta*x.x[i]*x.x[i];
      b=Kdelta*x.y[i]*x.y[i];
      ab=Kdelta*(x.x[i]-x.y[i])*(x.x[i]+x.y[i]);
      t=ab*ab/6.;
      xx=1.+a+b+ab*(5.*a-b)/6.;
      yy=1.-a-b+ab*(a-5.*b)/6.;
      h=2.*Kdelta*x.x[i]*x.y[i]*(1.-ab/3.);
      f=1./(xx*yy+h*h);
      tmpx=(yy*x.px[i]+h*x.py[i])*f;
      x.py[i]=(xx*x.py[i]-h*x.px[i])*f;
      x.px[i]=tmpx;
      tmpx=x.x[i]*(a/3.+b+t);
      tmpy=x.y[i]*(a+b/3.-t);
      x.z[i]+=E*((tmpy-t*x.y[i])*x.py[i]-(tmpx+t*x.x[i])*x.px[i]);
      x.x[i]+=tmpx;
      x.y[i]-=tmpy;
    }
  }

  for(int j=0;j<=N_div;j++){
    if(N_div==0 || (j>0 && j<N_div)) { 
      a21=-a12*rk;
      b21=b12*rk; 

      if(K1>0.) {
	for(i=0;i<x.np;i++) {
	  x_f=a11*x.x[i]+a12*x.px[i];
	  x.px[i]=a21*x.x[i]+a11*x.px[i];
	  x.x[i]=x_f;
	  y_f=b11*x.y[i]+b12*x.py[i];
	  x.py[i]=b21*x.y[i]+b11*x.py[i];
	  x.y[i]=y_f;
	}
      }
      else{
	for(i=0;i<x.np;i++) {
	  x_f=b11*x.x[i]+b12*x.px[i];
	  x.px[i]=b21*x.x[i]+b11*x.px[i];
	  x.x[i]=x_f;
	  y_f=a11*x.y[i]+a12*x.py[i];
	  x.py[i]=a21*x.y[i]+a11*x.py[i];
	  x.y[i]=y_f;
	}
      }
    }
    else {
      a21h=-a12h*rk; 
      b21h=b12h*rk; 

      if(K1>0.) {
	for(i=0;i<x.np;i++) {
	  x_f=a11h*x.x[i]+a12h*x.px[i];
	  x.px[i]=a21h*x.x[i]+a11h*x.px[i];
	  x.x[i]=x_f;
	  y_f=b11h*x.y[i]+b12h*x.py[i];
	  x.py[i]=b21h*x.y[i]+b11h*x.py[i];
	  x.y[i]=y_f;
	}
      }
      else {
	for(i=0;i<x.np;i++) {
	  x_f=b11h*x.x[i]+b12h*x.px[i];
	  x.px[i]=b21h*x.x[i]+b11h*x.px[i];
	  x.x[i]=x_f;
	  y_f=a11h*x.y[i]+a12h*x.py[i];
	  x.py[i]=a21h*x.y[i]+a11h*x.py[i];
	  x.y[i]=y_f;
	}
      }
    }

   // nonlinear kick
    if(j!=N_div){
      for(i=0;i<x.np;i++) {
	E=1.+x.pz[i];
	pxy2=x.px[i]*x.px[i]+x.py[i]*x.py[i];
	pz=sqrt(E*E-pxy2);
	dl=(pxy2-x.pz[i]*(2.+x.pz[i]))/(pz*(1.+pz))*ds0;
	x.x[i]+=x.px[i]*dl;
	x.y[i]+=x.py[i]*dl;
	x.z[i]-=pxy2/(pz*(pz+E))*ds0;
      }
    }
  }

  if(exit_edge) { // exit edge effect
    for(i=0;i<x.np;i++) {
      double a,b,ab,t,h,Kdelta,f,xx,yy,tmpx,tmpy;
      E=1./(1.+x.pz[i]);
      Kdelta=-0.25*K1/length*E;
      a=Kdelta*x.x[i]*x.x[i];
      b=Kdelta*x.y[i]*x.y[i];
      ab=Kdelta*(x.x[i]-x.y[i])*(x.x[i]+x.y[i]);
      t=ab*ab/6.;
      xx=1.+a+b+ab*(5.*a-b)/6.;
      yy=1.-a-b+ab*(a-5.*b)/6.;
      h=2.*Kdelta*x.x[i]*x.y[i]*(1.-ab/3.);
      f=1./(xx*yy+h*h);
      tmpx=(yy*x.px[i]+h*x.py[i])*f;
      x.py[i]=(xx*x.py[i]-h*x.px[i])*f;
      x.px[i]=tmpx;
      tmpx=x.x[i]*(a/3.+b+t);
      tmpy=x.y[i]*(a+b/3.-t);
      x.z[i]+=E*((tmpy-t*x.y[i])*x.py[i]-(tmpx+t*x.x[i])*x.px[i]);
      x.x[i]+=tmpx;
      x.y[i]-=tmpy;

      ef1=exp(-f1);
      x.x[i]=ef1*x.x[i]+f2*x.px[i];
      x.y[i]=x.y[i]/ef1-f2*x.py[i];
      x.px[i]=x.px[i]/ef1;
      x.py[i]=x.py[i]*ef1;

    }
  }
}
	 
