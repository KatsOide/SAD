#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
using std::exit;

void spline(double x[],double y[],int n,double yp1,double ypn,double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;
  u=new double[n];
  if(yp1>.99e33) {
    y2[0]=0.;
    u[0]=0.;}
  else {
    y2[0]=-0.5;
    u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for(i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.;
    y2[i]=(sig-1.0)/p;
    u[i]=(6.*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))
	  /(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if(ypn>.99e30) {
    qn=0.;
    un=0.; }
  else {
    qn=0.5;
    un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.);
  for(k=n-2;k>=0;k--) {
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
  delete u;
}


void splint(double xa[],double ya[],double y2a[],int n,double x,double *y)
{
  int klo,khi,k;
  double h,b,a;
  klo=0;
  khi=n-1;
  while(khi-klo > 1) {
    k=(khi+klo) >>1;
    if(xa[k]>x)  khi=k;
    else  klo=k;
  }

  h=xa[khi]-xa[klo];
  if(h==0.) cout << "Bad XA input." << endl;
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+
    ((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;
}


void bcucof(double y[],double y1[],double y2[],double y12[],double d1,
	    double d2,double *c)
{
  static int wt[16][16]=
    { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
      { 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
      { 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
      { 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
      {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
      { 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
      {-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
      { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
      {-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
      { 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1} };
  int i,j;
  double xx,d1d2,x[16];

  d1d2=d1*d2;
  for(i=0;i<4;i++) {
    x[i]=y[i];
    x[i+4]=y1[i]*d1;
    x[i+8]=y2[i]*d2;
    x[i+12]=y12[i]*d1d2;
  }
  for(i=0;i<16;i++) {
    xx=0.;
    for(j=0;j<16;j++) xx+=wt[i][j]*x[j];
    c[i]=xx;
  }
}


void bcuint(double y[],double y1[],double y2[],double y12[],double x1l,
	    double x1u,double x2l,double x2u,double x1,double x2,
	    double *ansy,double *ansy1,double *ansy2)
{
  double c[16];
  int i,ii;
  double t,u,d1,d2;

  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if(x1u==x1l || x2u==x2l) {cout << "bad input" << endl; exit(1);}
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  *ansy=0.;
  *ansy2=0.;
  *ansy1=0.;
  for(i=3;i>=0;i--) {
    ii=i*4;
    *ansy=t*(*ansy)+((c[ii+3]*u+c[ii+2])*u+c[ii+1])*u+c[ii+0];
    *ansy2=t*(*ansy2)+(3.*c[ii+3]*u+2.*c[ii+2])*u+c[ii+1];
    *ansy1=u*(*ansy1)+(3.*c[12+i]*t+2.*c[8+i])*t+c[4+i];
  }
  *ansy1/=d1;
  *ansy2/=d2;
}


void locate(double xx[],unsigned long n,double x,unsigned long *j)
{
  unsigned long ju,jm,jl;
  int ascnd;
  jl=0;
  ju=n+1;
  ascnd=(xx[n-1] > xx[0]);
  while(ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if((x > xx[jm]) == ascnd ) jl=jm;
    else ju=jm;
  }
  *j=jl;
}

