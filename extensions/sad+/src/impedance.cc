#include <iostream>
using std::cout;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <phys_const.h>
#include <impedance.h>
void spline(double x[],double y[],int n,double yp1,double ypn,double y2[]);
void splint(double xa[],double ya[],double y2a[],int n,double x,double *y);


Impedance::Impedance(char* s) : Element(s)
{
  int i;
  double b1,c1;

  niz=(int)get_parm(s,"NzSlice");
  Zw=get_parm(s,"ZWidth");
  if(niz%2) niz++;
  dZ=Zw/niz;

  //kfac is set at beam_set (impedance.h)
  b1=get_parm(s,"W1x1");
  c1=get_parm(s,"W1y1");

  Z=new double [niz];
  W0P=new double [niz*3];
  W1x=W0P+niz;
  W1y=W0P+2*niz;

  F0=new double [niz*6];
  F1x=F0+niz;
  F1y=F0+2*niz;
  G0=F0+3*niz;   //second order derivative of F0
  G1x=F0+4*niz;
  G1y=F0+5*niz;

  rho0=new double [niz*3];
  rho1x=rho0+niz;
  rho1y=rho0+2*niz;

  for(i=0;i<6*niz;i++) F0[i]=0.;
  for(i=0;i<niz;i++) W0P[i]=0.;
  for(i=0;i<niz;i++) W1x[i]=b1*i*dZ;
  for(i=0;i<niz;i++) W1y[i]=c1*i*dZ;
  for(i=0;i<3*niz;i++) rho0[i]=0.;
  for(i=0;i<niz;i++) Z[i]=Zw*2.*double(i-niz/2+0.5)/double(niz);
}

Impedance::Impedance(double* K) : Element(K)
{
  int i;
  double b1,c1;

  niz=(int)K[1];
  Zw=K[2];
  if(niz%2) niz++;
  dZ=2.*Zw/niz;

  //kfac is set at beam_set (impedance.h)
  b1=K[3];
  c1=K[4];

  Z=new double [niz];
  W0P=new double [niz*3];
  W1x=W0P+niz;
  W1y=W0P+2*niz;

  F0=new double [niz*6];
  F1x=F0+niz;
  F1y=F0+2*niz;
  G0=F0+3*niz;   //second order derivative of F0
  G1x=F0+4*niz;
  G1y=F0+5*niz;

  rho0=new double [niz*3];
  rho1x=rho0+niz;
  rho1y=rho0+2*niz;

  for(i=0;i<9*niz;i++) F0[i]=0.;
  for(i=0;i<niz;i++) W0P[i]=0.;
  for(i=0;i<niz;i++) W1x[i]=b1*i*dZ;
  for(i=0;i<niz;i++) W1y[i]=c1*i*dZ;
  for(i=0;i<3*niz;i++) rho0[i]=0.;
  for(i=0;i<niz;i++) Z[i]=Zw*2.*double(i-niz/2+0.5)/double(niz);
  print();
}

Impedance::~Impedance(void)
{
  delete rho0;
  delete F0;
  delete W0P;
  delete Z;
}

void Impedance::print(void)
{
  cout << "Ring wake founction\n";
  cout << "Nzslice=" << niz << "   Zrange=" << -Zw << " to " << Zw
       << "    dZ=" << dZ << "   kfac=" << kfac << "\n";
  cout << "Wx=" << W1x[1]/dZ << " z    Wy=" << W1y[1]/dZ << " z\n";
  cout.flush();
}

void Impedance::Mapping(map_double& x)
{
}

void Impedance::Mapping(map_da& x)
{
}

void Impedance::Mapping(map_p_da& x)
{
}

void Impedance::Mapping(pBeam& x)
{
  int i,j,iz,niz2;
  
// Calculate dipole momentum
  for(i=0;i<niz;i++) {
    rho0[i]=0.;
    rho1x[i]=0.;
    rho1y[i]=0.;
    F0[i]=0.;
    F1x[i]=0.;
    F1y[i]=0.;
  }

  niz2=niz/2;
  for(i=0;i<x.np;i++) {
    iz=(int)(x.z[i]/dZ)+niz2;
    if(iz<0) iz=0;
    if(iz>=niz) iz=niz-1;
    rho0[iz]+=1.;
    rho1x[iz]+=x.x[i];
    rho1y[iz]+=x.y[i];
  }

  for(i=0;i<niz;i++) {
    rho0[i]/=x.np;
    rho1x[i]/=x.np;
    rho1y[i]/=x.np;
  }
  

// Convolution with Wake force

  for(i=0;i<niz;i++) {
    for(j=i;j<niz;j++) {
      F0[i]+=rho0[j]*W0P[j-i];
      F1x[i]+=rho1x[j]*W1x[j-i];
      F1y[i]+=rho1y[j]*W1y[j-i];
    }
  }
  //  printrho(1);

// spline function
  spline(Z,F1x,niz,0.,0.,G1x);
  spline(Z,F1y,niz,0.,0.,G1y);
// Momentum kick

  for(i=0;i<x.np;i++) {
    splint(Z,F1x,G1x,niz,x.z[i],&K1x); // use spline fit
    splint(Z,F1y,G1y,niz,x.z[i],&K1y); 
    //cout << x.z[i] << " K " << kfac << " " 
    //      << K1y << " " << kfac*K1y << "\n";
    x.px[i]+=kfac*K1x;
    x.py[i]+=kfac*K1y;
  }
  //cout << "\n rhox & Fx " << rho1x[niz2] << " " << kfac*F1x[niz2] 
  //     << " rhoy & Fy " << rho1y[niz2] << " " << kfac*F1y[niz2] << "\n";
}

void Impedance::printrho(int id)
{
  cout << "\nrho\n";
  for(int i=0;i<niz;i++) {
    if(id==0) cout << Z[i] << " " << rho0[i] << " " << F0[i] << "\n";
    if(id==1) cout << Z[i] << " " << rho0[i] << " " 
		   << rho1x[i] << " " << F1x[i] << "\n";
    if(id==2) cout << Z[i] << " " << rho0[i] << " "
		   << rho1y[i] << " " << F1y[i] << "\n";
  }
}


//******** Chromaticity *******************************************

Chromaticity::Chromaticity(char* s) : Element(s)
{
  Kxp=get_parm(s,"Kxp");
  Kyp=get_parm(s,"Kyp");
}

Chromaticity::Chromaticity(double* K) : Element(K)
{
  Kxp=(int)K[1];
  Kyp=(int)K[2];
  print();
}


void Chromaticity::print(void)
{
  cout << "Artificial Chromaticity\n";
  cout << "Kxdelata=" << Kxp << "   Kydelta=" << Kyp  << "\n";
  cout.flush();
}

void Chromaticity::Mapping(map_double& x)
{
}

void Chromaticity::Mapping(map_da& x)
{
}

void Chromaticity::Mapping(map_p_da& x)
{
}

void Chromaticity::Mapping(pBeam& x)
{
  int i;
  
  for(i=0;i<x.np;i++) {
    x.px[i]-=Kxp*x.pz[i]*x.x[i];
    x.py[i]-=Kyp*x.pz[i]*x.y[i];
    x.z[i]+=0.5*(Kxp*x.x[i]*x.x[i]+Kyp*x.y[i]*x.y[i]);
  }
}









