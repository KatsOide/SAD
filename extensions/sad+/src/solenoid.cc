#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <solenoid.h>

// Solenoid magnet

Solenoid::Solenoid(char* s) : Element(s) 
{
   Ks=get_parm(s,"Ks");
   rk=Ks/length;
   //   Ks=Bs/(0.3*enerygy);
}

Solenoid::Solenoid(double* K) : Element(K) 
{
  Ks=K[1];
}

void Solenoid::print(void)
{
   Element::print(); 
   cout<< ", Ks="<< Ks << ")\n";
}

void Solenoid::Mapping(map_double& x)
{
  double x0,x3,x1;
  double pz,E,pxy2,kpz,vx,vy,cosks,sinks;
  double t11,t21,t31,t41,t12,t32;

  E=1.+x[5];
  vx=x[3]+0.5*rk*x[1];
  vy=x[4]-0.5*rk*x[0];
  pxy2=vx*vx+vy*vy;
  pz=sqrt(E*E-pxy2);
  kpz=Ks/pz;
  cosks=cos(kpz);
  sinks=sin(kpz);
  t11=0.5*(1.+cosks);
  t21=-0.25*rk*sinks;
  t31=-0.5*sinks;
  t41=rk*0.25*(1.-cosks);
  t12=sinks/rk;
  t32=(1.-cosks)/(-rk);

  x0  =  t11*x[0]+t12*x[3]-t31*x[1]-t32*x[4];
  x3  =  t21*x[0]+t11*x[3]-t41*x[1]-t31*x[4];
  x1  =  t31*x[0]+t32*x[3]+t11*x[1]+t12*x[4];
  x[4]=  t41*x[0]+t31*x[3]+t21*x[1]+t11*x[4];

  x[0]=x0; x[1]=x1; x[3]=x3;
  x[2]-=pxy2/(pz*(pz+E))*length;

}

void Solenoid::Mapping(map_da& x)
{
  da x0,x3,x1;
  da pz,E,pxy2,kpz,vx,vy,cosks,sinks;
  da t11,t21,t31,t41,t12,t32;

  E=1.+x[5];
  vx=x[3]+0.5*rk*x[1];
  vy=x[4]-0.5*rk*x[0];
  pxy2=vx*vx+vy*vy;
  pz=sqrt(E*E-pxy2);
  kpz=Ks/pz;
  cosks=cos(kpz);
  sinks=sin(kpz);
  t11=0.5*(1.+cosks);
  t21=-0.25*rk*sinks;
  t31=-0.5*sinks;
  t41=rk*0.25*(1.-cosks);
  t12=sinks/rk;
  t32=(1.-cosks)/(-rk);
   
  x0=  t11*x[0]+t12*x[3]-t31*x[1]-t32*x[4];
  x3=  t21*x[0]+t11*x[3]-t41*x[1]-t31*x[4];
  x1=  t31*x[0]+t32*x[3]+t11*x[1]+t12*x[4];
  x[4]=t41*x[0]+t31*x[3]+t21*x[1]+t11*x[4];

  x[0]=x0; x[1]=x1; x[3]=x3;
  x[2]-=pxy2/(pz*(pz+E))*length;

}

void Solenoid::Mapping(map_p_da& x)
{
  p_da x0,x3,x1;
  p_da pz,E,pxy2,kpz,vx,vy,cosks,sinks;
  p_da t11,t21,t31,t41,t12,t32;

  E=1.+x[5];
  vx=x[3]+0.5*rk*x[1];
  vy=x[4]-0.5*rk*x[0];
  pxy2=vx*vx+vy*vy;
  pz=sqrt(E*E-pxy2);
  kpz=Ks/pz;
  cosks=cos(kpz);
  sinks=sin(kpz);
  t11=0.5*(1.+cosks);
  t21=-0.25*rk*sinks;
  t31=-0.5*sinks;
  t41=rk*0.25*(1.-cosks);
  t12=sinks/rk;
  t32=(1.-cosks)/(-rk);
   
  x0=  t11*x[0]+t12*x[3]-t31*x[1]-t32*x[4];
  x3=  t21*x[0]+t11*x[3]-t41*x[1]-t31*x[4];
  x1=  t31*x[0]+t32*x[3]+t11*x[1]+t12*x[4];
  x[4]=t41*x[0]+t31*x[3]+t21*x[1]+t11*x[4];

  x[0]=x0; x[1]=x1; x[3]=x3;
  x[2]-=pxy2/(pz*(pz+E))*length;

}

void Solenoid::Mapping(pBeam& x)
{
}



