#include <cmath>

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <track.h>
#include <c_da.h>
#include <c_p_da.h>
#include <element.h>
#include <thin.h>

	 
// multipole

void Thin::Mapping(map_double& x)
{
  int Nkai=n_MP;
  if(K_n==0.) return;

  for(int i=1;i<n_MP;i++) Nkai=Nkai*i;
  double K_i=K_n/((double)Nkai);
  Complex z(x[0],-x[1]);
  z=pow(z,n_MP);
  x[3]-=K_i*real(z);
  x[4]-=K_i*imag(z);
}

void Thin::Mapping(map_da& x)
{
  int Nkai=n_MP;
  if(K_n==0.) return;

  for(int i=1;i<n_MP;i++) Nkai=Nkai*i;
  double K_i=K_n/((double)Nkai);
  c_da z(x[0],-x[1]);
  z=pow(z,n_MP);
  x[3]-=K_i*real(z);
  x[4]-=K_i*imag(z);
}

void Thin::Mapping(map_p_da& x)
{
  int Nkai=n_MP;
  if(K_n==0.) return;

  for(int i=1;i<n_MP;i++) Nkai=Nkai*i;
  double K_i=K_n/((double)Nkai);
  c_p_da z(x[0],-x[1]);
  z=pow(z,n_MP);
  x[3]-=K_i*real(z);
  x[4]-=K_i*imag(z);
}

void Thin::Mapping(pBeam& x)
{
  int i;
  int Nkai=n_MP;
  if(K_n==0.) return;

  for(i=1;i<n_MP;i++) Nkai=Nkai*i;
  double K_i=K_n/((double)Nkai);
  for(i=0;i<x.np;i++) {
    Complex z0(x.x[i],-x.y[i]);
    z0=pow(z0,n_MP);
    x.px[i]-=K_i*real(z0);
    x.py[i]-=K_i*imag(z0);
  }
}

