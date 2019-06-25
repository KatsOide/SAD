#include <cmath>
using std::sqrt;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <element.h>


// Drift

void DriftMapping(double length,map_double& x)
{
   double pz,E,pxy2;

   pxy2=x[3]*x[3]+x[4]*x[4];
   E=1.+x[5];	// x0=(1+\delta)
   pz=sqrt(E*E-pxy2);	// pz=\sqrt{(1+\delta)^2-px^2-py^2}

   double pz_i=1./pz;
   
   x[0]+=x[3]*pz_i*length;	// x=x0+px0*length
   x[1]+=x[4]*pz_i*length;	// y=y0+py0*length
   x[2]-=pxy2/(pz*(pz+E))*length;
}

void DriftMapping(double length,pBeam& x)
{
  double pz,E,pxy2;

  for(int i=0;i<x.np;i++) {
    pxy2=x.px[i]*x.px[i]+x.py[i]*x.py[i];
    E=1.+x.pz[i];	        // x0=(1+\delta)
    pz=sqrt(E*E-pxy2);	// pz=\sqrt{(1+\delta)^2-px^2-py^2}

    double pz_i=1./pz;
   
    x.x[i]+=x.px[i]*pz_i*length;	// x=x0+px0*length
    x.y[i]+=x.py[i]*pz_i*length;	// y=y0+py0*length
    x.z[i]-=pxy2/(pz*(pz+E))*length;
  }
}

// Drift

void DriftMapping(double length,map_da& x)
{
   da pz,E,pxy2;
   pxy2=x[3]*x[3]+x[4]*x[4];
   E=1.+x[5];	// x0=(1+\delta)
   pz=sqrt(E*E-pxy2);	// pz=\sqrt{(1+\delta)^2-px^2-py^2}

   da pz_i=1./pz;
   
   x[0]+=x[3]*pz_i*length;	// x=x0+px0*length
   x[1]+=x[4]*pz_i*length;	// y=y0+py0*length
   x[2]-=pxy2/(pz*(pz+E))*length;
}

// Drift

void DriftMapping(double length,map_p_da& x)
{
  p_da pz,E,pxy2;
  pxy2=x[3]*x[3]+x[4]*x[4];
  E=1.+x[5];	// x0=(1+\delta)
  pz=sqrt(E*E-pxy2);	// pz=\sqrt{(1+\delta)^2-px^2-py^2}

  p_da pz_i=1./pz;

  x[0]+=x[3]*pz_i*length;	// x=x0+px0*length
  x[1]+=x[4]*pz_i*length;	// y=y0+py0*length
  x[2]-=pxy2/(pz*(pz+E))*length;
}
