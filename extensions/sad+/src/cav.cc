#include <cmath>
using std::sin;
using std::cos;

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <element.h>
#include <cavity.h>



// Cavity

Cavity::Cavity(char* s) : Element(s)
{
     harm=(int)get_parm(s,"harm");
     freq=get_parm(s,"freq");
     phi=get_parm(s,"phi");
     voltage=get_parm(s,"volt");

     w=pi2*freq/C;
     sinp=sin(phi);
     cosp=cos(phi);

}

Cavity::Cavity(double* K)  : Element(K) 
{
  voltage=K[1];
  harm=(int)K[2];
  freq=K[4];

  phi=K[3];
  w=pi2*freq/C;
  sinp=sin(phi);
  cosp=cos(phi);
}

void Cavity::Mapping(map_double& x)
{
  if(VE==0.) return;

  double sz=sin(w*0.5*x[2]);
  x[5]+=VE*(sin(w*x[2])*cosp-2.*sinp*sz*sz);
}

void Cavity::Mapping(map_da& x)
{
  if(VE==0.) return;

  da sz=sin(w*0.5*x[2]);
  x[5]+=VE*(sin(w*x[2])*cosp-2.*sinp*sz*sz);
}

void Cavity::Mapping(map_p_da& x)
{
  if(VE==0.) return;

  p_da sz=sin(w*0.5*x[2]);
  x[5]+=VE*(sin(w*x[2])*cosp-2.*sinp*sz*sz);
}


void Cavity::Mapping(pBeam& x)
{
  if(VE==0.) return;

  for(int i=0;i<x.np;i++) {
    double sz=sin(w*0.5*x.z[i]);
    x.pz[i]+=VE*(sin(w*x.z[i])*cosp-2.*sinp*sz*sz);
  }
}

