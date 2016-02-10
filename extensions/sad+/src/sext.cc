#include <cmath>

#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <c_da.h>
#include <c_p_da.h>
#include <element.h>
#include <sext.h>

// ----------------------------------------------------------------------
// sextupole
// ----------------------------------------------------------------------

void Sext::Mapping(map_double& x)
{
  if(K_2==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=0.5*K_2*(x[0]*x[0]-x[1]*x[1]);
  x[4]+=K_2*x[0]*x[1];

  if(length!=0.) DriftMapping(length*0.5,x);
}
	 

void Sext::Mapping(map_da& x)
{
  if(K_2==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=0.5*K_2*(x[0]*x[0]-x[1]*x[1]);
  x[4]+=K_2*x[0]*x[1];

  if(length!=0.) DriftMapping(length*0.5,x);
}
	 

void Sext::Mapping(map_p_da& x)
{
  if(K_2==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=0.5*K_2*(x[0]*x[0]-x[1]*x[1]);
  x[4]+=K_2*x[0]*x[1];

  if(length!=0.) DriftMapping(length*0.5,x);
}


void Sext::Mapping(pBeam& x)
{
  if(K_2==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  for(int i=0;i<x.np;i++) {
    x.px[i]-=0.5*K_2*(x.x[i]*x.x[i]-x.y[i]*x.y[i]);
    x.py[i]+=K_2*x.x[i]*x.y[i];
  }

  if(length!=0.) DriftMapping(length*0.5,x);
}

	 
// -----------------------------------------------------------------------
// Octupole
// -----------------------------------------------------------------------

void Octu::Mapping(map_double& x)
{
  double f=K_3/6.;
  if(K_3==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=f*x[0]*(x[0]*x[0]-3.*x[1]*x[1]);
  x[4]+=f*x[1]*(3.*x[0]*x[0]-x[1]*x[1]);

  if(length!=0.) DriftMapping(length*0.5,x);
}
	 

void Octu::Mapping(map_da& x)
{
  double f=K_3/6.;
  if(K_3==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=f*x[0]*(x[0]*x[0]-3.*x[1]*x[1]);
  x[4]+=f*x[1]*(3.*x[0]*x[0]-x[1]*x[1]);

  if(length!=0.) DriftMapping(length*0.5,x);
}
	 

void Octu::Mapping(map_p_da& x)
{
  double f=K_3/6.;
  if(K_3==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  x[3]-=f*x[0]*(x[0]*x[0]-3.*x[1]*x[1]);
  x[4]+=f*x[1]*(3.*x[0]*x[0]-x[1]*x[1]);

  if(length!=0.) DriftMapping(length*0.5,x);
}


void Octu::Mapping(pBeam& x)
{
  double f=K_3/6.;
  if(K_3==0.) {
    if(length!=0.) DriftMapping(length,x);
    return;
  }

  if(length!=0.) DriftMapping(length*0.5,x);

  for(int i=0;i<x.np;i++) {
    x.px[i]-=f*x.x[i]*(x.x[i]*x.x[i]-3.*x.y[i]*x.y[i]);
    x.py[i]+=f*x.y[i]*(3.*x.x[i]*x.x[i]-x.y[i]*x.y[i]);
  }

  if(length!=0.) DriftMapping(length*0.5,x);
}
	 


