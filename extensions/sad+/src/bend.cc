#include <cmath>
using std::fabs;
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::atan;


#include <Complex.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <element.h>
#include <bend.h>

Bend::Bend(char* s) : Element(s)
{
  phi_0=get_parm(s,"phi0");
  dphi=get_parm(s,"dphi");
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
  sin_omega=sin(omega);
  cos_omega=cos(omega);
  sin_p2=sin(phi_0*0.5-e2);
  cos_p2=cos(phi_0*0.5-e2);
  if(phi_0!=0.) sin_pp=2.*sin(phi_0*0.5)/phi_0;
  else sin_pp=1.;
   // cout<< "initialize Bend\n";
}

Bend::Bend(double* K) : Element(K)
{
  phi_0=K[1];
  dphi=K[9];
  phi=phi_0+dphi;
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
  sin_omega=sin(omega);
  cos_omega=cos(omega);
  sin_p2=sin(phi_0*0.5-e2);
  cos_p2=cos(phi_0*0.5-e2);
  if(phi_0!=0.) sin_pp=2.*sin(phi_0*0.5)/phi_0;
  else sin_pp=1.;

      // cout<< "initialize Bend\n";
}

// Bending magnet

void Bend::Mapping(map_double& x)
{
  map_double x0;

  double pxy2,E,E2,pz,pz_i,Ep2,Ep_i,px_f,pz_f,lx_y,f,ff;
  double dpx,dpz,dpx_f,dpz_f,py2,psi1,psi2,psi11,psi22,psi3;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }
/*  if(phi_0==0.) {
    DriftMapping(length*0.5);
    
    DriftMapping(length*0.5);    
  }
*/
   // Entrance edge angle correction (coordinate transformation) 

  x[3]+=dphix;
  x[4]+=dphiy;

  E=1.+x[5];
  E2=E*E;
  if(e1!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    pz=sqrt(E2-pxy2);
    pz_i=1./pz;
    x[0]=x[0]/(cos_e1-x[3]*pz_i*sin_e1);
    x[1]+=x[4]*pz_i*x[0]*sin_e1;
    x[2]-=E*pz_i*x[0]*sin_e1;
    dpx=x[3]*cos_e1-pxy2/(pz+E)*sin_e1;
    x[3]=dpx+E*sin_e1;
   //   pz=pz*cos_e1-x[3]*sin_e1;
  }
  else dpx=x[3];

  // Entrance edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]+=E2*ff;
    x[4]-=x[3]*f;
    x[2]-=x[3]*E*ff;
  }
   // Bend body
  {
    py2=x[4]*x[4];
    Ep2=E2-py2;
   
    pz=sqrt(Ep2-x[3]*x[3]);
    dpz=(-dpx*(2.*E*sin_e1+dpx)-py2)/(pz+E*cos_e1);
    //     Ep_i=pow(Ep2,-0.5);
    dpx_f=cos_omega*dpx+sin_omega*dpz-sin_omega/rho*x[0]
      +(phi_0*x[5]-dphi)*sin_pp*cos_p2;
    px_f=-E*sin_e2+dpx_f;
    pz_f=sqrt(Ep2-px_f*px_f);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=cos_omega*x[0]+rho*(sin_omega*dpx-cos_omega*dpz
       +(phi_0*x[5]-dphi)*sin_pp*sin_p2+dpz_f);
    psi1=x[3]/pz-tan_e1;
    psi2=px_f/pz_f+tan_e2;
/*
   if(e1==0.) psi1=dpx/pz; else {
     psi1=(dpx*(dpx+2.*E*sin_e1)+py2*sin_e1*sin_e1)/
        (cos_e1*cos_e1*pz*(x[3]+pz*tan_e1));}
   if(e2==0.) psi2=px_f/pz_f; else {
     psi2=(dpx_f*(dpx_f-2.*E*sin_e2)+py2*sin_e2*sin_e2)/
       (cos_e2*cos_e2*pz_f*(px_f-pz_f*tan_e2));}
*/
    psi11=psi1/(1.+tan_e1*psi1+tan_e1*tan_e1);
    psi22=psi2/(1.-tan_e2*psi2+tan_e2*tan_e2);
    psi3=atan((psi11-psi22)/(1.+psi11*psi22));
    lx_y=rho*(phi_0+psi3);
    x[1]+=x[4]*lx_y;
    x[2]+=rho*(dphi-psi3-x[5]*(phi_0+psi3));
    x[3]=-E*sin_e2+dpx_f;
  }
   // Exit edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]-=ff*E2;
    x[4]+=x[3]*f;
    x[2]+=x[3]*E*ff;

//   x0.is_symplectic(" B exit edge kick");
  }
   // Exit edge angle transformation
  if(e2!=0.) {
    py2=x[4]*x[4];

    pz_f=sqrt(Ep2-py2);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=x[0]/(cos_e2-x[3]/pz_f*sin_e2);
    x[1]+=x[4]/pz_f*x[0]*sin_e2;
    x[2]-=E/pz_f*x[0]*sin_e2;
    x[3]=dpx_f*cos_e2+dpz_f*sin_e2;

   //   pz=pz*cos_e2-x[3]*sin_e2;
  }
  x[3]+=dphix;
  x[4]+=dphiy;
}


// Bending magnet

void Bend::Mapping(map_da& x)
{
  map_da x0;

  da pxy2,E,E2,pz,pz_i,Ep2,Ep_i,px_f,pz_f,lx_y,f,ff;
  da dpx,dpz,dpx_f,dpz_f,py2,psi1,psi2,psi11,psi22,psi3;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }

   // Entrance edge angle correction (coordinate transformation) 

  E=1.+x[5];
  E2=E*E;
  if(e1!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    pz=sqrt(E2-pxy2);
    pz_i=1./pz;
    x[0]=x[0]/(cos_e1-x[3]*pz_i*sin_e1);
    x[1]+=x[4]*pz_i*x[0]*sin_e1;
    x[2]-=E*pz_i*x[0]*sin_e1;
    dpx=x[3]*cos_e1-pxy2/(pz+E)*sin_e1;
    x[3]=dpx+E*sin_e1;
   //   pz=pz*cos_e1-x[3]*sin_e1;
  }
  else dpx=x[3];

   // Entrance edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]+=E2*ff;
    x[4]-=x[3]*f;
    x[2]-=x[3]*E*ff;
  }
   // Bend body
  {
    py2=x[4]*x[4];
    Ep2=E2-py2;

    pz=sqrt(Ep2-x[3]*x[3]);
    dpz=(-dpx*(2.*E*sin_e1+dpx)-py2)/(pz+E*cos_e1);
//     Ep_i=pow(Ep2,-0.5);
    dpx_f=cos_omega*dpx+sin_omega*dpz-sin_omega/rho*x[0]
      +(phi_0*x[5]-dphi)*sin_pp*cos_p2;;
    px_f=-E*sin_e2+dpx_f;
    pz_f=sqrt(Ep2-px_f*px_f);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=cos_omega*x[0]+rho*(sin_omega*dpx-cos_omega*dpz
       +(phi_0*x[5]-dphi)*sin_pp*sin_p2+dpz_f);
    psi1=x[3]/pz-tan_e1;
    psi2=px_f/pz_f+tan_e2;
/*
     psi1=(dpx*(dpx+2.*E*sin_e1)+py2*sin_e1*sin_e1)/
       (cos_e1*cos_e1*pz*(x[3]+pz*tan_e1));
     psi2=(dpx_f*(dpx_f-2.*E*sin_e2)+py2*sin_e2*sin_e2)/
       (cos_e2*cos_e2*pz_f*(px_f-pz_f*tan_e2));
*/
    psi11=psi1/(1.+tan_e1*psi1+tan_e1*tan_e1);
    psi22=psi2/(1.-tan_e2*psi2+tan_e2*tan_e2);
    psi3=atan((psi11-psi22)/(1.+psi11*psi22));
    lx_y=rho*(phi_0+psi3);
    x[1]+=x[4]*lx_y;
    x[2]+=rho*(dphi-psi3-x[5]*(phi_0+psi3));
    x[3]=-E*sin_e2+dpx_f;
  }
//  cout << "Out body" << x;
//  is_symplectic(x);
   // Exit edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]-=ff*E2;
    x[4]+=x[3]*f;
    x[2]+=x[3]*E*ff;

//   x0.is_symplectic(" B exit edge kick");
  }
   // Exit edge angle transformation
  if(e2!=0.) {
    py2=x[4]*x[4];

    pz_f=sqrt(Ep2-py2);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=x[0]/(cos_e2-x[3]/pz_f*sin_e2);
    x[1]+=x[4]/pz_f*x[0]*sin_e2;
    x[2]-=E/pz_f*x[0]*sin_e2;
    x[3]=dpx_f*cos_e2+dpz_f*sin_e2;

   //   pz=pz*cos_e2-x[3]*sin_e2;
  }
//  cout << "Fin";
//  is_symplectic(x);
}


// Bending magnet

void Bend::Mapping(map_p_da& x)
{
  map_p_da x0;

  p_da pxy2,E,E2,pz,pz_i,Ep2,Ep_i,px_f,pz_f,lx_y,f,ff;
  p_da dpx,dpz,dpx_f,dpz_f,py2,psi1,psi2,psi11,psi22,psi3;

  if(phi==0.) {
    DriftMapping(length,x);
    return;
  }
  if(length==0.) { 
    x[3]+=phi_0*x[5]-dphi; x[2]-=phi_0*x[0];
    return;
  }

   // Entrance edge angle correction (coordinate transformation) 

  E=1.+x[5];
  E2=E*E;
  if(e1!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    pz=sqrt(E2-pxy2);
    pz_i=1./pz;
    x[0]=x[0]/(cos_e1-x[3]*pz_i*sin_e1);
    x[1]+=x[4]*pz_i*x[0]*sin_e1;
    x[2]-=E*pz_i*x[0]*sin_e1;
    dpx=x[3]*cos_e1-pxy2/(pz+E)*sin_e1;
    x[3]=dpx+E*sin_e1;
   //   pz=pz*cos_e1-x[3]*sin_e1;
  }
  else dpx=x[3];

   // Entrance edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]+=E2*ff;
    x[4]-=x[3]*f;
    x[2]-=x[3]*E*ff;
  }
   // Bend body
  {
    py2=x[4]*x[4];
    Ep2=E2-py2;

    pz=sqrt(Ep2-x[3]*x[3]);
    dpz=(-dpx*(2.*E*sin_e1+dpx)-py2)/(pz+E*cos_e1);
//     Ep_i=pow(Ep2,-0.5);
    dpx_f=cos_omega*dpx+sin_omega*dpz-sin_omega/rho*x[0]
      +(phi_0*x[5]-dphi)*sin_pp*cos_p2;;
    px_f=-E*sin_e2+dpx_f;
    pz_f=sqrt(Ep2-px_f*px_f);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=cos_omega*x[0]+rho*(sin_omega*dpx-cos_omega*dpz
       +(phi_0*x[5]-dphi)*sin_pp*sin_p2+dpz_f);
    psi1=x[3]/pz-tan_e1;
    psi2=px_f/pz_f+tan_e2;
/*
     psi1=(dpx*(dpx+2.*E*sin_e1)+py2*sin_e1*sin_e1)/
       (cos_e1*cos_e1*pz*(x[3]+pz*tan_e1));
     psi2=(dpx_f*(dpx_f-2.*E*sin_e2)+py2*sin_e2*sin_e2)/
       (cos_e2*cos_e2*pz_f*(px_f-pz_f*tan_e2));
*/
    psi11=psi1/(1.+tan_e1*psi1+tan_e1*tan_e1);
    psi22=psi2/(1.-tan_e2*psi2+tan_e2*tan_e2);
    psi3=atan((psi11-psi22)/(1.+psi11*psi22));
    lx_y=rho*(phi_0+psi3);
    x[1]+=x[4]*lx_y;
    x[2]+=rho*(dphi-psi3-x[5]*(phi_0+psi3));
    x[3]=-E*sin_e2+dpx_f;
  }
   // Exit edge field map
  {
    Ep2=E2-x[3]*x[3];
    Ep_i=pow(Ep2,-0.5);
    f=x[1]/rho*Ep_i;
    ff=f*x[1]*0.5/Ep2;
    x[0]-=ff*E2;
    x[4]+=x[3]*f;
    x[2]+=x[3]*E*ff;

//   x0.is_symplectic(" B exit edge kick");
  }
   // Exit edge angle transformation
  if(e2!=0.) {
    py2=x[4]*x[4];

    pz_f=sqrt(Ep2-py2);
    dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
    x[0]=x[0]/(cos_e2-x[3]/pz_f*sin_e2);
    x[1]+=x[4]/pz_f*x[0]*sin_e2;
    x[2]-=E/pz_f*x[0]*sin_e2;
    x[3]=dpx_f*cos_e2+dpz_f*sin_e2;

   //   pz=pz*cos_e2-x[3]*sin_e2;
  }
}



// Bending magnet

void Bend::Mapping(pBeam& x)
{
  int i;
  double pxy2,E,E2,pz,pz_i,Ep2,Ep_i,px_f,pz_f,lx_y,f,ff;
  double dpx,dpz,dpx_f,dpz_f,py2,psi1,psi2,psi11,psi22,psi3;

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

   // Entrance edge angle correction (coordinate transformation) 

  for(i=0;i<x.np;i++) {
    E=1.+x.pz[i];
    E2=E*E;
    if(e1!=0.) {
      pxy2=x.px[i]*x.px[i]+x.py[i]*x.py[i];
      pz=sqrt(E2-pxy2);
      pz_i=1./pz;
      x.x[i]=x.x[i]/(cos_e1-x.px[i]*pz_i*sin_e1);
      x.y[i]=x.y[i]+x.py[i]*pz_i*x.x[i]*sin_e1;
      x.z[i]=x.z[i]-E*pz_i*x.x[i]*sin_e1;
      dpx=x.px[i]*cos_e1-pxy2/(pz+E)*sin_e1;
      x.px[i]=dpx+E*sin_e1;
    }
    else {
      dpx=x.px[i];
    }

   // Entrance edge field map
    {
      Ep2=E2-x.px[i]*x.px[i];
      Ep_i=pow(Ep2,-0.5);
      f=x.y[i]/rho*Ep_i;
      ff=f*x.y[i]*0.5/Ep2;
      x.x[i]+=E2*ff;
      x.py[i]-=x.px[i]*f;
      x.z[i]-=x.px[i]*E*ff;
   
    }
   // Bend body
    {
      py2=x.py[i]*x.py[i];
      Ep2=E2-py2;
   
      pz=sqrt(Ep2-x.px[i]*x.px[i]);
      dpz=(-dpx*(2.*E*sin_e1+dpx)-py2)/(pz+E*cos_e1);
//     Ep_i=pow(Ep2,-0.5);
      dpx_f=cos_omega*dpx+sin_omega*dpz-sin_omega/rho*x.x[i]
	+(phi_0*x.pz[i]-dphi)*sin_pp*cos_p2;
      px_f=-E*sin_e2+dpx_f;
      pz_f=sqrt(Ep2-px_f*px_f);
      dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
      x.x[i]=cos_omega*x.x[i]+rho*(sin_omega*dpx-cos_omega*dpz
	 +(phi_0*x.pz[i]-dphi)*sin_pp*sin_p2+dpz_f);
      psi1=x.px[i]/pz-tan_e1;
      psi2=px_f/pz_f+tan_e2;
/*
      psi1=(dpx*(dpx+2.*E*sin_e1)+py2*sin_e1*sin_e1)/
	(cos_e1*cos_e1*pz*(x.px[i]+pz*tan_e1));
      psi2=(dpx_f*(dpx_f-2.*E*sin_e2)+py2*sin_e2*sin_e2)/
	(cos_e2*cos_e2*pz_f*(px_f-pz_f*tan_e2));
*/
      psi11=psi1/(1.+tan_e1*psi1+tan_e1*tan_e1);
      psi22=psi2/(1.-tan_e2*psi2+tan_e2*tan_e2);
      psi3=atan((psi11-psi22)/(1.+psi11*psi22));
      lx_y=rho*(phi_0+psi3);
      x.y[i]+=x.py[i]*lx_y;
      x.z[i]+=rho*(dphi-psi3-x.pz[i]*(phi_0+psi3));
      x.px[i]=-E*sin_e2+dpx_f;
    }
   // Exit edge field map
    {
      Ep2=E2-x.px[i]*x.px[i];
      Ep_i=pow(Ep2,-0.5);
      f=x.y[i]/rho*Ep_i;
      ff=f*x.y[i]*0.5/Ep2;
      x.x[i]-=ff*E2;
      x.py[i]+=x.px[i]*f;
      x.z[i]+=x.px[i]*E*ff;

//   x.is_symplectic(" B exit edge kick");
   }
    // Exit edge angle transformation
    if(e2!=0.) {
      py2=x.py[i]*x.py[i];

      pz_f=sqrt(Ep2-py2);
      dpz_f=(dpx_f*(2.*E*sin_e2-dpx_f)-py2)/(pz_f+E*cos_e2);
      x.x[i]=x.x[i]/(cos_e2-x.px[i]/pz_f*sin_e2);
      x.y[i]+=x.py[i]/pz_f*x.x[i]*sin_e2;
      x.z[i]-=E/pz_f*x.x[i]*sin_e2;
      x.px[i]=dpx_f*cos_e2+dpz_f*sin_e2;
   }
   //   pz=pz*cos_e2-x.m[3]*sin_e2;
// End np loop
  }
}
