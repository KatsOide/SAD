#include <cmath>
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
#include <wiggler.h>

// Wiggler

Wiggler::Wiggler(char* s) : Element(s)
{
  By=get_parm(s,"By");
  Bx=get_parm(s,"Bx");
  N_pole=(int)get_parm(s,"N_pole");
  N_div=(int)get_parm(s,"N_div");
  dphase=get_parm(s,"dphase")*pi/180.;
 
  Kz=pi2*(double)N_pole/length;
  Kx=get_parm(s,"Kx");
  Ky=sqrt(Kz*Kz+Kx*Kx);
  Qy=get_parm(s,"Qy");
  Qx=sqrt(Kz*Kz+Qy*Qy);
  lambda=length/(double)N_pole;

  ds=length/(double)N_div; 
}

Wiggler::Wiggler(double* K) : Element(K)
{
}

void Wiggler::Mapping(map_double& x)
{

  double sinkz,sinkzp;
   
  // Following decralations should be changed for da calculation
  double kx1,kx2,kx3,kx4,f,g,fx,fy,gx,gy,dpds;
  double det1a,fgx,fgy,tmp;
   
  for(int i=0;i<N_div;i++) {
    kx1=Kx*x[0];
    kx2=Ky*x[1];
    kx3=Qx*x[0];
    kx4=Qy*x[1];
    sinkz=sin(Kz*(double)i*ds);
    sinkzp=sin(Kz*(double)i*ds+dphase);
    {
      f=F_0/Kz*cos(kx1)*cosh(kx2)*sinkz
	 -Qy*G_0/Qx/Kz*sinh(kx3)*sin(kx4)*sinkzp;
      g=-G_0/Kz*cos(kx4)*cosh(kx3)*sinkzp
	 +Kx*F_0/Ky/Kz*sinh(kx2)*sin(kx1)*sinkz;

      fx=-Kx*F_0/Kz*sin(kx1)*cosh(kx2)*sinkz
	 -Qy*G_0/Kz*cosh(kx3)*sin(kx4)*sinkzp;
      fy=Ky*F_0/Kz*cos(kx1)*sinh(kx2)*sinkz
	 -Qy*Qy*G_0/Qx/Kz*sinh(kx3)*cos(kx4)*sinkzp;

      gx=-G_0*Qx/Kz*cos(kx4)*sinh(kx3)*sinkzp
	 +Kx*Kx*F_0/Ky/Kz*sinh(kx2)*cos(kx1)*sinkz;
      gy=G_0*Qy/Kz*sin(kx4)*cosh(kx3)*sinkzp
	 +Kx*F_0/Kz*cosh(kx2)*sin(kx1)*sinkz;
    }
    {
      dpds=ds/(1.+x[5]);

      fx=dpds*fx;
      gx=dpds*gx;
      fy=dpds*fy;
      gy=dpds*gy;
      det1a=1.-fx-gy+fx*gy-fy*gx;

      fgx=dpds*(f*fx+g*gx);
      fgy=dpds*(f*fy+g*gy);

      tmp=fx;
      fx=(1.-gy)/det1a;
      gx=gx/det1a;
      fy=fy/det1a;
      gy=(1.-tmp)/det1a;
      
      tmp=x[3];
      x[3]=fx*(tmp-fgx)+gx*(x[4]-fgy);
      x[4]=fy*(tmp-fgx)+gy*(x[4]-fgy);

      fgx=(x[3]-f);
      fgy=(x[4]-g);
      x[0]+=dpds*fgx;
      x[1]+=dpds*fgy;
      x[2]-=0.5*(fgx*fgx+fgy*fgy)*(dpds*dpds)/ds;
    }
  }
}

void Wiggler::Mapping(map_da& x)
{
  double sinkz,sinkzp;

  // Following decralations should be changed for da calculation
  da kx1,kx2,kx3,kx4,f,g,fx,fy,gx,gy,dpds;
  da det1a,fgx,fgy,tmp;
   
  for(int i=0;i<N_div;i++) {
    kx1=Kx*x[0];
    kx2=Ky*x[1];
    kx3=Qx*x[0];
    kx4=Qy*x[1];
    sinkz=sin(Kz*(double)i*ds);
    sinkzp=sin(Kz*(double)i*ds+dphase);
    {
      f=F_0/Kz*cos(kx1)*cosh(kx2)*sinkz-
	 Qy*G_0/Qx/Kz*sinh(kx3)*sin(kx4)*sinkzp;
      g=-G_0/Kz*cos(kx4)*cosh(kx3)*sinkzp+
	 Kx*F_0/Ky/Kz*sinh(kx2)*sin(kx1)*sinkz;

      fx=-Kx*F_0/Kz*sin(kx1)*cosh(kx2)*sinkz-
	 Qy*G_0/Kz*cosh(kx3)*sin(kx4)*sinkzp;
      fy=Ky*F_0/Kz*cos(kx1)*sinh(kx2)*sinkz-
	 Qy*Qy*G_0/Qx/Kz*sinh(kx3)*cos(kx4)*sinkzp;

      gx=-G_0*Qx/Kz*cos(kx4)*sinh(kx3)*sinkzp+
	 Kx*Kx*F_0/Ky/Kz*sinh(kx2)*cos(kx1)*sinkz;
      gy=G_0*Qy/Kz*sin(kx4)*cosh(kx3)*sinkzp+
	 Kx*F_0/Kz*cosh(kx2)*sin(kx1)*sinkz;
    }
    {
      dpds=ds/(1.+x[5]);

      fx=dpds*fx;
      gx=dpds*gx;
      fy=dpds*fy;
      gy=dpds*gy;
      det1a=1.-fx-gy+fx*gy-fy*gx;

      fgx=dpds*(f*fx+g*gx);
      fgy=dpds*(f*fy+g*gy);

      tmp=fx;
      fx=(1.-gy)/det1a;
      gx=gx/det1a;
      fy=fy/det1a;
      gy=(1.-tmp)/det1a;
      
      tmp=x[3];
      x[3]=fx*(tmp-fgx)+gx*(x[4]-fgy);
      x[4]=fy*(tmp-fgx)+gy*(x[4]-fgy);

      fgx=(x[3]-f);
      fgy=(x[4]-g);
      x[0]+=dpds*fgx;
      x[1]+=dpds*fgy;
      x[2]-=0.5*(fgx*fgx+fgy*fgy)*(dpds*dpds)/ds;
    }
  }
}

void Wiggler::Mapping(map_p_da& x)
{
  double sinkz,sinkzp;

  // Following decralations should be changed for da calculation
  p_da kx1,kx2,kx3,kx4,f,g,fx,fy,gx,gy,dpds;
  p_da det1a,fgx,fgy,tmp;
   
  for(int i=0;i<N_div;i++) {
    kx1=Kx*x[0];
    kx2=Ky*x[1];
    kx3=Qx*x[0];
    kx4=Qy*x[1];
    sinkz=sin(Kz*(double)i*ds);
    sinkzp=sin(Kz*(double)i*ds+dphase);
    {
      f=F_0/Kz*cos(kx1)*cosh(kx2)*sinkz-
	 Qy*G_0/Qx/Kz*sinh(kx3)*sin(kx4)*sinkzp;
      g=-G_0/Kz*cos(kx4)*cosh(kx3)*sinkzp+
	 Kx*F_0/Ky/Kz*sinh(kx2)*sin(kx1)*sinkz;

      fx=-Kx*F_0/Kz*sin(kx1)*cosh(kx2)*sinkz-
	 Qy*G_0/Kz*cosh(kx3)*sin(kx4)*sinkzp;
      fy=Ky*F_0/Kz*cos(kx1)*sinh(kx2)*sinkz-
	 Qy*Qy*G_0/Qx/Kz*sinh(kx3)*cos(kx4)*sinkzp;

      gx=-G_0*Qx/Kz*cos(kx4)*sinh(kx3)*sinkzp+
	 Kx*Kx*F_0/Ky/Kz*sinh(kx2)*cos(kx1)*sinkz;
      gy=G_0*Qy/Kz*sin(kx4)*cosh(kx3)*sinkzp+
	 Kx*F_0/Kz*cosh(kx2)*sin(kx1)*sinkz;
    }
    {
      dpds=ds/(1.+x[5]);

      fx=dpds*fx;
      gx=dpds*gx;
      fy=dpds*fy;
      gy=dpds*gy;
      det1a=1.-fx-gy+fx*gy-fy*gx;

      fgx=dpds*(f*fx+g*gx);
      fgy=dpds*(f*fy+g*gy);

      tmp=fx;
      fx=(1.-gy)/det1a;
      gx=gx/det1a;
      fy=fy/det1a;
      gy=(1.-tmp)/det1a;
      
      tmp=x[3];
      x[3]=fx*(tmp-fgx)+gx*(x[4]-fgy);
      x[4]=fy*(tmp-fgx)+gy*(x[4]-fgy);

      fgx=(x[3]-f);
      fgy=(x[4]-g);
      x[0]+=dpds*fgx;
      x[1]+=dpds*fgy;
      x[2]-=0.5*(fgx*fgx+fgy*fgy)*(dpds*dpds)/ds;
    }
  }
}



void Wiggler::Mapping(pBeam& x)
{
  double sinkz,sinkzp;
  int j;
  // Following decralations should be changed for da calculation
  double kx1,kx2,kx3,kx4,f,g,fx,fy,gx,gy,dpds;
  double det1a,fgx,fgy,tmp;
   
  for(int i=0;i<N_div;i++) {
    for(j=0;i<x.np;j++) {
      kx1=Kx*x.x[j];
      kx2=Ky*x.y[j];
      kx3=Qx*x.x[j];
      kx4=Qy*x.y[j];
      sinkz=sin(Kz*(double)i*ds);
      sinkzp=sin(Kz*(double)i*ds+dphase);
      {
      f=F_0/Kz*cos(kx1)*cosh(kx2)*sinkz
	 -Qy*G_0/Qx/Kz*sinh(kx3)*sin(kx4)*sinkzp;
      g=-G_0/Kz*cos(kx4)*cosh(kx3)*sinkzp
	 +Kx*F_0/Ky/Kz*sinh(kx2)*sin(kx1)*sinkz;

      fx=-Kx*F_0/Kz*sin(kx1)*cosh(kx2)*sinkz
	 -Qy*G_0/Kz*cosh(kx3)*sin(kx4)*sinkzp;
      fy=Ky*F_0/Kz*cos(kx1)*sinh(kx2)*sinkz
	 -Qy*Qy*G_0/Qx/Kz*sinh(kx3)*cos(kx4)*sinkzp;

      gx=-G_0*Qx/Kz*cos(kx4)*sinh(kx3)*sinkzp
	 +Kx*Kx*F_0/Ky/Kz*sinh(kx2)*cos(kx1)*sinkz;
      gy=G_0*Qy/Kz*sin(kx4)*cosh(kx3)*sinkzp
	 +Kx*F_0/Kz*cosh(kx2)*sin(kx1)*sinkz;
      }
      {
      dpds=ds/(1.+x.pz[j]);

      fx=dpds*fx;
      gx=dpds*gx;
      fy=dpds*fy;
      gy=dpds*gy;
      det1a=1.-fx-gy+fx*gy-fy*gx;

      fgx=dpds*(f*fx+g*gx);
      fgy=dpds*(f*fy+g*gy);

      tmp=fx;
      fx=(1.-gy)/det1a;
      gx=gx/det1a;
      fy=fy/det1a;
      gy=(1.-tmp)/det1a;
      
      tmp=x.px[j];
      x.px[j]=fx*(tmp-fgx)+gx*(x.py[j]-fgy);
      x.py[j]=fy*(tmp-fgx)+gy*(x.py[j]-fgy);

      fgx=(x.px[j]-f);
      fgy=(x.py[j]-g);
      x.x[j]+=dpds*fgx;
      x.y[j]+=dpds*fgy;
      x.z[j]-=0.5*(fgx*fgx+fgy*fgy)*(dpds*dpds)/ds;
      }
    }
  }
}

