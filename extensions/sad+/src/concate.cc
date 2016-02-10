#include <c_da.h>
#include <p_da.h>
#include <map_da.h>
#include <map_c_da.h>
#include <map_p_da.h>
#include <track.h>

da concatenate(const da& x,const map_da& y)
{
  da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  double xi;
  da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1*=y.m[j1];
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

da concatenate(const da& x,const map_da& y,int No)
{
  da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  double xi;
  da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1=muln(y1,y.m[j1],No);
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}


map_da concatenate(const map_da& x,const map_da& y)
{
  map_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  double xi;
  da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1*=y.m[j1];
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

map_da concatenate(const map_da& x,const map_da& y,int No)
{
  map_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  double xi;
  da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1=muln(y1,y.m[j1],No);
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}
//
// da*double  used for DA tracking
//
double concatenate(const da& x,const map_double& y)
{
  double z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  double xi;
  double y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1*=y.m[j1];
	  }
        }
	xi=vget(x,i);
	if(xi!=0.) z+=xi*y1*y2;
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

map_double concatenate(const map_da& x,const map_double& y)
{
  map_double z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  double xi;
  double y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1*=y.m[j1];
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

pBeam concatenate(const map_da& x,const pBeam& y)
{
  pBeam z(y.N_particle,y.rn_particle);
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A,ip;
  int j1,j2,k1,k2,icv;
  double xi;
  double *y2,*y1;
  int N_cv2=N_cv/2;

  y1=new double [y.np];
  y2=new double [y.np];
  
  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(ip=0;ip<y.np;ip++) {
         y2[ip]=1.;
         for(j2=0;j2<N_cv2;j2++) {
	    for(k2=0;k2<iv2[j2];k2++) {
	       y2[ip]*=y.x[(j2+N_cv2)*y.N_particle+ip];
	    }
         }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	for(ip=0;ip<y.np;ip++) {
	   y1[ip]=1.;
           for(j1=0;j1<N_cv2;j1++) {
	     for(k1=0;k1<iv1[j1];k1++) {
	       y1[ip]*=y.x[j1*y.N_particle+ip];
	     }
           } 
	   for(icv=0;icv<N_cv;icv++) {
	     xi=vget(x.m[icv],i);
	     if(xi!=0.) z.x[icv*y.N_particle+ip]+=xi*y1[ip]*y2[ip];
	   }
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  delete [] y2;
  delete [] y1;
  z.np=y.np;
  return z;
}
//---------------------------------------------------------------
// Real*Complex 
//---------------------------------------------------------------

c_da concatenate(const da& x,const map_c_da& y)
{
  c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  double xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1*=y.m[j1];
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}


c_da concatenate(const da& x,const map_c_da& y,int No)
{
  c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  double xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1=muln(y1,y.m[j1],No);
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}


map_c_da concatenate(const map_da& x,const map_c_da& y)
{
  map_c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  double xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1*=y.m[j1];
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

map_c_da concatenate(const map_da& x,const map_c_da& y,int No)
{
  map_c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  double xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1=muln(y1,y.m[j1],No);
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

//---------------------------------------------------------------
// Complex 
//---------------------------------------------------------------

c_da concatenate(const c_da& x,const map_c_da& y)
{
  c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  Complex xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1*=y.m[j1];
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}


c_da concatenate(const c_da& x,const map_c_da& y,int No)
{
  c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2;
  Complex xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=vget(x,i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);
	  y1=1.;
          for(j1=0;j1<N_cv2;j1++) {
	    for(k1=0;k1<iv1[j1];k1++) {
	      y1=muln(y1,y.m[j1],No);
	    }
          }
	  z+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}


map_c_da concatenate(const map_c_da& x,const map_c_da& y)
{
  map_c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  Complex xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2*=y.m[j2+N_cv2];
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1*=y.m[j1];
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

map_c_da concatenate(const map_c_da& x,const map_c_da& y,int No)
{
  map_c_da z;
  int C2M,C1M;
  int i,i1A,i2A,N2A,N1A;
  int j1,j2,k1,k2,icv;
  Complex xi;
  c_da y2,y1;
  int N_cv2=N_cv/2;

  z=0.;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      y2=1.;
      for(j2=0;j2<N_cv2;j2++) {
	 for(k2=0;k2<iv2[j2];k2++) {
	    y2=muln(y2,y.m[j2+N_cv2],No);
	 }
      }
      
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);
	y1=1.;
        for(j1=0;j1<N_cv2;j1++) {
	  for(k1=0;k1<iv1[j1];k1++) {
	    y1=muln(y1,y.m[j1],No);
	  }
        }
	for(icv=0;icv<N_cv;icv++) {
	  xi=vget(x.m[icv],i);
	  if(xi!=0.) z.m[icv]+=xi*y1*y2;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return z;
}

map_p_da concatenate(const map_p_da& x,const map_p_da& y)
{
   return x;
}

