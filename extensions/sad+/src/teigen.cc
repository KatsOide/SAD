#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
using std::fabs;
using std::sqrt;

#include <Complex.h>
void teigen(double*,double*,int,int);
void tbal(double*,double*,double*,int,int);
void thess(double*,double*,double*,int,int);
void awe(double*,double*,double*,int,int);
void tqr(double*,double*,double*,int*,double*,int,int);
inline double max(double x,double y)
{
  if(x>y) return x; else return y;
}
inline double min(double x,double y)
{
  if(x<y) return x; else return y;
}
inline int max(int x,int y)
{
  if(x>y) return x; else return y;
}
inline int min(int x,int y)
{
  if(x<y) return x; else return y;
}
inline double sign(double x,double y)
{ if(y>0) return fabs(x); else return -fabs(x); }

void teigen(double* a,double* eig,int n,int ndim)
{
/*

   Subroutine to obtain eigen values and eigen vectors of a real
   matrix a.
					 7-Dec-1987   K. Oide
   Modified				27-Oct-1989   K. Oide
   Modified				17-Nov-1989   K. Oide
   Added Balancing			 1-Nov-1990   K. Oide
   Modified Convergence			20-Nov-1990   K. Oide
   Modified Shifting			 4-Dec-1990   K. Oide
   Fast Givens' Transformation		21-Jan-1992   K. Oide
   C++ version				K. Ohmi

   Usage:       call teigen(a,w,eig,n,ndim)

		where

		a     is the input n*n matrix and the eigen vectors
		      are returned here.
		w     is a work area of n*n words.
		eig   is a (2,n) real array where the eigen values
		      are returned.   eig(1,i) and eig(2,i) contain
		      the real and imaginary part of the i-th eigen-
		      valeu respectively.   When the eigenvalue is
		      a Complex, the output are stored as

		      eig(1,i  ) =      er
		      eig(2,i  ) =      ei
		      eig(1,i+1) =      ec 
		      eig(2,i+1) =     -ei

		      where er and ei are the real and imaginary
		      part of the i-th eigenvalue respectively.
		n     is the matrix dimension.
		ndim  is the size of the first dimension of a.

   Restriction: When the matrix cannot be diagonalized (Jordan
		standard form), the returned eigen vectors are not
		correct.

*/
   int n1;
   double* w=new double[n*n];
   int* ibtab=new int[n];
//      real*8 a(ndim,n),w(n,n),eig(2,n)
   n1=(n+1)/2;
   tbal(a,w,eig+2*n1,n,ndim);
   thess(w,a,eig+2*n1,n,ndim);
   tqr(w,a,eig,ibtab,eig+2*n1,n,ndim);
   delete [] ibtab;
   delete [] w;
}

void tbal(double* w,double* a,double* v,int n,int ndim)
{
   enum logical {F=0,T=1};
    const double radix=16.,sqrdx=256.;
    //   radix is the machine radix of floating.
    //  parameter (radix=16.d0,sqrdx=radix**2)
   int i,j;
   double c,r,g,f,s,u;
   logical last=F; 
   for(i=0;i<n;i++){
      for(j=0;j<n;j++) {
	  a[i*n+j]=w[i*ndim+j];
	  w[i*ndim+j]=0.;
      }
      w[i*ndim+i]=1.;
      v[i]=1.;
   }
   do {
      last=F;
      for(i=0;i<n;i++) {
	c=0.;  r=0.;
	for(j=0;j<n;j++) {
	  if(j!=i) {
	    u=v[j]/v[i];
	    c=c+fabs(a[i*n+j]*u);
	    r=r+fabs(a[j*n+i]/u);
	  }
	}
	if(c!=0. && r!= 0.) {
	  g=r/radix;
	  f=1.;
	  s=c+r;
	  while(c < g) {
	    f=f*radix;
	    c=c*sqrdx;
	  }
	  g=r*radix;
	  while(c > g) {
	    f=f/radix;
	    c=c/sqrdx;
	  }
	  if((c+r)/f < .95*s) {
	    last=T;
	    v[i]=v[i]/f;
	  }
	}
      }
   } while(last);
}



void thess(double* a,double* w,double* v,int n,int ndim)
{
   int i,i1,j,k;
   double h1,aa,ww,p,q;

   for(i=1;i<n;i++) v[i]=v[i]*v[i];

   for(i=0;i<n-2;i++) {
      i1=i+1;
      for(j=i1+1;j<n;j++) {
	  if(fabs(a[i*n+i1]) > fabs(a[i*n+j])) {
	    p=a[i*n+j]/a[i*n+i1];
	    h1=v[i1]+v[j]*p*p;
	    q=v[j]*p/h1;
	    v[j]=v[i1]*v[j]/h1;
	    v[i1]=h1;
	    a[i*n+j]=0.;
	    for(k=i1;k<n;k++) {
	      a[k*n+j]=a[k*n+j]-p*a[k*n+i1];
	      a[k*n+i1]=a[k*n+i1]+q*a[k*n+j];
	    }
	    for(k=0;k<n;k++) {
	      a[i1*n+k]=a[i1*n+k]+p*a[j*n+k];
	      a[j*n+k]=a[j*n+k]-q*a[i1*n+k];
	      w[i1*ndim+k]=w[i1*ndim+k]+p*w[j*ndim+k];
	      w[j*ndim+k]=w[j*ndim+k]-q*w[i1*ndim+k];
	    }
	  }
	  else if(a[i*n+j] != 0.) {
	    p=a[i*n+i1]/a[i*n+j];
	    h1=v[j]+v[i1]*p*p;
	    q=v[i1]*p/h1;
	    v[j]=v[i1]*v[j]/h1;
	    v[i1]=h1;
	    for(k=i1;k<n;k++) {
	      aa=a[j+ndim*k];
	      a[k*n+j]=p*aa-a[k*n+i1];
	      a[k*n+i1]=aa-q*a[k*n+j];
	    }
	    a[i*n+i1]=a[i*n+j];
	    a[i*n+j]=0.;
	    for(k=0;k<n;k++) {
	      aa=a[i1*n+k];
	      a[i1*n+k]=p*aa+a[j*n+k];
	      a[j*n+k]=q*a[i1*n+k]-aa;
	      ww=w[i1*ndim+k];
	      w[i1*ndim+k]=p*ww+w[j*ndim+k];
	      w[j*ndim+k]=q*w[i1*ndim+k]-ww;
	    }
	  }
      }
   }
   for(i=1;i<n;i++) v[i]=sqrt(v[i]);
}


//
//
void awe(double* a,double* w,double* eig,int n,int ndim)
{
   int j,k;
   for(k=0;k<n;k++) {
      for(j=0;j<n;j++) cout << a[j+ndim*k] << " ";
      cout << endl;
   }
      cout << endl;
   for(k=0;k<n;k++) {
      for(j=0;j<n;j++) cout << w[j+n*k] << " ";
      cout << endl;
   }
      cout << endl;
   for(k=0;k<n;k++) {
      for(j=0;j<2;j++) cout << eig[j+2*k] << " ";
      cout << endl;
   }
      cout << endl;
}
//

void tqr(double* a,double* w,double* eig,int* ibtab,double* vx,int n,int ndim)
{
   const int itmax=30;
   const double vmax=1.e10,vmin=1./vmax,bconv=1.7e-8;
   enum logical {F=0,T=1};
   logical paired,jordan;
   Complex ca,cc,cu,cr,ck,cm;
   int i,j,k,i1,i2,i3,is,is1,is2,j1,jm,ie1,iter;
   double a1,a2,w1,w2,w3,u,v,am,c,s,v1,v2,p,q,aa,ww,r,am1;
   double x,y,ee;

   jordan=F;
   int ibp=0;
   ibtab[0]=0;
   int ib=0;
   int ie=n-1;
   for(i=0;i<=ie-2;i++) {
	a[i*n+i+2]=0.;
	a[i*n+min(i+3,ie)]=0.;
   }
   double anorm=0.;
   for(i=0;i<n;i++) {
      for(j=max(0,i-1);j<n;j++) anorm+=fabs(vx[i]/vx[j]*a[j*n+i]);
   }
   anorm=anorm*2/n/n;
L1: 
   if(ie <= ib+1) {
      if(ib==0) {
	 goto L2000;
      }
      else {
	 ie=ib-1;
	 ibp=ibp-1;
	 ib=ibtab[ibp];
	 goto L1;
      }
   }
   ie1=ie-1;
   iter=0;
L2:
//       write(*,'(1P4G15.7)')((vx(i)/vx(j)*a(i,j),j=1,4),i=1,4)
   for(i=ie-1;i>=ib;i--) {
      i1=i+1;
      double s=max(anorm,fabs(a[i1+ndim*i1])+fabs(a[i+ndim*i]));
      if(fabs(vx[i1]/vx[i]*a[i1+ndim*i])+s == s) {
	  a[i1+ndim*i]=0.;
	  if(i1==ie) {
	    ie--;
	    goto L1;
	  }
	  ibp=ibp+1;
	  ibtab[ibp]=ib;
	  ib=i1;
	  goto L1;
      }
   }
   for(i=ib;i<=ie;i++) {
      if(fabs(vx[i])< vmin || fabs(vx[i])> vmax) {
	 for(j=max(i-1,ib);j<n;j++) 
	    a[i+j*n]=a[i+j*ndim]*vx[i];
	 for(j=0;j<min(ie,i+1);j++) 
	    a[j+i*ndim]=a[j+i*ndim]/vx[i];
	 for(j=0;j<n;j++) w[j+i*n]=w[j+i*n]/vx[i];
	 vx[i]=1.;
      }
   }
   w1=0.;w2=0.;w3=0.;
   for(is=ie-2,is1=ie-1,is2=ie;is>=ib;is--,is1--,is2--) {
      a1=a[ie+ie*n]-a[is+is*n];
      a2=a[ie1+ie1*n]-a[is+is*n];
      w1=((a1*a2-a[ie1+ie*n]*a[ie+ie1*n])/a[is1+ndim*is]+a[is+is1*n])
	   *vx[is]/vx[is1];
      w2=a[is1+ndim*is1]-a[is+ndim*is]-a1-a2;
      w3=a[is2+ndim*is1]*vx[is2]/vx[is1];
      if(is== ib) break;
      u=fabs(a[is+ndim*(is-1)]*vx[is]/vx[is-1])*(fabs(w2)+fabs(w3));
      v=fabs(w1)*(fabs(a[is-1+ndim*(is-1)])+fabs(a[is+ndim*is])+fabs(a[is1+ndim*is1]));
      if(u+v == v){
	  a[is1+ndim*(is-1)]=0.;
	  a[is2+ndim*(is-1)]=0.;
	  break;
      }
   }
//4511  continue
   a[is2+ndim*is]=0.;
   a[min(is2+1,ie)+ndim*is]=0.;
   am =abs(Complex(w1,w2));
   am1=abs(Complex(am,w3));
   if(is > ib) a[is+ndim*(is-1)]=a[is+ndim*(is-1)]*w1/am1*vx[is];

   for(i=is-1;i<=ie-2;i++){ // do 210 i=is-1,ie-2
	i1=i+1;
	i2=i1+1;
	i3=i2+1;
	if(i1> is) {
	  w1=vx[i1]*a[i1+ndim*i];
	  w2=vx[i2]*a[i2+ndim*i];
	  a[i2+ndim*i]=0.;
	  am =abs(Complex(w1,w2));
	  if(i3 <= ie) {
	    w3=vx[i3]*a[i3+ndim*i];
	    a[i3+ndim*i]=0.;
	    am1=abs(Complex(am,w3));
	  }
	  else {
	    w3=0.;
	    am1=0.;
	  }
	}
	if(am != 0.) {
	  c=w1/am;
	  s=w2/am;
	  if(fabs(c)>fabs(s)) {
	    v1=vx[i1]/c;
	    v2=c*vx[i2];
	    p=s*vx[i1]/v2;
	    q=s*vx[i2]/v1;
	    vx[i1]=v1;
	    vx[i2]=v2;
	    for(j=i1;j<n;j++) {
	      a[i2+ndim*j]=a[i2+ndim*j]-p*a[i1+ndim*j];
	      a[i1+ndim*j]=a[i1+ndim*j]+q*a[i2+ndim*j];
	    }
	    for(j=0;j<=min(ie,i3);j++) {
	      a[j+ndim*i1]=a[j+ndim*i1]+p*a[j+ndim*i2];
	      a[j+ndim*i2]=a[j+ndim*i2]-q*a[j+ndim*i1];
	    }
	    for(j=0;j<n;j++) {
	      w[j+n*i1]=w[j+n*i1]+p*w[j+n*i2];
	      w[j+n*i2]=w[j+n*i2]-q*w[j+n*i1];
	    }
	  }
	  else {
	    v1=vx[i2]/s;
	    v2=vx[i1]*s;
	    if(i1 > is) a[i1+ndim*i]=w2/vx[i2];
	    p=c*vx[i2]/v2;
	    q=c*vx[i1]/v1;
	    vx[i1]=v1;
	    vx[i2]=v2;
	    for(j=i1;j<n;j++) {
	      aa=a[i2+ndim*j];
	      a[i2+ndim*j]=p*aa-a[i1+ndim*j];
	      a[i1+ndim*j]=aa-q*a[i2+ndim*j];
	    }
	    for(j=0;j<=min(ie,i3);j++) {
	      aa=a[j+ndim*i1];
	      a[j+ndim*i1]= p*aa+a[j+ndim*i2];
	      a[j+ndim*i2]= q*a[j+ndim*i1]-aa;
	    }
	    for(j=0;j<n;j++) {
	      ww=w[j+n*i1];
	      w[j+n*i1]= p*ww+w[j+n*i2];
	      w[j+n*i2]= q*w[j+n*i1]-ww;
	    }
	  }
	}
	if(am1 != 0.) {
	  c=am/am1;
	  s=w3/am1;
	  if(fabs(c) > fabs(s)) {
	    v1=vx[i1]/c;
	    v2=c*vx[i3];
	    p=s*vx[i1]/v2;
	    q=s*vx[i3]/v1;
	    vx[i1]=v1;
	    vx[i3]=v2;
	    for(j=i1;j<n;j++) {
	      a[i3+ndim*j]=a[i3+ndim*j]-p*a[i1+ndim*j];
	      a[i1+ndim*j]=a[i1+ndim*j]+q*a[i3+ndim*j];
	    }
	    for(j=0;j<=min(ie,i3+1);j++) {
	      a[j+ndim*i1]=a[j+ndim*i1]+p*a[j+ndim*i3];
	      a[j+ndim*i3]=a[j+ndim*i3]-q*a[j+ndim*i1];
	    }
	    for(j=0;j<n;j++) {
	      w[j+n*i1]=w[j+n*i1]+p*w[j+n*i3];
	      w[j+n*i3]=w[j+n*i3]-q*w[j+n*i1];
	    }
	  } 
	  else {
	    v1=vx[i3]/s;
	    v2=vx[i1]*s;
	    if(i1 > is) a[i1+ndim*i]=w3/vx[i3];

	    p=c*vx[i3]/v2;
	    q=c*vx[i1]/v1;
	    vx[i1]=v1;
	    vx[i3]=v2;
	    for(j=i1;j<n;j++) {
	      aa=a[i3+ndim*j];
	      a[i3+ndim*j]=p*aa-a[i1+ndim*j];
	      a[i1+ndim*j]=aa-q*a[i3+ndim*j];
	    }
	    for(j=0;j<=min(ie,i3+1);j++) {
	      aa=a[j+ndim*i1];
	      a[j+ndim*i1]= p*aa+a[j+ndim*i3];
	      a[j+ndim*i3]= q*a[j+ndim*i1]-aa;
	    }
	    for(j=0;j<n;j++) {
	      ww=w[j+n*i1];
	      w[j+n*i1]= p*ww+w[j+n*i3];
	      w[j+n*i3]= q*w[j+n*i1]-ww;
	    }
	  }
	}
//c       write(*,'(1P4G15.7)')((vx(ii)/vx(jj)*a(ii,jj),jj=1,4),ii=1,4)
   } // 210   continue
   if(is > ib) a[is+ndim*(is-1)]=a[is+ndim*(is-1)]/vx[is];

   iter=iter+1;
   if(iter > itmax) {
        cout << " TEIGEN convergence failed. Range = " << ib << ie << endl;
	cout << "	Lower right corner = " << a[ie-1+ndim*(ie-1)]
	  //	     << " " << vx[ie]/vx[ie-1]*a[ie+ndim*(ie-1)]
	  //	     << " " << a[ie+ndim*ie]
	     << endl;
	a[ie+ndim*(ie-1)]=0.;
	iter=0;
   }
   goto L2;
L2000:  
   for(i=0;i<n;i++) {
	for(j=max(i-1,0);j<n;j++) a[i+n*j]=a[i+n*j]*vx[i];
	for(j=0;j<=min(n-1,i+1);j++) a[j+n*i]=a[j+n*i]/vx[i];
	for(j=0;j<n;j++) w[j+ndim*i]=w[j+ndim*i]/vx[i];
    }
//     write(*,'(1X,1P6G12.5)')a,w
    paired=F;
      
      
    for(i=0;i<n-1;i++) {  // 2010
	if(paired) {
	  paired=F;
	  goto L2010;
	}
	i1=i+1;
	if(a[i1+ndim*i] != 0.) {
	  if(fabs(a[i1+ndim*i]) > fabs(a[i+ndim*i1])) {
	    for(j=i;j<n;j++) {
	      aa=a[i+ndim*j];
	      a[i+ndim*j]=a[i1+ndim*j];
	      a[i1+ndim*j]=aa;
	    }
	    for(j=0;j<=i1;j++) {
	      aa=a[j+ndim*i];
	      a[j+ndim*i]=a[j+ndim*i1];
	      a[j+ndim*i1]=aa;
	    }
	    for(j=0;j<n;j++) {
	      ww=w[j+n*i];
	      w[j+n*i]=w[j+n*i1];
	      w[j+n*i1]=ww;
	    }
	  }
	  if(a[i1+ndim*i] == 0.) {
	    eig[2*i]=a[i+ndim*i];
	    eig[1+2*i]=0.;
	    goto L2010;
	  }
	  s=a[i+ndim*i]-a[i1+ndim*i1];
	  double decc=s*s+4.*a[i1+ndim*i]*a[i+ndim*i1];
	  if(decc >= 0.) {
	    double sqrd=sqrt(decc);
	    p=2.*a[i1+ndim*i]/(s+sign(sqrd,s));
	    for(j=i;j<n;j++) a[i1+ndim*j]=a[i1+ndim*j]-p*a[i+ndim*j];
	    for(j=0;j<=i1;j++) a[j+ndim*i]=a[j+ndim*i]+p*a[j+ndim*i1];
	    for(j=0;j<n;j++) w[j+n*i]=w[j+n*i]+p*w[j+n*i1];

	    eig[2*i]=a[i+ndim*i];
	    eig[2*i1]=a[i1+ndim*i1];
	    eig[1+2*i]=0.;
	    eig[1+2*i1]=0.;
	    a[i1+ndim*i]=0.;
	  } 
	  else {
	    p=-s/a[i+ndim*i1]*.5;
// ??
	    r=sign(sqrt(fabs((a[i1+ndim*i]-
	       p*(a[i+ndim*i]-a[i1+ndim*i1]+p*a[i+ndim*i1]))/a[i+ndim*i1])),
	       a[i+ndim*i1]);
	    for(j=i;j<n;j++) {
	      a[i1+ndim*j]=a[i1+ndim*j]-p*a[i+ndim*j];
	      a[i +ndim*j]=a[i +ndim*j]*r;
	    }
	    for(j=0;j<=i1;j++) a[j+ndim*i]=(a[j+ndim*i]+p*a[j+ndim*i1])/r;
	    for(j=0;j<n;j++) w[j+n*i]=(w[j+n*i]+p*w[j+n*i1])/r;
	    eig[2*i ]=a[i +ndim*i ];
	    eig[2*i1]=a[i1+ndim*i1];
	    eig[1+2*i ]=a[i +ndim*i1];
	    eig[1+2*i1]=a[i1+ndim*i ];
	  }
	  paired=T;
	} 
	else {
	  eig[2*i]=a[i+ndim*i];
	  eig[1+2*i]=0.;
	}
L2010: ;
   }
      //2010  continue

   if(!paired) {
	eig[2*(n-1)]=a[(n-1)+ndim*(n-1)];
	eig[1+2*(n-1)]=0.;
   }
   for(i=1;i<n;i++) {   //2110
	if(a[i+ndim*(i-1)] != 0.) goto L2110;

	i1=i+1;
	if(eig[1+2*i] == 0.) {
	  jm=0;
	  for(j=i-1;j>=0;j--) {   //  2120
	    if(a[j+1+ndim*j] != 0.) goto L2120;
	    if(eig[1+2*j] == 0.) {
	      double da=a[j+ndim*j]-a[i+ndim*i];
	      if(da != 0.) {
		if(fabs(da) > (fabs(a[i+ndim*i])+fabs(a[j+ndim*j]))*bconv) {
		  p=a[j+ndim*i]/da;
		  for(k=i;k<n;k++) 
		    a[j+ndim*k]=a[j+ndim*k]+p*a[i+ndim*k];
		  for(k=0;k<=j;k++) 
		    a[k+ndim*i]=a[k+ndim*i]-p*a[k+ndim*j];
		  for(k=0;k<n;k++) 
		    w[k+n*i]=w[k+n*i]-p*w[k+n*j];
		} 
		else {
		  a[j+ndim*j]=a[i+ndim*i];
		  eig[2*j]=eig[2*i];
		  jordan=T;
		}
	      } 
	      else {
		jordan=T;
	      }
	    } 
	    else {
	      j1=j-1;
	      u=a[i+ndim*i]-a[j1+ndim*j1];
	      v=a[j1+ndim*j];
	      r=u*u+v*v;
	      ck=Complex((a[j1+ndim*i]*u+a[j +ndim*i]*v)/r,
		  (a[j +ndim*i]*u-a[j1+ndim*i]*v)/r);
	      for(k=i;k<n;k++) {
		a[j1+ndim*k]=a[j1+ndim*k]-real(ck)*a[i+ndim*k];
		a[j +ndim*k]=a[j +ndim*k]-imag(ck)*a[i+ndim*k];
	      }
	      for(k=0;k<=j;k++)
		a[k+ndim*i]=a[k+ndim*i]+real(ck)*a[k+ndim*j1]+imag(ck)*a[k+ndim*j];
	      for(k=0;k<n;k++)
		w[k+n*i]=w[k+n*i]+real(ck)*w[k+n*j1]+imag(ck)*w[k+n*j];
	    }
L2120: ;
	  } //2120      continue
	} 
	else {
	  for(j=i-1;j>=0;j--) {    //2210
	    if(a[j+1+ndim*j] != 0.) goto L2210;
	    if(eig[1+2*j] == 0.) {
	      x=a[i+ndim*i]-a[j+ndim*j];
	      y=a[i+ndim*i1];
	      r=x*x+y*y;
	      ck=Complex((a[j+ndim*i ]*x+a[j+ndim*i1]*y)/r,
		(a[j+ndim*i1]*x-a[j+ndim*i ]*y)/r);
	      for(k=i;k<n;k++)
		a[j+ndim*k]=a[j+ndim*k]-real(ck)*a[i+ndim*k]-imag(ck)*a[i1+ndim*k];
	      for(k=0;k<=j;k++) {
		a[k+ndim*i ]=a[k+ndim*i ]+real(ck)*a[k+ndim*j];
		a[k+ndim*i1]=a[k+ndim*i1]+imag(ck)*a[k+ndim*j];
	      }
	      for(k=0;k<n;k++) {
		w[k+n*i ]=w[k+n*i ]+real(ck)*w[k+n*j];
		w[k+n*i1]=w[k+n*i1]+imag(ck)*w[k+n*j];
	      }
	    } 
	    else {
	      j1=j-1;
	      cu=Complex(a[j1+ndim*j1]-a[i+ndim*i],-a[i+ndim*i1]);
	      v =a[j1+ndim*j];
#if 0
	      cr=cu*cu+v*v;
#else
	      cr=cu*cu;
	      cr+=v*v;
#endif
	      if(cr != Complex(0.,0.)) {
		ca=Complex(a[j1+ndim*i],a[j1+ndim*i1]);
		cc=Complex(a[j +ndim*i],a[j +ndim*i1]);
		ck=(cu*ca-v*cc)/cr;
		cm=(cu*cc+v*ca)/cr;
		for(k=i;k<n;k++) {
		  a[j1+ndim*k]=a[j1+ndim*k]+real(ck)*a[i+ndim*k]+imag(ck)*a[i1+ndim*k];
		  a[j+ndim*k ]=a[j +ndim*k]+real(cm)*a[i+ndim*k]+imag(cm)*a[i1+ndim*k];
		}
		for(k=0;k<=j;k++) {
		  a[k+ndim*i ]=a[k+ndim*i ]-real(ck)*a[k+ndim*j1]-real(cm)*a[k+ndim*j];
		  a[k+ndim*i1]=a[k+ndim*i1]-imag(ck)*a[k+ndim*j1]-imag(cm)*a[k+ndim*j];
		}
		for(k=0;k<n;k++) {
		  w[k+n*i ]=w[k+n*i ]-real(ck)*w[k+n*j1]-real(cm)*w[k+n*j];
		  w[k+n*i1]=w[k+n*i1]-imag(ck)*w[k+n*j1]-imag(cm)*w[k+n*j];
		}
	      }
	    }
L2210:;
	  }  //2210      continue
	}
L2110:;
   }  // 2110  continue
   if(jordan) {
	for(i=n-1;i>=1;i--) { 
	  if(eig[1+2*i] == 0.) {
	    jm=-1;
	    for(j=i-1;j>=0;j--) {
	      if(eig[2*j] == eig[2*i] && eig[1+2*j] == 0.) {
		r=a[j+ndim*i];
		if(fabs(r) > 1.e-30) {
		  if(jm == -1) {
		    for(k=j;k<n;k++) a[j+ndim*k]=a[j+ndim*k]/r;
		    for(k=0;k<=j;k++) a[k+ndim*j]=a[k+ndim*j]*r;
		    for(k=0;k<n;k++) w[k+n*j]=w[k+n*j]*r;
		    jm=j;
		  } 
		  else {
		    for(k=jm;k<n;k++) a[j+ndim*k]=a[j+ndim*k]-r*a[jm+ndim*k];
		    for(k=0;k<=j;j++) a[k+ndim*jm]=a[k+ndim*jm]+r*a[k+ndim*j];
		    for(k=0;k<n;k++) w[k+n*jm]=w[k+n*jm]+r*w[k+n*j];
		  }
		} 
		else a[j+ndim*i]=0.;

	      }
	    }
	  }
	}
   }
//
//
// cout << "tqr check" << endl;
// awe(a,w,eig,n,ndim);
//
   if(n%2 == 0) {
	int ii=0;
	for(i=0;i<n;i++) {
	  if(eig[1+2*i] > 0. && ii != 0) {
	    for(j=0;j<n;j++) {
	      ww=w[j+n*(i-1)];
	      w[j+n*(i-1)]=w[j+n*i  ];
	      w[j+n*i  ]=w[j+n*(i+1)];
	      w[j+n*(i+1)]=ww;
	    }
	    ee=eig[2*(i-1)];
	    eig[2*(i-1)]=eig[2*i  ];
	    eig[2*i  ]=eig[2*(i+1)];
	    eig[2*(i+1)]=ee;
	    ee=eig[1+2*(i-1)];
	    eig[1+2*(i-1)]=eig[1+2*i  ];
	    eig[1+2*i  ]=eig[1+2*(i+1)];
	    eig[1+2*(i+1)]=ee;
	  }
	  ii=1-ii;
	}
   }
}
	
