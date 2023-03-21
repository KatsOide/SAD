/*************************************************************
*                                                            *
* Duplicated Differential Algebra Library at KEK.            *
*            C++ version 1            93/07/12               *
*                                     K.Ohmi                 *
*************************************************************/
#define Nresist 6
#include <iostream>
using std::ostream;
using std::cout;
#include <cstdlib>
using std::exit;
#include <cstring>
using std::memset;
#include <cmath>

#include <dacpp.h>
#include <p_da.h>
#define N_dim_max 3       // affect only asin table

void p_addr_tbl(void);

static int Nvar,Nord,Lvec;
static int *d1,*d2,*Id1,*NMS,*NM1,*ixref,*iv1,*iv2; 
static int Nvr2,Naddr,MaxC,Lvec2;

//
// copy constructor
//
/*da::da(const da& x)
{
   v=new double[Lvec];
   for(int i=0;i<Lvec;i++) v[i]=x.v[i];
}*/
p_da::p_da(const p_da& x){
   pv=new da[Lp_vec]; 
   for(int i=0;i<Lp_vec;i++) pv[i]=x.pv[i];
}
 

 
void p_da_init(int Npv,int Npo, const char* var_name[])
{
  int ialoc,i;
  
  cout << "\n\n Differential algebra of Duplicated variables\n\n";
  cout << "variables name : ";  
  for(i=0;i<Npv;i++) cout << var_name[i] << ' ';
  cout << '\n';
  
  Nord=Npo;
  Nvar=Npv;
  Nvr2=Nvar/2;
  Naddr=1;  for(i=0;i<Nvr2;i++) Naddr=Naddr*(Nord+1);
  MaxC=Nord; for(i=0;i<Nvr2-1;i++) MaxC=MaxC*(Nord+1);
  Lvec=icomb(Nord+Nvar,Nvar);
  Lvec2=icomb(Nord+Nvr2,Nvr2);
  Np_ord=Nord;
  Np_var=Nvar;
  Lp_vec=Lvec;

  ialoc=2*(Naddr+1)+Lvec2+2*(Nord+1)+2*Nvar;

  d1=new int[ialoc];
  cout << "*** " << ialoc << '*' << sizeof(int) 
	<< " Byte allocation for DA table \n";
  memset(d1,0,ialoc*sizeof(int));
/*	Make address table	*/
  p_addr_tbl();

}

void p_da_init(int Npv,int Npo)
{
  int ialoc,i;
  
  cout << "\n\n Differential algebra of Duplicated variables\n\n";
  cout << '\n';
  
  Nord=Npo;
  Nvar=Npv;
  Nvr2=Nvar/2;
  Naddr=1;  for(i=0;i<Nvr2;i++) Naddr=Naddr*(Nord+1);
  MaxC=Nord; for(i=0;i<Nvr2-1;i++) MaxC=MaxC*(Nord+1);
  Lvec=icomb(Nord+Nvar,Nvar);
  Lvec2=icomb(Nord+Nvr2,Nvr2);
  Np_ord=Nord;
  Np_var=Nvar;
  Lp_vec=Lvec;

  ialoc=2*(Naddr+1)+Lvec2+2*(Nord+1)+2*Nvar;

  d1=new int[ialoc];
  cout << "*** " << ialoc << '*' << sizeof(int) 
	<< " Byte allocation for DA table \n";
  memset(d1,0,ialoc*sizeof(int));
/*	Make address table	*/
  p_addr_tbl();

}
   

void p_addr_tbl(void)
{
/*
**************************************************************
*                                                            *
* Differential Algebra Library at KEK.                       *
*            sad version 1            92/07/19               *
*                                     K.Ohmi                 *
**************************************************************
*/

  int C1,D10,D20;
  int Nc,II,Nth,ICP;
  int i,j;

  d2=d1+Naddr+1;
  Id1=d2+Naddr+1;
  NMS=Id1+Lvec2;
  NM1=NMS+Nord+1;
  ixref=NM1+Nord+1;
  iv1=ixref+Nvar;
  iv2=iv1+Nvr2;

  //int* iv1=new int[Nvr2];

  cout << " \n \n Differential Algebra Library \n";
  cout <<" Order and number of variables = " 
	<< Nord << ' ' << Nvar << '\n';
  cout << " Total address is " << Naddr 
	<< "  :  Length of DA vector = " << Lvec << '\n';
  Nvr2=Nvar/2;

  for(i=0;i<Naddr;i++) {
    ICP=Naddr-1-i;
    ivcalp(ICP,Nord,Nvr2,iv1);

    Nc=0;
    for(j=0;j<Nvr2;j++) Nc+=*(iv1+j);

    if(Nc<=Nord) {

      C1=ccal(Nord,Nvr2,iv1);

      *(NM1+Nc)=(*(NM1+Nc))+1;
      if(C1<=MaxC) *(d1+C1)=*(NM1+Nc)-1;
/*      if(C1<=MaxC) *(d1+C1)=*(NM1+Nc); fortran version */
    }
  }

  *NMS=1;
  cout << "\n\n    Order      NM1       NMS  \n";
  for(i=0;i<Nord;i++) {
    *(NMS+i+1)=(*(NMS+i))+(*(NM1+i+1));
    cout.width(10); cout << i+1;
    cout.width(10); cout << *(NM1+i+1);
    cout.width(10); cout << *(NMS+i+1) << '\n';
  }


  if(Nord<4 && Nvar==4) 
    cout<<"\n        iv(1)     iv(2)       c1         d1        d2 \n";

  for(i=0;i<Naddr;i++) {
    ICP=Naddr-i-1;
    ivcalp(ICP,Nord,Nvr2,iv1);

    Nc=0;
    for(j=0;j<Nvr2;j++) Nc+=*(iv1+j);

    if(Nc<=Nord) {

      C1=ccal(Nord,Nvr2,iv1);

      D10=0;
      if(Nc>0) D10=*(NMS+Nc-1);
      D20=0;
      for(j=0;j<Nc;j++) D20+=*(NM1+j)**(NMS+Nord-j);

      if(C1<=MaxC) {

        *(d2+C1)=(*(d1+C1))**(NMS+Nord-Nc)+D20;
        *(d1+C1)=*(d1+C1)+D10;

/*c err FIXED 12/DEC/91 */

        II=*(d1+C1);

	if(II<=*(NMS+Nord)) {
	  *(Id1+II)=C1;
	}

	if(Nord<=4 && Nvar==4) 
    	  cout.width(10); cout << *iv1;
    	  cout.width(10); cout << *(iv1+1);
    	  cout.width(10); cout << C1;
    	  cout.width(10); cout << *(d1+C1);
    	  cout.width(10); cout << *(d2+C1) << '\n';
      }
    }
  }

  Nth=0;
  cout<< "\n\n      Order  No. till the order \n";
  for(i=0;i<=Nord;i++){
    for(j=0;j<=i;j++) Nth+=*(NM1+j)**(NM1+i-j);
    cout.width(10); cout << i;
    cout.width(10); cout << Nth << '\n';
  }
  cout<< "\n Dinv\n";
  for(i=0;i<Lvec2;i++) {
    cout.width(10); cout << i;
    cout.width(10); cout << *(Id1+i) << '\n';
  }
   int icx=1;
   for(i=0;i<Nvr2;i++) {
      if(i>0) icx=icx*(Nord+1); else icx=1;
      *(ixref+i)=*(d1+icx);
      *(ixref+i+Nvr2)=*d1+(*(d2+icx));
   }
   cout << "\n  1st derivative reference \n";
   for(i=0;i<Nvar;i++) {
    cout.width(10); cout << i;
    cout.width(10); cout << *(ixref+i) << '\n';
   }
}

p_da& p_da::operator=(const char* s)
{
   pv[0]=s;
   for(int i=1;i<Lvec;i++) pv[i]=0.;
   return *this;
}
//
//
//  *******   DA Library   ***************************
//
//

p_da operator*(const p_da& x, const p_da& y)
{
  p_da z0;
  da xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int No,N1A0;
  int is=0;

  z0=0.;
  No=Nord;
  if(x.pv[0]==0.) { is=1;}
  if(y.pv[0]==0.) { No=Nord-1; }

/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=No;N2A++) {
    if(is) { N1A0=1; } else {N1A0=0;}
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=N1A0;
      C2M=*(Id1+i2A+i2A0);
      for(i1A=is;i1A<*(NMS+No-N2A);i1A++){
	i=i1A+i1A0;
	xi=x.pv[i];
	C1M=*(Id1+i1A);

	if(i1A==*(NMS+N1A)) N1A++;
	NA=N2A+N1A;

	i1B0=0;  i2B0=0;
	for(N2B=0;N2B<=Nord-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+Nord-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      z0.pv[IMN]+=xi*y.pv[j];
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
    is=0;
  }
  return z0;
}


p_da muln(const p_da& y,const p_da& x,int No)
{
   p_da z0;
  da xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int is,N1A0;

   z0=0.;

/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=No;N2A++) {
    N1A0=No-N2A; 
    if(N1A0) is=*(NMS+N1A0-1); else is=0;

    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
       N1A=N1A0;
       C2M=*(Id1+i2A+i2A0);
       for(i1A=is;i1A<*(NMS+No-N2A);i1A++){
          i=i1A+i1A0;
          xi=x.pv[i];
          C1M=*(Id1+i1A);

	 if(i1A<*(NMS+N1A)) {	//N1A++;
	 NA=N2A+N1A;

	 i1B0=0;  i2B0=0;
	 for(N2B=0;N2B<=Nord-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+Nord-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      z0.pv[IMN]+=xi*y.pv[j];
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	 }
	 }
       }
       i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z0;
}


p_da da_mul(const p_da& x, const p_da& y, int No, int Mx, int My)
{
   p_da z0;
  da xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int Nox,is,N1A0;

  z0=0.;
  Nox=No-My;
  //  Noy=No-Mx;
/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nox;N2A++) {
    if(N2A>=Mx) { N1A0=0; is=0; }
    else {N1A0=Mx-N2A; is=*(NMS+N1A0-1);}

    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
       N1A=N1A0;
       C2M=*(Id1+i2A+i2A0);
       for(i1A=is;i1A<*(NMS+Nox-N2A);i1A++){
          i=i1A+i1A0;
          xi=x.pv[i];
          C1M=*(Id1+i1A);

	if(i1A==*(NMS+N1A)) N1A++;
	NA=N2A+N1A;

	i1B0=0;  i2B0=0;
	for(N2B=0;N2B<=No-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+No-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      z0.pv[IMN]+=xi*y.pv[j];
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z0;
}

//template <class T>
p_da operator/(double c,const p_da& x)
{
  p_da z1,z3,z;
  /*  if(x.pv[0]==0.) {
    cout <<" p_DA inverse is not usable because of da x(1)=0 \n";
    exit(1);
  }*/

  z3=x/x.pv[0];
  z3.ci(0,0.);

  //     z3  : Expanding variable

  z1=1.+z3;
  for(int i=2;i<=Nord;i++) {
     z1=1.+da_mul(z3,z1,i,1,0);
  }
  return ((c/x.pv[0])*z1);
}

p_da operator/(const p_da& x,const p_da& y)
{
   p_da z;
   z=x*(1./y);
   return z;
}

p_da pow(const p_da& x,double k)
{

  p_da z1,z3,z4;

  if(x.pv[0]==0.) {
    cout<< " p_DA pow is not usable because of x(0)=0 \n";
    exit(1);
  }

   // Bug ??? can not set break point hereafter but result is OK

  z3=x/(x.pv[0]);
  z3.pv[0]=0.;

  //     z3  : Expanding variable

  double f=(k-(double)(Nord-1))/(double)(Nord);

  if(N_ord>1) {
     z4=1.+z3*f;
     for(int i=1;i<(Nord-1);i++) {
       f=(k-(double)(Nord-i-1))/(double)(Nord-i);
       z4=da_mul(z3,z4,i+1,1,0);
       z4=1.+z4*f;
     }
  }
  else{ z4=1.; }
  z4=1.+z3*z4*k;

  da y1=pow(x.pv[0],k);
  return z4*y1;
}   


p_da sin(const p_da& x)
{
   p_da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne;
   
  da y1s=sin(x.pv[0]);
  da y1c=cos(x.pv[0]);
  if(x.pv[0]==0.) Mx=1;
  z1=x;
  z1.pv[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2=z1*z1;

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;
  }
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    f=-1./(double)(No*(No-1));
    zs=1.+z2*f;
    for(int i=2;i<No-2;i+=2) {
       f=-1./(double)((No-i)*(No-i-1));
       zs=da_mul(z2,zs,i,2,0); 
       zs=1.+zs*f;
    }
    zs=da_mul(z1,zs,Nord,1,0);
  //
  //cos like da
    if(!Mx) {
      f=-1./(double)(Ne*(Ne-1));
      zc=1.+z2*f;
      for(int i=2;i<=Ne-2;i+=2) {
        f=-1./(double)((Ne-i)*(Ne-i-1));
        zc=da_mul(z2,zc,i,2,0); 
        zc=1.+zc*f;
      }
    }
  }

  z=zs*y1c;
  if(!Mx) z+=zc*y1s;
  return z;
}



p_da cos(const p_da& x)
{
   p_da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne;
   
  da y1s=sin(x.pv[0]);
  da y1c=cos(x.pv[0]);
  if(x.pv[0]==0.) Mx=1;
  z1=x;
  z1.pv[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2=z1*z1;

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;
  }
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    if(!Mx){
    f=-1./(double)(No*(No-1));
    zs=1.+z2*f;
    for(int i=2;i<No-2;i+=2) {
       f=-1./(double)((No-i)*(No-i-1));
       zs=da_mul(z2,zs,i,2,0); 
       zs=1.+zs*f;
    }
    zs=da_mul(z1,zs,Nord,1,0);
    }
  //
  //cos like da
//    if(!Mx) {
      f=-1./(double)(Ne*(Ne-1));
      zc=1.+z2*f;
      for(int i=2;i<=Ne-2;i+=2) {
        f=-1./(double)((Ne-i)*(Ne-i-1));
        zc=da_mul(z2,zc,i,2,0); 
        zc=1.+zc*f;
      }
  //}
  }

  
  z=zc*y1c;
  if(!Mx) z-=zs*y1s;
  return z;
}


p_da sinh(const p_da& x)
{
  p_da z1,z2,zc,zs,z;
  double f;
  int Mx=0;
  int No,Ne;

  da y1s=sinh(x.pv[0]);
  da y1c=cosh(x.pv[0]);
  if(x.pv[0]==0.) Mx=1;
  z1=x;
  z1.pv[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2=z1*z1;

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;
  }
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    f=1./(double)(No*(No-1));
    zs=1.+z2*f;
    for(int i=2;i<No-2;i+=2) {
       f=1./(double)((No-i)*(No-i-1));
       zs=da_mul(z2,zs,i,2,0); 
       zs=1.+zs*f;
    }
    zs=da_mul(z1,zs,Nord,1,0);
  //
  //cos like da
    if(!Mx) {
      f=1./(double)(Ne*(Ne-1));
      zc=1.+z2*f;
      for(int i=2;i<=Ne-2;i+=2) {
        f=1./(double)((Ne-i)*(Ne-i-1));
        zc=da_mul(z2,zc,i,2,0); 
        zc=1.+zc*f;
      }
    }
  }

  z=zc*y1s;
  if(!Mx) z+=(zs*y1c);
  return z;
}

p_da cosh(const p_da& x)
{
   p_da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne;
   
  da y1s=sinh(x.pv[0]);
  da y1c=cosh(x.pv[0]);
  if(x.pv[0]==0.) Mx=1;
  z1=x;
  z1.pv[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2=z1*z1;

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;
  }
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    if(!Mx){
    f=1./(double)(No*(No-1));
    zs=1.+z2*f;
    for(int i=2;i<No-2;i+=2) {
       f=1./(double)((No-i)*(No-i-1));
       zs=da_mul(z2,zs,i,2,0); 
       zs=1.+zs*f;
    }
    zs=da_mul(z1,zs,Nord,1,0);
    }
  //
  //cos like da
      f=1./(double)(Ne*(Ne-1));
      zc=1.+z2*f;
      for(int i=2;i<=Ne-2;i+=2) {
        f=1./(double)((Ne-i)*(Ne-i-1));
        zc=da_mul(z2,zc,i,2,0); 
        zc=1.+zc*f;
      }
  //}
  }

  z=zc*y1c;
  if(!Mx) z+=(zs*y1s);
  return z;
}



p_da exp(const p_da& x)
{  
  p_da z3,z4; 

  da y1=exp(x.pv[0]);
  z3=x;
  z3.pv[0]=0.;
  
  //     z3  : Expanding variable

  z4=z3/((double)(Nord));
  z4.pv[0]=1.;
  for(int i=1;i<Nord;i++) {
     z4=da_mul(z3,z4,i+1,1,0);
     z4=1.+z4/((double)(Nord-i));
  }
  return z4*y1;
}




p_da log(const p_da& x)
{
   p_da z3,z4,z;
/*   if(x.v[0]<=0.) {
    cout <<" DA log is not usable because of x(1)<=0 \n";
    exit(1);
   }*/
   da y1=log(x.pv[0]);

   z3=x/(x.pv[0]);
   z3.pv[0]=0.;

  //     z3  : Expanding variable
  if(Nord>1){
     double f=-(double)(Nord-1)/(double)Nord;
     z4=1.+z3*f;
     for(int i=1;i<Nord-1;i++){
        f=-(double)(Nord-i-1)/(double)(Nord-i);
        z4=da_mul(z3,z4,i+1,1,0);
        z4=1.+z4*f;
     }
  }
  else {
     z4=1.;
  }
  z=z4*z3;
  z.pv[0]=y1;
  return z;
}

p_da sqrt(const p_da& x)
{
  p_da z1,z4,z3;
/*  if(x.v[0]<=0.) {
    cout <<" DA sqrt is not usable because of x(1)<=0 \n";
    exit(1);
  }*/

  da y1=sqrt(x.pv[0]);

  z3=x/(x.pv[0]);
  z3.pv[0]=0.;

  //     z3  : Expanding variable

  double f=-(double)(2*Nord-3)/(double)(2*Nord);

  if(Nord>1) {
     z4=1.+z3*f;
     for(int i=1;i<(Nord-1);i++) {
       f=-(double)(2*(Nord-i)-3)/(double)(2*(Nord-i));
       z4=da_mul(z3,z4,i+1,1,0);
       z4=1.+z4*f;
     }
  }
  else{ z4=1.; }
  z4=1.+z3*z4*0.5;
  return z4*y1;
}
//
//  Problem 94.1.24
//
p_da asin(const p_da& x)
{
  int i;
   da c[N_dim_max];
   
   // Make taylor coefficient table
   if(Nord>=N_dim_max) {
    cout<< " Nord is larger than max value (asin) "
	<< Nord << '>' << N_dim_max;
    exit(1);
  }
  da cr,cn;
  for(int n=1;n<=N_dim_max;n++) {
    cn=0.;
   double sign;
   for(int r=0;r<=n-1;r++) {
     if(r%2) sign=-1.; else sign=1.;
     int m=2*r-1;
     int k1=m;
     if(m<1) k1=1; else{
     while(m>=3){ 
        m=m-2;
        k1=k1*m;
     }}
     m=2*n-2*r-3;
     int k2=m;
     if(m<1) k2=1; else{
     while(m>=3){ 
        m=m-2;
        k2=k2*m;
     }}
     da x1=pow((1.+x.pv[0]),-0.5-r)*pow((1.-x.pv[0]),0.5-n+r);
     if(n!=1) 
         cr=sign*icomb(n-1,r)*k1*k2*x1;
     else cr=sign*k1*k2*x1;
     cn+=cr;
   }
   double c1=1.; for(i=0;i<n-1;i++) { c1=c1*0.5;}
   int nki=n;
   c[0]=asin(x.pv[0]);
   for(i=1;i<n;i++) {
      nki=nki*i; }
   c[n]=cn/(double)(nki)*c1;
   //cout << " n = " << n << " c[n]= " << cn '\n';
  }

  // main of asin da
  p_da z4;
  p_da z1=x;
  da f;
  
  z1.pv[0]=0.;
  
  if(x.pv[0]!=0.){
    if(Nord>1) {
     f=c[Nord]/c[Nord-1];
     // z4=1+c_n/c_{n-1} z1
     z4=1.+z1*f;
     for(i=1;i<(Nord-1);i++) {
       f=c[Nord-i]/c[Nord-i-1];
       z4=da_mul(z1,z4,i+1,1,0);
       z4=1.+z4*f;
     }
    }
    else{ z4=1.; }
    z4=c[0]+z1*z4*c[1];
  }
  else {
     int No;
    if(Nord%2) {No=Nord;} else {No=Nord-1;}
    if(Nord!=1) {
    p_da z2;
    z2=z1*z1;
    f=c[No]/c[No-2];
    z4=1.+z2*f;
    for(i=2;i<No-2;i+=2) {
       f=c[No-i]/c[No-i-2];
       z4=da_mul(z2,z4,i+1,1,0);
       z4=1.+z4*f;
    }
    z4=z1*z4;
    } else {z4=z1;}
    z4=z4*c[1];

  }
  return z4;
}

p_da atan(const p_da& x)
{
#if 0
  p_da z1,z4;
  da z2;
  double f;
  int No,i;

  z1=x;
  z1.pv[0]=0.;
  if(x.pv[0]!=0.) {
    z1=z1/(1.+x.pv[0]*x.pv[0]+x.pv[0]*z1);
  }

  if(Nord%2) {No=Nord;} else {No=Nord-1;}
  if(Nord!=1) {
    z2.mul(z1,z1);
    f=-double(No-2)/double(No);
    z4.mul(z2,f);
    *z4.pv=1.;
    for(i=2;i<No-2;i+=2) {
      f=-(double)(No-i-2)/(double)(No-i);
      z4=mul(z2,z4,i+1,1,0);
      z4.mul(z4,f);
      z4.pv[0]=1.;
    }
    z4.mul(z1,z4);
  } 
  else {z4=z1;}
  if(x.pv[0]!=0.) {
    z4+=atan(x.pv[0]);
  }
  return z4;
#else
  p_da z4;

  z4=x;
  return x;
#endif
}

p_da dif(const p_da& x,int k)
{
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A;
  p_da z;
  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	//if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  z.pv[i]=dif(x.pv[i],k);
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z;
}

p_da line_itg(const p_da& x)
{
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A;
  p_da z;
  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	//if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  z.pv[i]=line_itg(x.pv[i]);
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z;
}

p_da itg(const p_da& x,int k)
{
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A;
  p_da z;
  da y;
  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	//if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  z.pv[i]=itg(x.pv[i],k);
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z;
}


p_da p_dif(const p_da& x,int k)
{
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A;
  int kk,nodk,Iaddr;
  p_da z;

  if(k>=Nvar) {
    cout<< " k(" << k << ") >= Nvar(" << Nvar << ") :  terminated \n";
    exit(1);
  }

  z=0.;

  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);

      // ****** difinition of k  
      nodk=0;
      if(k>=Nvr2) {
	kk=k-Nvr2;
	nodk=*(iv2+kk);
	if(nodk>0) {
	  (*(iv2+kk))--;
	  C2M=ccal(Nord,Nvr2,iv2);
	}
      }

      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;

	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);

        if(k<Nvr2) {
          kk=k;
          nodk=*(iv1+kk);
          if(nodk>0) {
            (*(iv1+kk))--;
            C1M=ccal(Nord,Nvr2,iv1);
          }
	}

        if(nodk>0) {
          Iaddr=*(d1+C1M)+*(d2+C2M);
          z.pv[Iaddr]=x.pv[i]*nodk;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z;
}

p_da poi(const p_da& f,const p_da& g)
{
  p_da z1,z2,z4;
  z4=0.;
  for(int i=0;i<Nvr2;i++){
    z1=dif(f,i);
    z2=dif(g,i+Nvr2);
    z4+=(z1*z2);

    z1=dif(f,i+Nvr2);
    z2=dif(g,i);
    z4-=(z1*z2);
  }
  return z4;
}

ostream& operator<<(ostream& s,const p_da& x)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A;
  da xi;
  
  s << "\n\n  p_daprint Npord=" << Nord << " Npvar=" << Nvar 
       << " Lpvec=" << Lvec << "\n\n";
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=x.pv[i];
	if(i1A==*(NMS+N1A)) N1A++;
	//if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  for(j=0;j<Nvar;j++) {
	     s << iv1[j];
	  }
	  s << xi;
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  s << '\n';
  return s;
}
