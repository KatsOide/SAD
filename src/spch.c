#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>



/* #define m_e 0.51099906e-3  */
#define m_e 0.51099906e6 /* unit: eV */
#define r_e 2.81794092e-15 

extern void create_env_psn(int, int, double, double);
extern void clear_env_psn(void);
extern void calc_potential_psn(double*, double*);


void rhocal(double*,double*,double*,int);
void makebcu(double*);
void epotkick(double*,double*,double*,double*,double*,
	      double*,int,double);
void splined(double x[],double y[],int n,double yp1,double ypn,double y2[]);








double *rhoxy,*rhoz;
double *rhoz_sp,*z_sp;
double *btab,*ftab;
double *phi;
double **SP;

int nx_psn;
int ny_psn;
int nz_psn;

int np_psn;

int CalcPotentialFlag=1;

double dx_psn,dy_psn,dz_psn;

int ixc_psn,jyc_psn,kzc_psn;
double xc_psn,yc_psn,zc_psn,xm_psn,ym_psn,zm_psn;

double Kcoeff,r_0;
  
  
  
  
  


void InitPotentialSolver(int nx,int ny,int nz,
				 double dx,double dy,double dz)
{
  int size,i;

  nx_psn=nx;
  ny_psn=ny;
  nz_psn=nz;
  dx_psn=dx;
  dy_psn=dy;
  dz_psn=dz;

  printf(" InitPotentialSolver nx= %d \n",nx);

  size = (nx_psn+1)*(ny_psn+1)*sizeof(double);
  rhoxy = (double*)malloc(size);
  phi= (double*)malloc(size);
  /*phixy= (double*)malloc(size);*/
  size = (nz_psn+1)*sizeof(double);
  rhoz = (double*)malloc(size);
  rhoz_sp = (double*)malloc(size);
  z_sp=(double*)malloc(size);
  size = 16*(nx_psn+1)*(ny_psn+1)*sizeof(double);
  btab =  (double*)malloc(size);
  ftab = (double*)malloc(size);


  ixc_psn=nx_psn/2;
  jyc_psn=ny_psn/2;
  kzc_psn=nz_psn/2;
  xc_psn=ixc_psn*dx_psn;
  yc_psn=jyc_psn*dy_psn;
  zc_psn=kzc_psn*dz_psn;
  xm_psn=dx_psn*(ixc_psn-1);
  ym_psn=dy_psn*(jyc_psn-1);
  zm_psn=dz_psn*(kzc_psn-1);

  create_env_psn(nx_psn, ny_psn, dx_psn, dy_psn);

  for(i=0;i<nz_psn+1;i++) { 
    z_sp[i]=(i-kzc_psn)*dz_psn;
  }

  printf("Potential solver is initialized\n");
}


void ClearPotentialSolver(int ispend)
{
  int i;

  free(rhoxy);
  free(phi);
  free(rhoz);
  free(rhoz_sp);
  free(z_sp);
  free(btab);
  free(ftab);
  for(i=0;i<ispend-1;i++) free(SP[i]);
  free(SP);
  clear_env_psn();
}


void SpaceChargeBPSet(double rnp,double mass,double Gamma)
{
  double Beta;
  r_0=r_e*m_e/mass;
  Beta=sqrt(1.-1./Gamma/Gamma);
  Kcoeff=2.*rnp*r_0/Beta/Beta/Gamma/Gamma/Gamma;
  printf("Np, beta, gamma, K  %e %e %e %e \n",rnp,Beta,Gamma,Kcoeff);
}

double **SpaceChargeParamSet(int nlatend,int isp)
{
  int n;
  int size;

  if(isp==0){
    /*printf(" SpaceChargeParamSet nlatend= %d \n",nlatend);*/

    SP=(double **)malloc(nlatend*sizeof(double *));
    if(SP==NULL){
      printf("failed in malloc SP\n");
      exit(1);
    }
  }

  size=sizeof(double)*16*(nx_psn+1)*(ny_psn+1);

  SP[isp]=(double*)malloc(size);
  if(SP[isp]==NULL){
    printf("failed in malloc SP[isp]\n");
    exit(1);
  }
  memset(SP[isp], 0, size);

  /*  printf(" SpaceChargeParamSet isp= %d, SP = %p %p %p \n",isp,SP,
      SP[0],SP[isp]);*/
  return SP;
}

void PotentialSolverPrint()
{
  printf("   Nx Ny Nz = %d x %d x %d\n",nx_psn,ny_psn,nz_psn);
  printf("    dx dy dz = %f x %f x %f\n",dx_psn,dy_psn,dz_psn);
  printf(" r_0 = %le    2 N_p r_0/beta^2/gamm^3 = %le\n",r_0,Kcoeff);
}


//void SpaceChargeMapping(double* x,double* px,double* y,double* py,
//			double* z,double* pz,int np,
//			double TravelLength)

//void SpaceChargeMapping(double** SP,int isp,double* x,int np,
//			double TravelLength)
void SpaceChargeMapping(double** SP,int isp,double* x,double* px,
			double* y,double* py, double* z, double* pz,
			int np,double TravelLength)
{
  int i,j;
  /* double *px,*y,*py,*z,*pz;*/
  /* FILE *fout; */
  /* px=x+np; y=px+np; py=y+np; z=py+np; pz=z+np; */

  btab=SP[isp];

  /* printf(" SpaceChargeMapping isp= %d btab = %p TravelL = %f \n",
     isp,btab,TravelLength); */
  /* printf(" SpaceChargeMapping x= %p px = %p z = %p pz = %p \n",
     x,px,z,pz); */

  if(CalcPotentialFlag) {
    rhocal(x,y,z,np);

    np_psn=np;

    calc_potential_psn(rhoxy,phi);

    /* if(isp==0){
      fout=fopen("rhoxy.dat","w");
      for(i=0; i<nx_psn+1; i++){
	for( j = 0; j<ny_psn+1; j++ ){
	  fprintf(fout,"%d %d %lf \n", i, j, rhoxy[i*(ny_psn+1)+j]);
	}
      }
      fclose(fout);
      fout=fopen("phixy.dat","w");
      for(i=0; i<nx_psn+1; i++){
	for( j = 0; j<ny_psn+1; j++ ){
	  fprintf(fout,"%d %d %lf \n", i, j, phi[i*(ny_psn+1)+j]);
	}
      }
      fclose(fout);
      } */

    /*memcpy(phi,phixy,sizeof(double)*(nx_psn+1)*(ny_psn+1));*/

    makebcu(btab);
  }
  /*else {
    memcpy(phixy,phi,sizeof(double)*(nx_psn+1)*(ny_psn+1));
  }*/

  /*  makebcu(phi);*/

  splined(z_sp,rhoz,nz_psn,0.,0.,rhoz_sp);

  epotkick(x,px,y,py,z,pz,np,TravelLength);
}

void SetCalcPotential(int it)
{
  CalcPotentialFlag=1;
}

void UnSetCalcPotential(int it)
{
  CalcPotentialFlag=0;
}

void PrintPot(FILE *fdo)
{
  int i,j;

  for(i=0; i<nx_psn+1; i++){
    for( j = 0; j<ny_psn+1; j++ ){
      fprintf(fdo,"%d %d %lf %le\n", i, j, rhoxy[i*(ny_psn+1)+j],
	     phi[i*(ny_psn+1)+j]);
    }
  }
}

void rhocal(double* x,double* y,double* z,int np)
{

  int i,j,ii,jj,nm,iijj;
  double xic,yjc,zic,dxic,dyjc,dzic,rxi,ryj,rzi,rxin,ryjn,rzin;

  memset(rhoxy, 0, sizeof(double)*(nx_psn+1)*(ny_psn+1));
  memset(rhoz, 0, sizeof(double)*(nz_psn+1));

  nm=0;
  

  for(i=0;i<np;i++) {
    if(fabs(x[i])>xm_psn || fabs(y[i])>ym_psn) continue; 
    ii=(x[i]+xc_psn)/dx_psn;
    jj=(y[i]+yc_psn)/dy_psn;

    xic=(ii-ixc_psn)*dx_psn;
    yjc=(jj-jyc_psn)*dy_psn;
    dxic=x[i]-xic;
    dyjc=y[i]-yjc;
    rxin=fabs(dxic/dx_psn);
    rxi=1.-rxin;
    ryjn=fabs(dyjc/dy_psn);
    ryj=1.-ryjn;

    iijj=ii*(ny_psn+1)+jj;
    rhoxy[iijj]+=rxi*ryj;
    rhoxy[iijj+ny_psn+1]+=rxin*ryj;
    rhoxy[iijj+ny_psn+2]+=rxin*ryjn;
    rhoxy[iijj+1]+=rxi*ryjn;

    nm=nm+1;
    
    if(fabs(z[i])>zm_psn) continue; 
    ii=(z[i]+zc_psn)/dz_psn;

    /*
    zic=(ii-kzc_psn)*dz_psn;
    dzic=z[i]-zic;
    rzin=fabs(dzic/dz_psn);
    rzi=1.-rzin;
    if(ii<=0 || ii>nz_psn) continue;
    rhoz[ii]+=rzi;
    rhoz[ii+1]+=rzin;
    */

    zic=(ii-kzc_psn)*dz_psn;
    dzic=(z[i]-zic)/dz_psn;

    //if(ii<=0 || ii>nz_psn) continue;
    if(dzic<0.5) {
      rzin=2.*dzic*dzic;
      rzi=1.-rzin;
      rhoz[ii]+=rzi;
      rhoz[ii]+=rzin;
    }
    else {
      rzi=2.*(1.-dzic)*(1.-dzic);
      rzin=1.-rzi;
      rhoz[ii]+=rzi;
      rhoz[ii+1]+=rzin;
    }
    //printf("%d %f %f dzic=%f %f %f\n",ii,z[i],zic,dzic,rzi,rzin);
  }

}

  static double wt[16][16]=
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
      2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
      0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
     -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
      9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
     -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
      2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
     -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
      4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1 };

void makebcu(double* btab)
{
  double ax,ay,b;
  int i,j,k,m,l,mx,my,ij,ijf,size;
  double dx4,dy4,d1d2,xx;

  ax=1./(2.*dx_psn);
  ay=1./(2.*dy_psn);
  b=ax*ay;

  for(i=1;i<nx_psn;i++) {
    for(j=1;j<ny_psn;j++) {
      ij=i*(ny_psn+1)+j;
      ijf=ij*16;
      
      ftab[ijf]=phi[ij];
      ftab[1+ijf]=phi[ij+ny_psn+1];
      ftab[2+ijf]=phi[ij+ny_psn+2];
      ftab[3+ijf]=phi[ij+1];

      ftab[4+ijf]=(phi[ij+ny_psn+1]-phi[ij-ny_psn-1])*0.5;
      ftab[5+ijf]=(phi[ij+2*ny_psn+2]-phi[ij])*0.5;
      ftab[6+ijf]=(phi[ij+2*ny_psn+3]-phi[ij+1])*0.5;
      ftab[7+ijf]=(phi[ij+ny_psn+2]-phi[ij-ny_psn])*0.5;

      ftab[8+ijf]=(phi[ij+1]-phi[ij-1])*0.5;
      ftab[9+ijf]=(phi[ij+ny_psn+2]-phi[ij+ny_psn])*0.5;
      ftab[10+ijf]=(phi[ij+ny_psn+3]-phi[ij+ny_psn+1])*0.5;
      ftab[11+ijf]=(phi[ij+2]-phi[ij])*0.5;
      
      ftab[12+ijf]=(phi[ij+ny_psn+2]-phi[ij+ny_psn]-
		 phi[ij-ny_psn]+phi[ij-ny_psn-2])*0.25;
      ftab[13+ijf]=(phi[ij+2*ny_psn+3]-phi[ij+2*ny_psn+1]-
		 phi[ij+1]+phi[ij-1])*0.25;
      ftab[14+ijf]=(phi[ij+2*ny_psn+4]-phi[ij+2*ny_psn+2]-
		 phi[ij+2]+phi[ij])*0.25;
      ftab[15+ijf]=(phi[ij+ny_psn+3]-phi[ij+ny_psn+1]-
		 phi[ij-ny_psn+1]+phi[ij-ny_psn-1])*0.25;
	  
    }
  }


      
  for(m=0;m<16;m++) {
    for(i=1;i<nx_psn;i++) {
      for(j=1;j<ny_psn;j++) {
	ij=i*(ny_psn+1)+j;
	ijf=ij*16;
	btab[m+ijf]=wt[m][0]*ftab[ijf]+wt[m][1]*ftab[1+ijf]+
             	    wt[m][2]*ftab[2+ijf]+wt[m][3]*ftab[3+ijf]+
	            wt[m][4]*ftab[4+ijf]+wt[m][5]*ftab[5+ijf]+
	            wt[m][6]*ftab[6+ijf]+wt[m][7]*ftab[7+ijf]+
	            wt[m][8]*ftab[8+ijf]+wt[m][9]*ftab[9+ijf]+
	            wt[m][10]*ftab[10+ijf]+wt[m][11]*ftab[11+ijf]+
	            wt[m][12]*ftab[12+ijf]+wt[m][13]*ftab[13+ijf]+
	            wt[m][14]*ftab[14+ijf]+wt[m][15]*ftab[15+ijf];
      }
    }
  }
  
  /*  for(i=1;i<nx_psn;i++) {
    j=64;
    for(m=0;m<16;m++) {
	ij=i*(ny_psn+1)+j;
	ijf=ij*16;
	printf("%d %d %d %f %f\n",i,j,m,btab[m+ijf],f[m+ijf]);
    }
  }
  printf("Compare wt\n");
    for(i=0;i<16;i++) {
      for(j=0;j<16;j++) {
	if(wtfor[i][j]!=wt[i][j]) {
	  printf("%d %d %f %f\n",i,j,wtfor[i][j],wt[i][j]);
	}
      }
    }
    printf("Compare wt end\n");*/

}




void epotkick(double* x,double* px,double* y,double* py,double* z,
	      double* pz,int np,double TravelLength)
{
  double xu,xl,yu,yl,zu,zl,t,u,v,ph,phx,phy,a4,a3,a2,a1,w;
  double lambda,dlambda_dz,Kfac,lambda0,dlambda_dz0;
  int i,ii,jj,kk,nn,it,ij,ijf;

  Kfac=Kcoeff*TravelLength/np_psn/dz_psn/np_psn;

  for(i=0;i<np;i++) {
    if(fabs(x[i])>xm_psn || fabs(y[i])>ym_psn) continue;

    ii=(x[i]+xc_psn)/dx_psn;
    jj=(y[i]+yc_psn)/dy_psn;
    kk=(z[i]+zc_psn)/dz_psn;
    
    if(ii<1 || ii>nx_psn-1) continue;
    if(jj<1 || jj>ny_psn-1) continue;
    xl=(ii-ixc_psn)*dx_psn;
    xu=xl+dx_psn;
    yl=(jj-jyc_psn)*dy_psn;
    yu=yl+dy_psn;
    zl=(kk-kzc_psn)*dz_psn;
    zu=zl+dz_psn;
    
    t=(x[i]-xl)/dx_psn;
    u=(y[i]-yl)/dy_psn;
    v=(z[i]-zl)/dz_psn;

    ph=0.;
    phx=0.;
    phy=0.;
    ij=ii*(ny_psn+1)+jj;
    ijf=ij*16;
    a4=((btab[15+ijf]*u+btab[14+ijf])*u+ 
	 btab[13+ijf])*u+btab[12+ijf];
    a3=((btab[11+ijf]*u+btab[10+ijf])*u+ 
	 btab[9+ijf])*u+btab[8+ijf];
    a2=((btab[7+ijf]*u +btab[6+ijf])*u+ 
	 btab[5+ijf])*u+btab[4+ijf];
    a1=((btab[3+ijf]*u +btab[2+ijf])*u+ 
	 btab[1+ijf])*u+btab[ijf];
    ph=a1+(a2+(a3+a4*t)*t)*t;

    a4=(3.*btab[15+ijf]*u+2.*btab[14+ijf])*u+btab[13+ijf];
    a3=(3.*btab[11+ijf]*u+2.*btab[10+ijf])*u+btab[9+ijf];
    a2=(3.*btab[7+ijf]*u+2.*btab[6+ijf])*u+btab[5+ijf];
    a1=(3.*btab[3+ijf]*u+2.*btab[2+ijf])*u+btab[1+ijf];
    phy=a1+(a2+(a3+a4*t)*t)*t;

    a4=(3.*btab[15+ijf]*t+2.*btab[11+ijf])*t+btab[7+ijf];
    a3=(3.*btab[14+ijf]*t+2.*btab[10+ijf])*t+btab[6+ijf];
    a2=(3.*btab[13+ijf]*t+2.*btab[9+ijf])*t+btab[5+ijf];
    a1=(3.*btab[12+ijf]*t+2.*btab[8+ijf])*t+btab[4+ijf];
    phx=a1+(a2+(a3+a4*u)*u)*u;

    phx=phx/dx_psn;
    phy=phy/dy_psn;

    if(kk<0 || kk>=nz_psn) {
      printf("Out of z slice\n");
      continue;
    }
    // The latest rhoz is used. Is it okay?
    lambda0=(rhoz[kk+1]*v+rhoz[kk]*(1.-v));  // /dz_psn/P.N(); Kfac
    dlambda_dz0=(rhoz[kk+1]-rhoz[kk])/dz_psn;       // /dz_psn/P.N();

    /*  Spline fit */
    w=1.-v;
    lambda=w*rhoz[kk]+v*rhoz[kk+1]+
      ((w*w*w-w)*rhoz_sp[kk]+(v*v*v-v)*rhoz_sp[kk+1])*(dz_psn*dz_psn)/6.;
    dlambda_dz=(rhoz[kk+1]-rhoz[kk])/dz_psn
      -((3.*w*w-1.)*rhoz_sp[kk]-(3.*v*v-1.)*rhoz_sp[kk+1])*dz_psn/6.;
    //printf("%f %e %e %e %e\n",z[i],lambda0,lambda,dlambda_dz0,dlambda_dz);
    
    px[i]-=Kfac*lambda*phx;
    py[i]-=Kfac*lambda*phy;
    
    pz[i]-=Kfac*dlambda_dz*ph;
    
    
  }
}

void splined(double x[],double y[],int n,double yp1,double ypn,double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;
  u=(double*)malloc(sizeof(double)*n);
  if(yp1>.99e33) {
    y2[0]=0.;
    u[0]=0.;}
  else {
    y2[0]=-0.5;
    u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for(i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.;
    y2[i]=(sig-1.0)/p;
    u[i]=(6.*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))
          /(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if(ypn>.99e30) {
    qn=0.;
    un=0.; }
  else {
    qn=0.5;
    un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.);
  for(k=n-2;k>=0;k--) {
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
  free(u);
}




	 
