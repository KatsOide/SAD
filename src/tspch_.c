#include <stdio.h>
#include <math.h>

#define pi2  6.28318530717958
extern void InitPotentialSolver(int nx,int ny,int nz,
                         double dx,double dy,double dz);
extern void ClearPotentialSolver(int ispend);
extern void SetCalcPotential(int it);
extern void UnSetCalcPotential(int it);
extern void UnSetCalcPotentialSolver(int it);
extern void SpaceChargeBPSet(double rnp,double Mass,double Gamma);
extern double **SpaceChargeParamSet(int nlatend,int isp);
//extern void SpaceChargeMapping(double** SP,int isp,double* x,int np,
//			       double TravelLength);
extern void SpaceChargeMapping(double** SP,int isp,double* x,double* px,
			       double* y, double* py, double* z,double* pz,
			       int np,double TravelLength);
extern void PotentialSolverPrint(void);
extern void epotkick(double* x,double* px,double* y,double* py,
                        double* z,double* pz,int np,
                        double TravelLength);


void tpspac_(int* nn, double* x, double *px, double* y, double* py,
	     double* z, double* g, double* dv, double* pz,
	     double* rnp_, double* Mass, double* P0_, 
	     double* Gamma, double* Length,
	     int* pspac_nx, int* pspac_ny, int* pspac_nz,
	     double* pspac_dx, double* pspac_dy, double* pspac_dz, 
	     int* pspac_nlat, int* pspac_nlatend,
	     int* pspac_nturnend, int* pspac_nturncalc)
{
  int i, np, nx, ny, nz;
  int nlat,nlatend,nturnend,nturncalc;
  double dx, dy, dz, rnp, Mp, TravelLength;
  double tpx,tpy,pr,p0,h0,h1;
  static double** SP;
  static int isp=0;
  static int ispend=0;
  static int nturn=1;

  np = *nn;
  rnp = *rnp_;
  Mp = *Mass; /* Unit eV */
  p0 = *P0_;
  h0 = *Gamma;
  TravelLength = *Length;

  nx = *pspac_nx;
  ny = *pspac_ny;
  nz = *pspac_nz;

  dx = *pspac_dx;
  dy = *pspac_dy;
  dz = *pspac_dz;

  nlat     = *pspac_nlat;
  nlatend  = *pspac_nlatend;
  nturnend = *pspac_nturnend;
  nturncalc= *pspac_nturncalc;

  if(nturncalc==0) nturncalc=nturnend;

  /* printf(" tpspac_ nturn %d nturnend %d nturncalc %d nlat %d nlatend %d\n",
     nturn,nturnend,nturncalc,nlat,nlatend);*/

  for(i=0;i<np;i++){
    tpx = px[i];
    tpy = py[i];
    pr = 1.0+g[i];

    px[i] = tpx*pr;
    py[i] = tpy*pr;
  }
  
  if(nturn==1){
    if(isp==0){
      InitPotentialSolver(nx,ny,nz,dx,dy,dz);

      SpaceChargeBPSet(rnp,Mp,h0);

      PotentialSolverPrint();

      SetCalcPotential(nturn);
    }
    SP=SpaceChargeParamSet(nlatend,isp);
  }
  else if(nturn>nturncalc){
    UnSetCalcPotential(nturn);
  }

  /*  SpaceChargeMapping(x,px,y,py,z,g,np,TravelLength);*/
  /*  printf(" tpspac_ SP %p SP[isp] %p isp %d\n",SP,SP[isp],isp);*/
  /*SpaceChargeMapping(SP,isp,x,np,TravelLength);*/
  SpaceChargeMapping(SP,isp,x,px,y,py,z,g,np,TravelLength);

  isp++;

  for(i=0;i<np;i++){
    tpx = px[i];
    tpy = py[i];
    pr = 1.0+g[i];
    h1 = sqrt(1.0+p0*p0*pr*pr);

    px[i] = tpx/pr;
    py[i] = tpy/pr;
    dv[i] = -g[i]*(1.0+pr)/h1/(h1+pr*h0);
  }    

  if(nturn==nturnend && isp==ispend){
    ClearPotentialSolver(ispend);
    printf("Freeing Memory at nturn %d nlat %d ispend %d\n",nturn,nlat,ispend);
  }

  if(nlat == nlatend){
    ispend = isp;
    isp = 0;
    nturn++;
  }
}

