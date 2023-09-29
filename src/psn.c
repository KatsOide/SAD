#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "psn.h"
#include "ft.h"

const double PI = 3.14159265358979;
const double TWO_PI = 6.28318530717958;

/* number of mesh */
int nh_count_psn;
int nv_count_psn;

double ratio_psn; 		/* cw/ch */


/* size: double[nhor_cell_count][nver_cell_count+1] */
double *g_facr_psn;

/* the following array are used to do dfst */
/* the max data length is 1024 */
int ip_dfst_psn[64];
double t_dfst_psn[1024];
double w_dfst_psn[1024];

/* the following two array are used to do rdft */
/* the max data length is 2048 */

void calc_g0_facr_psn(double* charge, double r, double* phi);
void calc_g1_facr_psn(double r);
void calc_evenusinft_facr_psn(double* phi, double r);
void calc_evenu_facr_psn(double* phi);
void calc_oddu_facr_psn(double* phi, double r);
void solve_3diag_facr_psn(double b, double ac, double* x, double* f, int n);

void create_env_psn(int nx, int ny, double cx, double cy)
{
     size_t size;
     double cw, ch;

     nh_count_psn = nx;
     nv_count_psn = ny;
     cw = cx;
     ch = cy;
     ratio_psn = cw/ch;

     ip_dfst_psn[0]=0;		/* */

     size = (nh_count_psn)*(nv_count_psn+1)*sizeof(double);
     g_facr_psn = (double*)malloc(size);
}

void clear_env_psn(void)
{
     free(g_facr_psn);
}


void calc_potential_psn(double *gridcharge, double *gridphi)
{
     static double u0_sinft[1024];
     static double uJ_sinft[1024];
     
     size_t size;
     int j, n;
     double r;
     double *g1, *g;

     size = (nh_count_psn+1)*(nv_count_psn+1)*sizeof(double);
     memset(gridphi, 0, size);
     /* Now potential on the boundary is zero, KOKOKO. */
     /*calc_sidepotential_psn(gridcharge, gridphi);*/

     /* facr poisson */
     
     r = ratio_psn;

     /* in the following, the cyclic reduction time is 1, i.e. k == 1 */
     calc_g0_facr_psn(gridcharge, r, gridphi);
     /* one time of cyclic reduction */
     calc_g1_facr_psn(r);

     /* sinft of g1 */
     g1 = g_facr_psn;
     for(j = 2; j < nh_count_psn; j+=2){
	  g1 += 2*(nv_count_psn+1);/* points to [j][0] */
	  g1[0] = 0;
	  dfst(nv_count_psn, g1, t_dfst_psn, ip_dfst_psn, w_dfst_psn);
     }

     /* sinft of u0, saved to u0_sinft */
     memcpy( u0_sinft, gridphi, sizeof(double)*nv_count_psn );
     u0_sinft[0] = 0;
     dfst(nv_count_psn,u0_sinft,t_dfst_psn,ip_dfst_psn,w_dfst_psn);

     /* sinft of uJ, saved to uJ_sinft */
     memcpy( uJ_sinft,
	     gridphi + nh_count_psn*(nv_count_psn+1),
	     sizeof(double)*nv_count_psn );
     uJ_sinft[0] = 0;
     dfst(nv_count_psn,uJ_sinft,t_dfst_psn,ip_dfst_psn,w_dfst_psn);

     /* g_facr[2][] -= u0_sinft[] */
     g = g_facr_psn + 2*(nv_count_psn+1);
     for( n = 1; n < nv_count_psn; n++)
	  g[n] -= u0_sinft[n];
     /* g_facr[J-2][] -= uJ_sinft[] */
     g = g_facr_psn + (nh_count_psn-2)*(nv_count_psn+1);
     for( n = 1; n < nv_count_psn; n++)
	  g[n] -= uJ_sinft[n];

     /* calc sinft of evenu */
     calc_evenusinft_facr_psn(gridphi, r);

     /* calc evenu by inverse-sinft */
     calc_evenu_facr_psn(gridphi);

     /* calc oddu by inverse-cyclic-reduction */
     calc_oddu_facr_psn(gridphi, r);

}

void calc_g0_facr_psn(double* charge, double r, double* phi)
{
     int j, l;
     double cc, cu;
     double *g0, *c, *u;

     cc = -TWO_PI*r;
     cu = r*r;

     g0 = g_facr_psn;
     c = charge; 
     u = phi;

     for(j = 1; j < nh_count_psn; j++){
	  g0 += nv_count_psn+1; /* g0 points to [j][0] */
	  c += nv_count_psn+1; /* c points to [j][0] */
	  u += nv_count_psn+1; /* u points to [j][0] */

	  g0[1] = c[1] * cc - u[0] * cu;
	  for(l = 2; l < nv_count_psn-1; l++)
	       g0[l] = c[l] * cc;
	  g0[l] = c[l] * cc - u[l+1] * cu;
     }
}

void calc_g1_facr_psn(double r)
{
     int j, l;
     double *g1, *left, *right;
     double temp1, temp2;
     double r2, diag;

     r2 = r*r;
     diag = -(1+r2)*2;
     g1 = g_facr_psn;

     for(j = 2; j < nh_count_psn; j+=2){
	  g1 += 2*(nv_count_psn+1); /* g1 points to [j][0] */
	  left = g1 - (nv_count_psn+1); 
	  right = g1 + (nv_count_psn+1);

	  temp1 = diag * g1[1] + r2 * g1[2];
	  for(l = 2; l < nv_count_psn-1; l++){
	       temp2 = r2*( g1[l-1] + g1[l+1] ) + diag*g1[l];
	       g1[l-1] = left[l-1] - temp1 + right[l-1];
	       temp1 = temp2;
	  }
	  temp2 = r2 * g1[l-1] + diag*g1[l];
	  g1[l-1] = left[l-1] - temp1 + right[l-1];
	  g1[l] = left[l] -temp2 + right[l];
     }
}

void calc_evenusinft_facr_psn(double* phi, double r)
{
     static double f[512];
     static double x[512];
     int i, j, n;
     double *g, *evenusinft;
     double b, r2;

     r2 = r*r;

     for(n = 1; n < nv_count_psn; n++){
	  i = 0;
	  g = g_facr_psn + n;/* points to [0][n] */
	  for(j = 2; j < nh_count_psn; j+=2){
	       g += 2*(nv_count_psn+1);/* points to [j][n] */
	       f[i++] = *g;
	  }

	  b = 2 * ( cos(n*PI/nv_count_psn) * r2 - r2 - 1 );
	  b *= -b;
	  b += 2;

	  solve_3diag_facr_psn(b, 1.0, x, f, nh_count_psn/2 - 1);

	  i = 0;
	  evenusinft = phi + n;/* points to [0][n] */
	  for(j = 2; j < nh_count_psn; j+=2){
	       evenusinft += 2*(nv_count_psn+1);/* points to [j][n] */
	       *evenusinft = x[i++];
	  }
     }
}

void calc_evenu_facr_psn(double* phi)
{
     int j,l;
     double u0;
     double invfactor;

     invfactor = 2.0/nv_count_psn;

     for(j = 2; j < nh_count_psn; j+=2){
	  phi += 2*(nv_count_psn+1);/* points to [rank*nhor_cell_count/4+2*j][0] */

	  u0 = phi[0];
	  phi[0] = 0;
	  dfst(nv_count_psn, phi, t_dfst_psn, ip_dfst_psn, w_dfst_psn);
	  phi[0] = u0;

	  for(l = 1; l < nv_count_psn; l++)
	       phi[l] *= invfactor;
     }
}

void calc_oddu_facr_psn(double* phi, double r)
{
     double *left, *right, *g0;
     double ac, b;
     int j,l;

     ac = r*r;
     b = -2*(1+ac);

     g0 = g_facr_psn - (nv_count_psn+1);/* points to [-1][0] */
     phi -= nv_count_psn + 1;/* points to [-1][0] */
     for( j = 1; j < nh_count_psn; j += 2 ){
	  g0 += 2 * (nv_count_psn+1);/* points to [j][0] */

	  phi += 2*(nv_count_psn+1);/* points to [j][0] */
	  left = phi - (nv_count_psn + 1);/* points to [j-1][0] */
	  right = phi + (nv_count_psn + 1);/* points to [j+1][0] */

	  for(l = 1; l < nv_count_psn; l++)
	       g0[l] -= left[l] + right[l];

	  solve_3diag_facr_psn(b, ac, phi+1, g0+1, nv_count_psn-1);
     }
}

void solve_3diag_facr_psn(double b, double ac, double* x, double* f, int n)
{
     static double u3diag[1024];
     int i;
     double l3diag;

     u3diag[0] = b;
     for(i = 1; i < n; i++){
	  l3diag = ac/u3diag[i-1];
	  u3diag[i] = b - l3diag*ac;
	  f[i] -= l3diag*f[i-1];
     }

     x[n-1] = f[n-1]/u3diag[n-1];
     for(i = n-2; i >= 0; i--)
	  x[i] = (f[i] - ac*x[i+1])/u3diag[i];
}

