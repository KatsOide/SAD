#ifndef _SBBE_MPI_FT
#define _SBBE_MPI_FT

/*
  cdft: Complex Discrete Fourier Transform
  rdft: Real Discrete Fourier Transform
  ddct: Discrete Cosine Transform
  ddst: Discrete Sine Transform
  dfct: Cosine Transform of RDFT (Real Symmetric DFT)
  dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
*/

void cdft(int, int, double *, int *, double *);
void rdft(int, int, double *, int *, double *);
void ddct(int, int, double *, int *, double *);
void ddst(int, int, double *, int *, double *);
void dfct(int, double *, double *, int *, double *);
void dfst(int, double *, double *, int *, double *);

/*
  -------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
  [definition]
  S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
  [usage]
  ip[0] = 0; // first time only
  dfst(n, a, t, ip, w);
  [parameters]
  n              :data length + 1 (int)
  n >= 2, n = power of 2
  a[0...n-1]     :input/output data (double *)
  output data
  a[k] = S[k], 0<k<n
  (a[0] is used for work area)
  t[0...n/2-1]   :work area (double *)
  ip[0...*]      :work area for bit reversal (int *)
  length of ip >= 2+sqrt(n/4)
  strictly, 
  length of ip >= 
  2+(1<<(int)(log(n/4+0.5)/log(2))/2).
  ip[0],ip[1] are pointers of the cos/sin table.
  w[0...n*5/8-1] :cos/sin table (double *)
  w[],ip[] are initialized if ip[0] == 0.
  [remark]
  Inverse of 
  dfst(n, a, t, ip, w);
  is 
  dfst(n, a, t, ip, w);
  for (j = 1; j <= n - 1; j++) {
  a[j] *= 2.0 / n;
  }
  .
*/


#endif
