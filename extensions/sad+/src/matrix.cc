#include <iostream>
using std::ostream;
using std::ios;
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>
using std::fabs;
using std::sqrt;
using std::sin;
using std::cos;
using std::sinh;
using std::cosh;
using std::exp;

#include <matrix.h>

void teigen(double*,double*,int ,int );

matrix::matrix(const matrix& x) 
{
   Nr=x.Nr;
   Nc=x.Nc;
   int n=Nr*Nc;
   t=new double [n];
   for(int i=0;i<n;i++) t[i]=x.t[i];
}

matrix& matrix::operator=(const matrix& x) 
{
   Nr=x.Nr;
   Nc=x.Nc;
   int n=Nr*Nc;
   for(int i=0;i<n;i++) t[i]=x.t[i];
   return *this;
}
matrix& matrix::operator=(double f) 
{
  int i;
  int n=Nr*Nc;
  for(i=0;i<n;i++) t[i]=0.;
  if(Nr!=Nc && f!=0.) {
    cout << "Error: N_row != N_column : return Null Matrix\n";
    return *this;
  }
  for(i=0;i<Nr;i++) t[i*Nc+i]=f;
  return *this;
}

matrix matrix::operator-(void)
{
   matrix x(Nr,Nc);
   for(int i=0;i<Nr*Nc;i++) x.t[i]=-t[i];
   return x;
}


ostream& operator<<(ostream& s, matrix& A)
{
   for(int i=0;i<A.Nr;i++) {
      for(int j=0;j<A.Nc;j++) {
        s.width(12);
	s << A.t[i*A.Nc+j] << ' ';}
      s << '\n';
   }
   return s;
}

void PrintList(const matrix& x)
{
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(15);
  cout << '{';
  for(int i=0;i<x.Nr;i++) {
    cout << '{';
    for(int j=0;j<x.Nc;j++) {
      if(j!=x.Nc-1) {
	cout << x.t[i*x.Nc+j] << ','; }
      else {
	cout << x.t[i*x.Nc+j] << '}'; }
      if(j==2) cout <<'\n';
    }
    if(i<x.Nr-1) cout << ",\n";
  }
  cout << "};\n";
  cout.unsetf(ios::scientific);
  cout.unsetf(ios::floatfield);
  cout.precision(6);
}


void teigen(const matrix& x,matrix& R,double* eig)
{
  if(x.Nr!=R.Nr || x.Nc!=R.Nc) {
    cout << " x and  have different dimensions\n";
    exit(1);
  }

  R=Transpose(x);
  teigen(R.t,eig,x.Nr,x.Nc);
  R=Transpose(R);
}

matrix diag(const matrix& x,double* eig)
{
   cout << "\n\n(teigen.f) Eigen value routine by Oide (fortran)\n";
   int Nr=x.Nr;
   int Nc=x.Nc;
   matrix R(Nr,Nc);

   int i;

   R=Transpose(x);
 
   teigen(R.t,eig,Nr,Nc);
   
   cout << "\nEigen vectors\n" << R;
   R=Transpose(R);
   
   cout << "\n Eigen values: Real and Imag part\n";
   for(i=0;i<Nr;i++) {
      cout.width(10);
      cout << eig[2*i];} 
   cout << '\n';
   for(i=0;i<Nr;i++) {
      cout.width(10);
      cout << eig[i*2+1]; } 
   cout << '\n';
   cout.flush();
   return R;
}

matrix operator+(const matrix& x,const matrix& y)
{
  if(x.Nc!=y.Nc || x.Nr!=y.Nr) {
    cout << "Error in matrix+matrix\n"; exit(1);
  }
  matrix z(x.Nr,x.Nc);
  for(int i=0;i<x.Nc*x.Nr;i++) z.t[i]=x.t[i]+y.t[i];
  return z;
}

matrix operator-(const matrix& x,const matrix& y)
{
  if(x.Nc!=y.Nc || x.Nr!=y.Nr) {
    cout << "Error in matrix+matrix\n"; exit(1);
  }
  matrix z(x.Nr,x.Nc);
  for(int i=0;i<x.Nc*x.Nr;i++) z.t[i]=x.t[i]-y.t[i];
  return z;
}


matrix operator*(const matrix& x,const matrix& y)
{
   if(x.Nc!=y.Nr) { 
      cout << " Error : cannot multiply the matrices \n";
      exit(1);
   }
   matrix z(x.Nr,y.Nc);
   z=0.;
   int i,j,k;
   
   for(i=0;i<x.Nr;i++) {
      for(j=0;j<y.Nc;j++) {
	 for(k=0;k<x.Nc;k++){
	    z.t[i*x.Nc+j]+=x.t[i*x.Nc+k]*y.t[k*y.Nc+j];
	 }
      }
   }
   return z;
}

matrix operator*(double a,const matrix& x)
{
   matrix z(x.Nr,x.Nc);
   for(int i=0;i<x.Nr*x.Nc;i++) z.t[i]=a*x.t[i];
   return z;
}

matrix operator*(const matrix& x,double a)
{
   matrix z(x.Nr,x.Nc);
   for(int i=0;i<x.Nr*x.Nc;i++) z.t[i]=a*x.t[i];
   return z;
}

matrix operator/(const matrix& x,double a)
{
   matrix z(x.Nr,x.Nc);
   for(int i=0;i<x.Nr*x.Nc;i++) z.t[i]=x.t[i]/a;
   return z;
}

void matrix::Symp(void)
{
  int i;
  if(Nr!=Nc) { cout << "Error : Nr!=Nc \n"; exit(1); }
  for(i=0;i<Nc*Nc;i++) t[i]=0.;
  for(i=0;i<Nr;i+=2){
    t[i*Nc+i+1]=1.;
    t[(i+1)*Nc+i]=-1.;
  }
}

matrix SymplecticMatrix(int N)
{
  int i;
  matrix x(N,N);
  if(N%2) { cout << "Error : N%2!=0 \n"; exit(1); }
  for(i=0;i<N*N;i++) x.t[i]=0.;
  for(i=0;i<N;i+=2){
    x.t[i*N+i+1]=1.;
    x.t[(i+1)*N+i]=-1.;
  }
  return x;
}

matrix Transpose(const matrix& x)
{
   matrix y(x.Nc,x.Nr);
   for(int i=0;i<x.Nr;i++) {
     for(int j=0;j<x.Nc;j++) y.t[i+x.Nr*j]=x.t[i*x.Nc+j];
   }
   return y;
}

void swap(double&,double&);

inline void swap(double& a,double& b)
{
   double tmp=a;
   a=b;
   b=tmp;
}

matrix Inverse_G(const matrix& x)
{
  int i,j,k;
  if(x.Nr!=x.Nc) { cout << "Error : Nr!=Nc \n"; exit(1); }
  matrix a=x,b(x.Nr,x.Nc);
  b=1.;
  for(k=0; k<x.Nc-1;k++) {
    int mx=k;
    for(i=k+1;i<x.Nc;i++) 
      if(fabs(a.t[i*x.Nc+k])>fabs(a.t[mx*x.Nc+k])) mx=i;
    if(mx!=k) {
      for(i=k;i<x.Nc;i++) swap(a.t[k*x.Nc+i],a.t[mx*x.Nc+i]);
      for(i=0;i<x.Nc;i++) swap(b.t[k*x.Nc+i],b.t[mx*x.Nc+i]);
    }
    for(i=k+1;i<x.Nc;i++) {
      double tmp=a.t[i*x.Nc+k]/a.t[k*x.Nc+k];
      for(j=0;j<x.Nc;j++) {
	a.t[i*x.Nc+j]-=tmp*a.t[k*x.Nc+j];
	b.t[i*x.Nc+j]-=tmp*b.t[k*x.Nc+j];
      }
    }
  }
  for(k=x.Nc-1;k>=0;k--) {
    for(i=0;i<x.Nc;i++) b.t[k*x.Nc+i]/=a.t[k*x.Nc+k];
    for(i=k-1;i>=0;i--) {
      for(j=0;j<x.Nc;j++) 
	b.t[i*x.Nc+j]-=b.t[k*x.Nc+j]*a.t[i*x.Nc+k];
    }
  }
  return b;
}

void is_symplectic(const matrix& V)
{
   if(V.Nr!=V.Nc) { cout << "Error : Nr!=Nc \n"; exit(1); }
   matrix S(V.Nr,V.Nc);
   matrix T(V.Nr,V.Nc);
   S.Symp();
   T=Transpose(V)*S*V;
   cout << "-----------------------------------------\n" <<
            "Symplecticity check : \n" << T;
}      

/* Numerical recipes in C */


#define TINY 1.e-20
void ludcmp(double*,int,int*,double*);
void lubksb(double*,int,int*,double*);


matrix Inverse(const matrix& x)
{
   int i,j;
   if(x.Nr!=x.Nc) { cout << "Error : Nr!=Nc \n"; exit(1); }
   matrix a=x,b(x.Nr,x.Nc);

   double d;
   int* indx=new int [a.Nr];
   double* col=new double [a.Nr];
   ludcmp(a.t,a.Nr,indx,&d);
   for(j=0;j<x.Nr;j++) {
      for(i=0;i<x.Nr;i++) col[i]=0.;
      col[j]=1.0;
      lubksb(a.t,a.Nr,indx,col);
      for(i=0;i<x.Nr;i++) b.t[x.Nr*i+j]=col[i];
   }
   delete [] col;
   delete [] indx;
   return b;
}

double Det(const matrix& x)
{
   int j,iz=0;
   if(x.Nr!=x.Nc) { cout << "Error : Nr!=Nc \n"; exit(1); }
   for(j=0;j<x.Nc*x.Nr;j++) {
     if(x.t[j]!=0.) iz=1;
   }
   if(iz==0) return 0.;
     
   matrix a=x;

   double d;
   int* indx=new int [x.Nr];
   ludcmp(a.t,a.Nr,indx,&d);
   for(j=0;j<a.Nr;j++) d*=a.t[j+x.Nr*j];
   delete [] indx;
   return d;
}   

void ludcmp(double *a,int n,int *indx,double *d)
{
   int i,imax,j,k;
   double big,dum,sum,temp;
   double *vv;
   
   vv=new double [n];
   *d=1.;
   for(i=0;i<n;i++) {
      big=0.;
      for(j=0;j<n;j++)
         if((temp=fabs(a[n*i+j]))>big) big=temp;
      if(big==0.) { 
          cout << "Singular matrix in routine ludcmp";
          exit(1);
      }
      vv[i]=1./big;
   }
   for(j=0;j<n;j++) {
      for(i=0;i<j;i++) {
         sum=a[n*i+j];
         for(k=0;k<i;k++) sum-=a[n*i+k]*a[n*k+j];
         a[n*i+j]=sum;
      }
      big=0.;
      for(i=j,imax=j;i<n;i++) {
         sum=a[n*i+j];
         for(k=0;k<j;k++) sum-=a[n*i+k]*a[n*k+j];
         a[n*i+j]=sum;
         if((dum=vv[i]*fabs(sum))>=big) {
            big=dum;
            imax=i;
         }
      }
      if(j!=imax) {
         for(k=0;k<n;k++) {
            dum=a[n*imax+k];
            a[n*imax+k]=a[n*j+k];
            a[n*j+k]=dum;
         }
         *d=-(*d);
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if(a[n*j+j]==0.) a[n*j+j]=TINY;
      if(j!=n) {
         dum=1./a[n*j+j];
         for(i=j+1;i<n;i++) a[n*i+j]*=dum;
      }
   }
   delete [] vv;
}

void lubksb(double *a,int n,int *indx,double *b)
{
   int i,ii=-1,ip,j;
   double sum;
   
   for(i=0;i<n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii!=-1) {
         for(j=ii;j<i;j++) sum-=a[n*i+j]*b[j];}
      else if(sum) ii=i;
      b[i]=sum;
   }
   for(i=n-1;i>=0;i--) {
      sum=b[i];
      for(j=i+1;j<n;j++) sum-=a[n*i+j]*b[j];
      b[i]=sum/a[n*i+i];
   }
}

matrix SInverse(const matrix& x)
{
   matrix y(x.Nr,x.Nc),S(x.Nr,x.Nc);
   S.Symp();
   y=-S*Transpose(x)*S;
   return y;
}

void matrix::DiagonalMatrix(double* a)
{
  int i;
   for(i=0;i<Nr*Nc;i++) t[i]=0.;
   int n=(Nr<Nc) ? Nr : Nc;
   for(i=0;i<n;i++) t[i+Nc*i]=a[i];
}

void matrix::I(void)
{
  int i;
   for(i=0;i<Nr*Nc;i++) t[i]=0.;
   int n=(Nr<Nc) ? Nr : Nc;
   for(i=0;i<n;i++) t[i+Nc*i]=1.;
}

void matrix::Clear(void)
{
   for(int i=0;i<Nr*Nc;i++) t[i]=0.;
}

matrix IdentityMatrix(int N)
{
  int i;
  matrix x(N,N);
  for(i=0;i<N*N;i++) x.t[i]=0.;
  for(i=0;i<N;i++) x.t[i+N*i]=1.;
  return x;
}

matrix Append_c(const matrix& x,const matrix& y)
{
   if(x.Nr!=y.Nr) {
      cout << "x.Nr!=y.Nr in Append_c " << x.Nr << y.Nr << '\n';
      exit(1);
   }
   matrix z(x.Nr,x.Nc+y.Nc);
   for(int i=0;i<z.Nr;i++) {
      for(int j=0;j<z.Nc;j++) {
         if(j<x.Nc) z.t[i*z.Nc+j]=x.t[i*x.Nc+j];
	 else z.t[i*z.Nc+j]=y.t[i*y.Nc+j-x.Nc];
      }   
   }
   return z;
}

matrix Append_r(const matrix& x,const matrix& y)
{
   if(x.Nc!=y.Nc) {
      cout << "x.Nc!=y.Nc in Append_r " << x.Nc << y.Nc << '\n';
      exit(1);
   }
   matrix z(x.Nr+y.Nr,x.Nc);
   for(int i=0;i<z.Nr;i++) {
      for(int j=0;j<z.Nc;j++) {
         if(i<x.Nr) z.t[i*z.Nc+j]=x.t[i*x.Nc+j];
	 else z.t[i*z.Nc+j]=y.t[(i-x.Nr)*y.Nc+j];
      }   
   }
   return z;
}

matrix SubMatrix(const matrix& x,int i1,int i2,int j1,int j2)
{
  int nr=i2-i1+1;
  int nc=j2-j1+1;
   matrix z(nr,nc);
   for(int i=i1;i<=i2;i++) {
      for(int j=j1;j<=j2;j++) {
         z.t[(i-i1)*z.Nc+(j-j1)]=x.t[i*x.Nc+j];
      }   
   }
   return z;
}

void Normalize(matrix& A,double* vec)
{
  int i,j,k;
  double rnorm,remax,absa;
  matrix B(A.Nc,A.Nc);
  int* ikeep=new int[A.Nc];
  int* naxis=new int[A.Nc];
  double* bvec=new double[2*A.Nc];
  for(i=0;i<A.Nc;i++) ikeep[i]=0;
  for(i=0;i<A.Nc;i++) {
    rnorm=0.;
    remax=0.;
    for(j=0;j<A.Nc;j++) {
      rnorm=rnorm+A.t[j*A.Nc+i]*A.t[j*A.Nc+i];
    }
    rnorm=1./sqrt(rnorm);
    for(j=0;j<A.Nc;j++) {
      B.t[j*A.Nc+i]=A.t[j*A.Nc+i]*rnorm;
      absa=fabs(B.t[j*A.Nc+i]);
      if(absa>remax && ikeep[j]==0) {
	remax=absa;
	naxis[i]=j;
      }
    }
    ikeep[naxis[i]]=1;
  }

  for(i=0;i<A.Nc;i++) {
    if(ikeep[i]==0) {
      cout << "matrix.f : Normal axis is not determined";
      exit(1);
    }
  }
  for(i=0;i<A.Nc;i++) {
    k=naxis[i];
    for(j=0;j<A.Nc;j++) {
      A.t[j*A.Nc+k]=B.t[j*A.Nc+i];
    }
    bvec[2*k]=vec[2*i];
    bvec[2*k+1]=vec[2*i+1];
  }
  for(i=0;i<A.Nc*2;i++) vec[i]=bvec[i];

  delete [] bvec;
  delete [] naxis;
  delete [] ikeep;
}


   

