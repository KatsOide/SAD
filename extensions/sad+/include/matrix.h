#ifndef MATRIX_H
#define MATRIX_H
// template <class T>      // T is double or complex.
#include <iostream>
using std::ostream;

class matrix;

extern matrix IdentityMatrix(int);

class matrix
{
   double* t;
   int Nr,Nc;
 public:
   matrix(int N,int M) {Nr=N; Nc=M; t=new double [Nr*Nc];}
   matrix(const matrix&);
   ~matrix(void) {delete [] t;}
   
   matrix& operator=(const matrix&);
   matrix& operator=(double);
   matrix operator-(void);
   double* operator[](int i) { return (t+i*Nc);}
//   double* operator[](int i) { return (t+(i-1)*Nc-1); }

   void ci(int i,double f) {t[i]=f;}
   friend int getNr(const matrix& x) { return x.Nr;}
   friend int getNc(const matrix& x) { return x.Nc;}
   friend double m_element(const matrix& x,int i,int j) {
      return x.t[i*x.Nc+j];}
   friend matrix operator+(const matrix&,const matrix&);
   friend matrix operator-(const matrix&,const matrix&);
   friend matrix operator*(const matrix&,const matrix&);
   friend matrix operator*(double,const matrix&);
   friend matrix operator*(const matrix&,double);
   friend matrix operator/(const matrix&,double);
   friend matrix Transpose(const matrix&);
   friend matrix Inverse_G(const matrix&);
   friend matrix Inverse(const matrix&);
   friend void teigen(const matrix&,matrix&,double*);
   void Symp(void);
   friend matrix SymplecticMatrix(int);
   friend ostream& operator<<(ostream& s, matrix& A);
   friend matrix diag(const matrix&,double*);
   friend void is_symplectic(const matrix&);
   friend matrix SInverse(const matrix&);
   friend void ReadMathList(char*,matrix&);
   void DiagonalMatrix(double*);
   void I(void);
   void Clear(void);
   friend matrix IdentityMatrix(int);
   friend double Det(const matrix&);
   friend void Normalize(matrix& A,double* vec);

   friend void PrintList(const matrix&);
   friend matrix SubMatrix(const matrix&,int,int,int,int);
   friend matrix Append_c(const matrix&,const matrix&);
   friend matrix Append_r(const matrix&,const matrix&);
   
   //   template <class T>
   //   friend vector<T> operator*(const matrix&,const vector<T>& x);
};

#endif
