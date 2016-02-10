#ifndef LIN_MAP_H
#define LIN_MAP_H
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//               lin_map.h
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <iostream>
using std::ostream;
#include <cmath>
using std::atan2;

#include <matrix.h>

extern void teigen(double*,double*,int ,int );
extern matrix set_NtoP(double*);
extern matrix set_U(double,double,double);
extern matrix set_B(double,double,double,double,double,double);
extern matrix set_R(double,double,double,double);
extern matrix set_H(double,double,double,double,
		    double,double,double,double);
extern matrix set_Bi(double,double,double,double,double,double);
extern matrix set_Ri(double,double,double,double);
extern matrix set_Hi(double,double,double,double,
		     double,double,double,double);

class lin_map : public matrix
{
public:
   int Ntwiss;
   double* eig;
   double* twiss;
   lin_map(void) : matrix(N_cv,N_cv) 
       {eig=new double [2*N_cv];    Ntwiss=N_cv*(N_cv+1)/2-N_cv/2;
         twiss=new double [Ntwiss]; }
   lin_map(const lin_map&);
   ~lin_map(void) { delete [] twiss; delete [] eig; }
   
   lin_map& operator=(const lin_map&);
   lin_map& operator=(const matrix&);
   void EigenSystem(void);
   
   friend matrix operator*(const lin_map& x,const lin_map& y)
      {return ((matrix&) x*(matrix&) y);}
   friend matrix operator*(const lin_map& x,const matrix& y)
      {return ((matrix&) x*y);}
   friend matrix operator*(const matrix& x,const lin_map& y)
      {return ( x*(matrix&) y);}
   friend matrix diag(lin_map&);
   friend lin_map Transpose(const lin_map&);
   friend lin_map Inverse(const lin_map&);
   friend lin_map Sym_trans(const lin_map&,const matrix&);
   friend lin_map SInverse(const lin_map&);
   friend void Normal_axis_sort(lin_map&,matrix&);
   friend void is_symplectic(const lin_map&);
   friend ostream& operator<<(ostream&,lin_map&);
   friend void PrintList(const lin_map&);
   friend void PrintTwiss(const lin_map&);
   
   friend void get_tune(const lin_map& x,double* rmu)
   {  for(int i=0;i<N_cv/2;i++) 
        rmu[i]=atan2(x.eig[4*i+1],x.eig[4*i]);
   }

   friend matrix set_NtoP(double*);
   friend matrix set_U(double,double,double);
   friend matrix set_B(double,double,double,double,double,double);
   friend matrix set_R(double,double,double,double);
   friend matrix set_H(double,double,double,double,
	      double,double,double,double);
   friend matrix set_Bi(double,double,double,double,double,double);
   friend matrix set_Ri(double,double,double,double);
   friend matrix set_Hi(double,double,double,double,
	      double,double,double,double);
   friend matrix set_H(const lin_map T) {
     return set_H(T.twiss[10],T.twiss[11],T.twiss[12],T.twiss[13],
		  T.twiss[14],T.twiss[15],T.twiss[16],T.twiss[17]);
   }
   friend matrix set_Hi(const lin_map T) {
     return set_Hi(T.twiss[10],T.twiss[11],T.twiss[12],T.twiss[13],
		  T.twiss[14],T.twiss[15],T.twiss[16],T.twiss[17]);
   }
   friend matrix set_R(const lin_map T) {
     return set_R(T.twiss[6],T.twiss[7],T.twiss[8],T.twiss[9]);
   }
   friend matrix set_Ri(const lin_map T) {
     return set_Ri(T.twiss[6],T.twiss[7],T.twiss[8],T.twiss[9]);
   }


};


#endif




