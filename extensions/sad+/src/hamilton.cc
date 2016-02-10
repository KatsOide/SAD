#include <iostream>
using std::cout;
using std::cerr;
using std::ios;
#include <fstream>
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <cmath>
using std::exp;

#include <map_double.h>
#include <map_da.h>
#include <lin_map.h>
#include <lie_da.h>

double max_x[]={0.001,0.001,0.001,0.001,0.001,0.001};
int nsl=5;
double ex=1.e-8;

void hamilton(const map_da& x)
{
   da f,f0;
   c_da fz;
   f=fac_drg_type3(x);
   fz=Normal_expression(f);
   f0=real(diag_part(fz));
   cout << "f_00" << f0;
   
}

da invariant_part(const da& f)
{
   da f0;
   c_da fz;
   fz=Normal_expression(f);
   fz=diag_part(fz);
   f0=real(Real_expression(fz));
   return f0;
}

void nonlinear_dist(const map_da& x)
{
   int i,j,k,l;
   c_da fz;
   da f,fc,Jx,Jy,Jz;
   map_double y;
   map_da x2;
   double z,Jxval;

   f=can_perturbation(x,fc);
   fz=Normal_expression(f);
   // comment out to avoid error  Nov.7.2003
   //   cout << "Diagonalized Hamiltonian \n" << real(diag_part(fz));
   cout << "\n\nCanonical transformation\n" << fc;
   //cout << "\n\nCanonical transformation\n" << Normal_expression(fc);
   x2=lie_exp(fc);
   Jx=(x2[0]*x2[0]+x2[3]*x2[3])*0.5;
   Jy=(x2[1]*x2[1]+x2[4]*x2[4])*0.5;
   Jz=(x2[2]*x2[2]+x2[5]*x2[5])*0.5;
   cout << " Jx \n" << Jx;
   cout << " Jy \n" << Jy;
   cout << " Jz \n" << Jz;

   y=0.;
   ofstream fout("nonlinear.dist",ios::out);
   if(!fout) {
      cerr << "cannot open output file : nonlinear.dist\n";
      exit(1); 
   }
   for(i=0;i<nsl;i++) {
      for(j=0;j<nsl;j++) {
	 for(k=0;k<nsl;k++) {
	    for(l=0;l<nsl;l++) {
	       y[0]=(double)i*max_x[0]/nsl;
	       y[1]=(double)j*max_x[1]/nsl;
	       y[3]=(double)k*max_x[3]/nsl;
	       y[4]=(double)l*max_x[4]/nsl;
	       Jxval=concatenate(Jx,y);
	       cout << Jxval << '\n';
	       z=exp(-Jxval/ex);
	       fout.width(10);
	       fout << y[0];
	       fout.width(10);
	       fout << y[3];
	       fout.width(10);
	       fout << y[1];
	       fout.width(10);
	       fout << y[4]<< " ";
	       fout.width(12);
	       fout << z << '\n';
	    }
	 }
      }
   }
   fout.close();
}

