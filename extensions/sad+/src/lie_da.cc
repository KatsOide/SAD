#include <iostream>
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>

#include <map_da.h>
#include <lin_map.h>
#include <lie_da.h>

static char label[7]="xyzpqe";

void is_symplectic(const map_da& x)
{
   da z;
   cout << "-------------------------------------------\n";
   cout << "Symplecticity check\n";
   for(int i=0;i<N_cv;i++) {
      for(int j=i+1;j<N_cv;j++) {
        z=poi(x.m[i],x.m[j]);
	cout << "\n[" << label[i] << ", " 
             << label[j] << "]=" ;
	z.daprint(N_ord-1);
      }
   }
}

//void Normalize(const map_da& x0,lin_map& Z,map_da& x)
map_da Normalize(const map_da& x0)
{
  lin_map T;
  T=lin_da(x0);
  matrix R(N_cv,N_cv);
  cout << " \n Linear map is \n" << T;
  R=diag(T);
  Normal_axis_sort(T,R);
  cout << "Normalized eigen vector matrix\n" << R;
  //lin_map Z;
  //Z=Sym_trans(T,R);	// R^{-1}*T*R
  //cout << "R^{-1}TR\n" << Z;
  //cout.flush();

  return Sym_trans(x0,R);		// R^{-1}*x2*R
}
 
map_da Normalize(const map_da& x0,lin_map& Z)
{
  lin_map T;
  T=lin_da(x0);
  matrix R(N_cv,N_cv);
  cout << " \n Linear map is \n" << T;
  R=diag(T);
  Normal_axis_sort(T,R);
  cout << "\n\nNormalized eigen vector matrix\n" << R;
  Z=Sym_trans(T,R);	// R^{-1}*T*R
  //cout << "R^{-1}TR\n" << Z;

  return Sym_trans(x0,R);		// R^{-1}*x2*R
}
	 
da fac_drg_type1(const map_da& x0)
{
  cout << " DF type 1 factrization as x0=R e^{:F_3:}....e^{:F_n:}\n";
  lin_map Z,Zi;
  map_da x;
  x=Normalize(x0);	//T which is got like this does not have eigen values.	
  Z=lin_da(x);
  Zi=SInverse(Z);
// x=Z^{-1} x
  x=x*Zi;
  //cout << " For debug \n" << x;

   da RFn,RFi,f;
   RFn=0.;
   for(int j=3;j<=N_ord;j++) {
      f=poi_itg(x);
      RFi=msk(f,j);

      RFn+=RFi;
      x=lie_exp(-RFi,x);
   }
   //cout << " For debug 2\n" << fac_map_type1(RFn);
   return RFn;
}

da fac_drg_type3(const map_da& x0)
{
  cout << " DF type 3 factrization as x0=R e^{:F:}\n";
  lin_map Z,Zi;
  map_da x,xn;
  xn=Normalize(x0,Z);      //T which is got like this does not have eigen values.
  //  Z=lin_da(xn);
  cout << "\n\n Normalized linear map \n" << Z;
  Zi=SInverse(Z);
// x=Z^{-1} x
  xn=xn*Zi;
  //cout << " For debug \n" << xn;
  x=xn;

   da RFn,RFi,f;
   RFn=0.;
   for(int j=3;j<=N_ord;j++) {
      f=poi_itg(x);
      RFi=msk(f,j);
      RFn+=RFi;
      x=lie_exp(-RFn,xn);
   }
   //cout << " For debug 2\n" << fac_map_type3(RFn);
   return RFn;
}

da fac_drg_type1(const map_da& x0,map_da& x)
{
  cout << " DF type 1 factrization as x0=R e^{:F_3:}....e^{:F_n:}\n";
  lin_map Z,Zi;
  x=Normalize(x0);      //T which is got like this does not have eigen values.
  Z=lin_da(x);
  Zi=SInverse(Z);
// x=Z^{-1} x
  x=x*Zi;

  da RFn,RFi,f;
   RFn=0.;
   for(int j=3;j<=N_ord;j++) {
      f=poi_itg(x);
      RFi=msk(f,j);
      RFn+=RFi;
      x=lie_exp(-RFi,x);  // x=e^{-:F_i:} x
   }
   return RFn;
}

da fac_drg_type2(const map_da& x0)
{
  cout << " DF type 2 factrization as x0=e^{:F_n:}....e^{:F_3:} R\n";
  lin_map Z,Zi;
  map_da x;
  x=Normalize(x0);      //T which is got like this does not have eigen values.
  Z=lin_da(x);
  Zi=SInverse(Z);
// x=Z^{-1} x
  x=Zi*x;
  //cout << " For debug \n" << x;

   da RFn,RFi,f;
   RFn=0.;
   for(int j=3;j<=N_ord;j++) {
      f=poi_itg(x);
      RFi=msk(f,j);

      RFn+=RFi;
      x=lie_exp(-RFi)*x;
   }
   //cout << " For debug 2\n" << fac_map_type2(RFn);
   return RFn;
}

da fac_drg_type3(const map_da& x0,map_da& x)
{
  cout << " DF type 3 factrization as x0=R e^{:F:}\n";
  lin_map Z,Zi;
  x=Normalize(x0,Z);      //T which is got like this does not have eigen values.
  //  Z=lin_da(x);
  cout << "\n\n Normalized linear map \n" << Z;
  Zi=SInverse(Z);
// x=Z^{-1} x
  x=x*Zi;

   da RFn,RFi,f;
   RFn=0.;
   for(int j=3;j<=N_ord;j++) {
      f=poi_itg(x);
      RFi=msk(f,j);
      RFn+=RFi;
      x=lie_exp(-RFn,x0);
   }
   return RFn;
}

map_da fac_map_type1(const da& f)
{
   map_da x;
   da RFn;
   x.dBase();
   for(int j=N_ord;j>=3;j--) {
      RFn=msk(f,j);
      x=lie_exp(RFn,x);	 
   }
   return x;
}

map_da fac_map_type2(const da& f)
{
   map_da x;
   da RFn;
   x.dBase();
   for(int j=3;j<=N_ord;j++) {
      RFn=msk(f,j);
      x=lie_exp(RFn,x);	 
   }
   return x;
}

map_da fac_map_type3(const da& f)
{
   map_da x;
   x.dBase();
   x=lie_exp(f,x);
   return x;
}

da poi_itg(const map_da& x)
{
   
   //  integrate poisson bracket through a straight path.

   da z1,z2,f;
   f=0.;
   
   int Ncv2=N_cv/2;
   for(int i=0;i<Ncv2;i++) {
      z1.dBase(Ncv2+i,1.);
      z2=line_itg(x.m[i]);
// sign should be checked.
      f-=muln(z2,z1,1);
      z1.dBase(i,1.);
      z2=line_itg(x.m[i+Ncv2]);
      f+=muln(z2,z1,1);
   }
   return f;
}



map_da lie_exp(const da& f,const map_da& x)
{

   map_da z;
   double a;
   int N=min_ord(f);
   if(N<=1) {
      cout << " Can not expand lie_exp : " << N << '\n';
      exit(1);
   }
   int Nitr=(N_ord-1)/(N-2);

   for(int j=0;j<N_cv;j++) {
      z.m[j]=x.m[j];
      da y=x.m[j];
      for(int i=0;i<Nitr;i++) {
         if(i==0) {
	    y=poi(f,y);
	 } 
	 else {
            a=1./((double)(i+1));
	    y=a*poi(f,y);
	 }
	 z.m[j]+=y;
      }
   }
   return z;
}

map_da lie_exp(const da& f)
{

   map_da x,z;
   double a;
   x.dBase();
   int N=min_ord(f);
   if(N<=1) {
      cout << " Can not expand lie_exp : " << N << '\n';
      exit(1);
   }
   int Nitr=(N_ord-1)/(N-2);

   for(int j=0;j<N_cv;j++) {
      z.m[j]=x.m[j];
      da y=x.m[j];
      for(int i=0;i<Nitr;i++) {
         if(i==0) {
	    y=poi(f,y);
	 } 
	 else {
            a=1./((double)(i+1));
	    y=a*poi(f,y);
	 }
	 z.m[j]+=y;
      }
   }
   return z;
}

da lie_exp(const da& f,const da& x)
{
   da z,z1;
   double a;
   int N=min_ord(f);
   if(N<=1) {
      cout << " Can not expand lie_exp : " << N << '\n';
      exit(1);
   }
   int Nitr=(N_ord-1)/(N-2);

   z=x;
   z1=x;
   for(int i=0;i<Nitr;i++) {
      a=1./((double)(i+1));
      if(i==0) {z1=poi(f,z1);} 
      else {
         a=1./((double)(i+1));
	 z1=a*poi(f,z1);
      }
      z+=z1;
   }
   return z;
}

#include <Complex.h>
#include <c_da.h>
#include <map_c_da.h>
da can_perturbation(const map_da& x)
{
  lin_map T,Ti,Z;
  map_da xm;

  xm=Normalize(x,T);

//  T=lin_da(xm);  T which is got like this does not have eigen values.
  Ti=SInverse(T);

  cout << "\n\n Normalized linear map \n" << T;
  //<< "\n Normalized map\n" << xm;
   map_da xf;
   c_da y,y2;
   map_c_da z,zi;
   da RFi,f,RFn;
   RFn=0.;
  double rmu[3];
  get_tune(T,rmu);

// z=((x+px)/2, (x-px)/2i )   zi=((x + ipx), (x-ipx) )
   z.NBase();
   zi.NiBase();

// xf=Ti xm 
  xf=xm*Ti;
//      cout << "\n xm" << xm;
//   cout << "\n Map : Non-Linear effect" << xf;
   f=poi_itg(xf);
//   cout << "\n 3rd order integral\n" << f;

   for(int j=3;j<=N_ord;j++) {
      RFi=msk(f,j);
//    cout << "\n RFi" << RFi;
      y=concatenate(RFi,z,1);
      y2=p_gen_f(rmu,y);
      y=concatenate(y2,zi,1);
//    xm=e^{:y:} xm e^{-:y:}
      xm=lie_exp(-real(y))*lie_exp(real(y),xm);
//      xf=Ti xm 
	xf=xm*Ti;
	for(int i=3;i<=j;i++) {
	   f=poi_itg(xf);
	   RFi=msk(f,i);
//	   cout << "H " << RFi;
	   xf=lie_exp(-RFi,xf);
	}
	RFn+=RFi;
	f=poi_itg(xf); 
//	cout << '\n' << j+1 << "th order integral\n" << f;
   }
// return diagonalized Hamiltonian
   return RFn;
}


da can_perturbation(const map_da& x,da& Fc)
{
  lin_map T,Ti,Z;
  map_da xm;

  Fc=0.;
  xm=Normalize(x,T);

//  T=lin_da(xm);  T which is got like this does not have eigen values.
  Ti=SInverse(T);

  cout << "\n\n Normalized linear map \n" << T << "\n\n";
  //<< "\n Normalized map\n" << xm;
   map_da xf;
   c_da y,y2;
   map_c_da z,zi;
   da RFi,f,RFn;
   RFn=0.;
  double rmu[3];
  get_tune(T,rmu);

// z=((x+px)/2, (x-px)/2i )   zi=((x + ipx), (x-ipx) )
   z.NBase();
   zi.NiBase();

// xf=Ti xm 
  xf=xm*Ti;
//      cout << "\n xm" << xm;
//   cout << "\n Map : Non-Linear effect" << xf;
   f=poi_itg(xf);
//   cout << "\n 3rd order integral\n" << f;

   for(int j=3;j<=N_ord;j++) {
      RFi=msk(f,j);
//    cout << "\n RFi" << RFi;
      y=concatenate(RFi,z,1);
      y2=p_gen_f(rmu,y);
      y=concatenate(y2,zi,1);
      Fc+=real(y);
//    xm=e^{:y:} xm e^{-:y:}
      xm=lie_exp(-real(y))*lie_exp(real(y),xm);
//      xf=Ti xm 
	xf=xm*Ti;
	for(int i=3;i<=j;i++) {
	   f=poi_itg(xf);
	   RFi=msk(f,i);
//	   cout << "H " << RFi;
	   xf=lie_exp(-RFi,xf);
	}
	RFn+=RFi;
	f=poi_itg(xf); 
//	cout << '\n' << j+1 << "th order integral\n" << f;
   }
// return diagonalized Hamiltonian
   return RFn;
}

c_da Normal_expression(const da& f)
{	   
   map_c_da z;
   c_da fz;
   z.NBase();
   fz=concatenate(f,z,1);
   return fz;
}

c_da Real_expression(const c_da& f)
{	   
   map_c_da z;
   c_da fz;
   z.NiBase();
   fz=concatenate(f,z,1);
   return fz;
}
