#include <iostream>
using std::ostream;
using std::cout;
using std::ios;
#include <cstdlib>
using std::exit;
#include <cmath>
using std::fabs;
using std::sqrt;
using std::sin;
using std::cos;
using std::tan;

#include <dacpp.h>
#include <lin_map.h>
#include <phys_const.h>
//
//--------------------------  linmap ------------------------------
//

lin_map::lin_map(const lin_map& x) : matrix(x)
{
  int i;
  Ntwiss=N_cv*(N_cv+1)/2-N_cv/2;
  eig=new double [2*N_cv];
  twiss=new double [Ntwiss];
  for(i=0;i<2*N_cv;i++) eig[i]=x.eig[i];
  for(i=0;i<Ntwiss;i++) twiss[i]=x.twiss[i];
}

lin_map& lin_map::operator=(const lin_map& x) 
{
  int i;
  matrix::operator=(x);
  Ntwiss=x.Ntwiss;
  for(i=0;i<2*N_cv;i++) eig[i]=x.eig[i];
  for(i=0;i<Ntwiss;i++) twiss[i]=x.twiss[i];
  return *this;
}
lin_map& lin_map::operator=(const matrix& x) 
{
  int i;
  if(getNr(x)!=N_cv || getNc(x)!=N_cv) {
    cout << " x.Nr!=N_cv or x.Nc!=N_cv " << getNr(x) << getNc(x) << N_cv;
    exit(1);
  }
  Ntwiss=N_cv*(N_cv+1)/2-N_cv/2;
  matrix::operator=(x);
  for(i=0;i<2*N_cv;i++) eig[i]=0.;
  for(i=0;i<Ntwiss;i++) twiss[i]=0.;
  return *this;
}


lin_map Transpose(const lin_map& x)
{
  int i;
  lin_map y;
  y=Transpose((matrix) x);
  for(i=0;i<N_cv;i++) { 
    y.eig[2*i]=x.eig[2*i];
    y.eig[2*i+1]=x.eig[2*i+1]; }
  for(i=0;i<x.Ntwiss;i++) y.twiss[i]=0.;
  y.Ntwiss=x.Ntwiss;
  return y;
}

void is_symplectic(const lin_map& x)
{
  is_symplectic((matrix&) x);
}

ostream& operator<<(ostream& s, lin_map& x)
{
  int i;
  s << "\nLinear map print out  \n\n";

  s << (matrix&) x;
  s << "\n exp(i mu) \n";
  for(i=0;i<N_cv;i++) {
    s.width(10);
    s << x.eig[2*i];}
  s << '\n';
  for(i=0;i<N_cv;i++) {
    s.width(10);
    s << x.eig[2*i+1];}
  s << '\n';
  return s;
}

matrix diag(lin_map& x)
{
  matrix R(N_cv,N_cv);

  int i,j;

  R=diag((matrix&) x,x.eig);
  // Normarization
  for(i=0;i<N_cv;i+=2) {
    double det=0.,f_nor;
    for(j=0;j<N_cv;j+=2) 
      det+=R[j][i]*R[j+1][i+1]
	-R[j+1][i]*R[j][i+1];
    f_nor=1./sqrt(fabs(det));
    for(j=0;j<N_cv;j+=2) {
      R[j][i]*=f_nor;
      R[j+1][i]*=f_nor;
      if(det>0) {
	R[j][i+1]*=f_nor;
	R[j+1][i+1]*=f_nor;
      }
      else {
	R[j][i+1]*=-f_nor;
	R[j+1][i+1]*=-f_nor;
	x.eig[2*i+1]=-x.eig[2*i+1];
	x.eig[2*(i+1)+1]=-x.eig[2*(i+1)+1];
      }
    }
  }
  return R;
}

void Normal_axis_sort(lin_map& x,matrix& R)
{
  lin_map y=x;
  matrix Rn(N_cv,N_cv);
  int i,j,k,Naxis[10],reserved;
  double det_max,det;
  for(i=0;i<N_cv;i+=2) {
    det_max=0.;
    Naxis[i/2]=0;
    for(j=0;j<N_cv;j+=2) {
      reserved=0;
      for(k=0;k<j/2;k++) {
	if(Naxis[k]==j/2)  reserved=1; }
      if(!reserved) {
	det=R[j][i]*R[j+1][i+1]
	  -R[j+1][i]*R[j][i+1];
	det=fabs(det);
	if(det_max<det) {
	  det_max=det;
	  Naxis[i/2]=j/2;
	}
      }
    }
  }

  for(i=0;i<N_cv;i+=2) {
    k=Naxis[i/2]*2;
    for(j=0;j<N_cv;j+=2) {
      Rn[j][k]=R[j][i];
      Rn[j+1][k]=R[j+1][i];
      Rn[j][k+1]=R[j][i+1];
      Rn[j+1][k+1]=R[j+1][i+1];
    }
    x.eig[2*k]=y.eig[2*i];
    x.eig[2*k+1]=y.eig[2*i+1];
    x.eig[2*(k+1)]=y.eig[2*(i+1)];
    x.eig[2*(k+1)+1]=y.eig[2*(i+1)+1];
  
  }
  //R=Rn;
  //Phase space rotation as R(1,2)=R(3,4)=R(5,6)=0

  double rot,sin_r,cos_r;
  for(i=0;i<N_cv;i+=2) {
    rot=atan2(Rn[i][i+1],Rn[i][i]);
    sin_r=sin(rot);
    cos_r=cos(rot);
    for(j=0;j<N_cv;j+=2) {
      R[j][i]=cos_r*Rn[j][i]+sin_r*Rn[j][i+1];
      R[j][i+1]=-sin_r*Rn[j][i]+cos_r*Rn[j][i+1];
      R[j+1][i]=cos_r*Rn[j+1][i]+sin_r*Rn[j+1][i+1];
      R[j+1][i+1]=-sin_r*Rn[j+1][i]+cos_r*Rn[j+1][i+1];
      }
   }
   //cout << "\n Rn and R\n\n" << Rn << "\n\n" << R;
}

lin_map Inverse(const lin_map& x)
{
   lin_map y;
   y=Inverse((matrix&) x);
   return y;
}

lin_map SInverse(const lin_map& x)
{
   lin_map y;
   matrix S(N_cv,N_cv);
   S.Symp();
   y=-S*Transpose((matrix&) x)*S;
   return y;
}

lin_map Sym_trans(const lin_map& x,const matrix& y)
{
   lin_map z;
   matrix S(N_cv,N_cv);
   S.Symp();
   z=-S*Transpose(y)*S*x*y;
   for(int i=0;i<N_cv;i++) { z.eig[2*i]=x.eig[2*i];
        z.eig[2*i+1]=x.eig[2*i+1]; }
   return z;
}

matrix set_U(double nux,double nuy,double nuz)
{
  matrix U(6,6);
  double mux,muy,muz;

  U=0.;
  mux=nux*pi2; muy=nuy*pi2; muz=nuz*pi2;

  U[0][0]=cos(mux);
  U[1][1]=cos(mux);
  U[2][2]=cos(muy);
  U[3][3]=cos(muy);
  U[4][4]=cos(muz);
  U[5][5]=cos(muz);

  U[0][1]=sin(mux);
  U[1][0]=-sin(mux);
  U[2][3]=sin(muy);
  U[3][2]=-sin(muy);
  U[4][5]=sin(muz);
  U[5][4]=-sin(muz);

  return U;
}
matrix set_B(double ax,double bx,double ay,double by,
		   double az,double bz)
{
  matrix B(N_cv,N_cv);
  B=0.;
  double sqrb=sqrt(bx);
  B[1][1]=sqrb;
  B[1][0]=ax/sqrb;
  B[0][0]=1./sqrb;
  sqrb=sqrt(by);
  B[3][3]=sqrb;
  B[3][2]=ay/sqrb;
  B[2][2]=1./sqrb;
  sqrb=sqrt(bz);
  B[5][5]=sqrb;
  B[5][4]=az/sqrb;
  B[4][4]=1./sqrb;
  return B;
}

matrix set_Bi(double ax,double bx,double ay,double by,
		   double az,double bz)
{
  matrix B(N_cv,N_cv);
  B=0.;

  double sqrb=sqrt(bx);
  B[0][0]=sqrb;
  B[1][0]=-ax/sqrb;
  B[1][1]=1./sqrb;
  sqrb=sqrt(by);
  B[2][2]=sqrb;
  B[3][2]=-ay/sqrb;
  B[3][3]=1./sqrb;
  sqrb=sqrt(bz);
  B[4][4]=sqrb;
  B[5][4]=-az/sqrb;
  B[5][5]=1./sqrb;
  return B;
}

matrix set_R(double r1,double r2,double r3,double r4)
{
// This is R^{-1} in Beam envelope paper (PRE49)
  matrix R(N_cv,N_cv);
  double rmu=sqrt(1.-r1*r4+r2*r3);

  R=0.;

  R[0][0]=rmu;
  R[1][1]=rmu;
  R[2][2]=rmu;
  R[3][3]=rmu;
  R[4][4]=1.;
  R[5][5]=1.;

  R[2][0]=r1;
  R[2][1]=r2;
  R[3][0]=r3;
  R[3][1]=r4;
  R[0][2]=-r4;
  R[0][3]=r2;
  R[1][2]=r3;
  R[1][3]=-r1;

  return R;
}

matrix set_Ri(double r1,double r2,double r3,double r4)
{
// This is R^{-1} in Beam envelope paper (PRE49)
  matrix R(N_cv,N_cv);
  double rmu=sqrt(1.-r1*r4+r2*r3);

  R=0.;

  R[0][0]=rmu;
  R[1][1]=rmu;
  R[2][2]=rmu;
  R[3][3]=rmu;
  R[4][4]=1.;
  R[5][5]=1.;

  R[2][0]=-r1;
  R[2][1]=-r2;
  R[3][0]=-r3;
  R[3][1]=-r4;
  R[0][2]=r4;
  R[0][3]=-r2;
  R[1][2]=-r3;
  R[1][3]=r1;

  return R;
}

matrix set_H(double ex,double epx,double ey,double epy,
		   double zx,double zpx,double zy,double zpy)
// This is H in Beam envelope paper (PRE49)
{
  matrix H(N_cv,N_cv);
  double a,ai;

  H=0.;

  double detHx=zx*epx-zpx*ex;
  double detHy=zy*epy-zpy*ey;
  double a2=1.-detHx-detHy;
  if(a2>0.) a=sqrt(a2);
  else { cout << " det Hx +det Hy >1 " << detHx << detHy;
	 exit(1); }
  ai=1./(1.+a);

  H[0][0]=1.-detHx*ai;
  H[1][1]=H[0][0];
  H[2][2]=1.-detHy*ai;
  H[3][3]=H[2][2];
  H[4][4]=a;
  H[5][5]=H[4][4];

  H[0][4]=-zx;
  H[0][5]=-ex;
  H[1][4]=-zpx;
  H[1][5]=-epx;
  H[2][4]=-zy;
  H[2][5]=-ey;
  H[3][4]=-zpy;
  H[3][5]=-epy;

  H[4][0]=epx;
  H[4][1]=-ex;
  H[5][0]=-zpx;
  H[5][1]=zx;
  H[4][2]=epy;
  H[4][3]=-ey;
  H[5][2]=-zpy;
  H[5][3]=zy;

  H[0][2]=(ex*zpy-epy*zx)*ai;
  H[0][3]=(ey*zx-ex*zy)*ai;
  H[1][2]=(epx*zpy-epy*zpx)*ai;
  H[1][3]=(ey*zpx-epx*zy)*ai;
  H[2][0]=(ey*zpx-epx*zy)*ai;
  H[2][1]=(-ey*zx+ex*zy)*ai;
  H[3][0]=(epy*zpx-epx*zpy)*ai;
  H[3][1]=(ex*zpy-epy*zx)*ai;

  return H;
}

matrix set_Hi(double ex,double epx,double ey,double epy,
		   double zx,double zpx,double zy,double zpy)
// This is H^{-1} in Beam envelope paper (PRE49)
{
  matrix H(N_cv,N_cv);
  double a,ai;

  H=0.;

  double detHx=zx*epx-zpx*ex;
  double detHy=zy*epy-zpy*ey;
  double a2=1.-detHx-detHy;
  if(a2>0.) a=sqrt(a2);
  else { cout << " det Hx +det Hy >1 " << detHx << detHy;
	 exit(1); }
  ai=1./(1.+a);

  H[0][0]=1.-detHx*ai;
  H[1][1]=H[0][0];
  H[2][2]=1.-detHy*ai;
  H[3][3]=H[2][2];
  H[4][4]=a;
  H[5][5]=H[4][4];

  H[0][4]=zx;
  H[0][5]=ex;
  H[1][4]=zpx;
  H[1][5]=epx;
  H[2][4]=zy;
  H[2][5]=ey;
  H[3][4]=zpy;
  H[3][5]=epy;

  H[4][0]=-epx;
  H[4][1]=ex;
  H[5][0]=zpx;
  H[5][1]=-zx;
  H[4][2]=-epy;
  H[4][3]=ey;
  H[5][2]=zpy;
  H[5][3]=-zy;

  H[0][2]=(ex*zpy-epy*zx)*ai;
  H[0][3]=(ey*zx-ex*zy)*ai;
  H[1][2]=(epx*zpy-epy*zpx)*ai;
  H[1][3]=(ey*zpx-epx*zy)*ai;
  H[2][0]=(ey*zpx-epx*zy)*ai;
  H[2][1]=(-ey*zx+ex*zy)*ai;
  H[3][0]=(epy*zpx-epx*zpy)*ai;
  H[3][1]=(ex*zpy-epy*zx)*ai;

  return H;
}

matrix set_NtoP(double* twiss)
{
   matrix Bi(N_cv,N_cv),Ri(N_cv,N_cv),Hi(N_cv,N_cv);
   Bi=set_Bi(twiss[0],twiss[1],twiss[2],twiss[3],twiss[4],twiss[5]);
   Ri=set_Ri(twiss[6],twiss[7],twiss[8],twiss[9]);
   Hi=set_Hi(twiss[10],twiss[11],twiss[12],twiss[13],
	   twiss[14],twiss[15],twiss[16],twiss[17]);
   //cout << Bi << Ri << Hi;
   //cout.flush();
   return (Hi*Ri*Bi);
}

matrix set_Cros(double x_angle)
{
  matrix T(6,6);
  double cosxi,tanx;
  cosxi=1./cos(x_angle);
  tanx=tan(x_angle);

  T=0.;
  T[0][0]=1.;
  T[1][1]=cosxi;
  T[2][2]=1.;
  T[3][3]=cosxi;
  T[4][4]=cosxi;
  T[5][5]=1.;
  T[0][4]=tanx;
  T[5][1]=-tanx;
  return T;
}
void lin_map::EigenSystem()
{
   matrix V0(N_cv,N_cv);
   V0=diag(*this);
   Normal_axis_sort(*this,V0);
   matrix U0(N_cv,N_cv),UV(N_cv,N_cv),VH(N_cv,N_cv);
   matrix S2(2,2),Hx(2,2),Hy(2,2);
   matrix H(6,6),B(6,6),R(6,6);

   double aa;
   {
   matrix UV11(2,2),UV12(2,2),UV13(2,2);
   matrix UV21(2,2),UV22(2,2),UV23(2,2);
   matrix UV31(2,2),UV32(2,2),UV33(2,2);
   matrix H11(2,2),H12(2,2),H13(2,2);
   matrix H21(2,2),H22(2,2),H23(2,2);
   matrix H31(2,2),H32(2,2),H33(2,2);
   matrix H1(2,6),H2(2,6),H3(2,6);
   U0=SInverse(V0)*(*this)*V0;
//   cout << "V0 U0\n" << V0 << U0 << SInverse(V0)*V0;
   V0=SInverse(V0);
   UV=SInverse(U0)*V0;
   UV11=SubMatrix(UV,0,1,0,1);
   UV12=SubMatrix(UV,0,1,2,3);
   UV13=SubMatrix(UV,0,1,4,5);
   UV21=SubMatrix(UV,2,3,0,1);
   UV22=SubMatrix(UV,2,3,2,3);
   UV23=SubMatrix(UV,2,3,4,5);
   UV31=SubMatrix(UV,4,5,0,1);
   UV32=SubMatrix(UV,4,5,2,3);
   UV33=SubMatrix(UV,4,5,4,5);
   S2.Symp();
   if((aa=Det(UV33))>0) {
      aa=sqrt(aa);
   }
   else { cout << "Not supported now.(H)\n" << aa; exit(1);}
   Hx=-S2*Transpose(UV31)*S2*UV33/aa;
   Hy=-S2*Transpose(UV32)*S2*UV33/aa;
//   cout << " Hx and Hy\n" << Hx << Hy;cout.flush();
   if((aa=1.-Det(Hx)-Det(Hy))>0) aa=sqrt(aa);
   else { cout << "Not supported now.(H)\n" << aa; exit(1);} 
   H11=(1.-Det(Hx)/(1.+aa))*IdentityMatrix(2);
   H22=(1.-Det(Hy)/(1.+aa))*IdentityMatrix(2);
   H12=Hx*S2*Transpose(Hy)*S2/(1.+aa);
   H21=Hy*S2*Transpose(Hx)*S2/(1.+aa);
   H13=-Hx;
   H23=-Hy;
   H31=-S2*Transpose(Hx)*S2;
   H32=-S2*Transpose(Hy)*S2;
   H33=aa*IdentityMatrix(2);
   
   H1=Append_c(Append_c(H11,H12),H13);
   H2=Append_c(Append_c(H21,H22),H23);
   H3=Append_c(Append_c(H31,H32),H33);
   H=Append_r(Append_r(H1,H2),H3);
//   cout << " H\n" << H;
   VH=V0*SInverse(H);
//   cout << "VH\n" << VH;
   }
   {
   matrix U11(2,2),U12(2,2);
   matrix U21(2,2),U22(2,2);
   matrix R11(2,2),R12(2,2),R13(2,2);
   matrix R21(2,2),R22(2,2),R23(2,2);
   matrix R31(2,2),R32(2,2),R33(2,2);
   matrix R1(2,6),R2(2,6),R3(2,6);
   U11=SubMatrix(VH,0,1,0,1);
   U12=SubMatrix(VH,0,1,2,3);
   U21=SubMatrix(VH,2,3,0,1);
   U22=SubMatrix(VH,2,3,2,3);
   if((aa=Det(U22))>0) aa=sqrt(aa);
   else { cout << "Not supported now.(R)\n"; exit(1);} 
   R21=-S2*Transpose(U22)*S2*U21/aa;
   if((aa=1.-Det(R21))>0) aa=sqrt(aa);
   else { cout << "Not supported now.(R)\n"; exit(1);}
   R11=aa*IdentityMatrix(2);
   R12=S2*Transpose(R21)*S2;
   R22=R11;
   R13.Clear();
   R23.Clear();
   R31.Clear();
   R32.Clear();
   R33=IdentityMatrix(2);
   R1=Append_c(Append_c(R11,R12),R13);
   R2=Append_c(Append_c(R21,R22),R23);
   R3=Append_c(Append_c(R31,R32),R33);
   R=Append_r(Append_r(R1,R2),R3);
   }
   {
   B=VH*SInverse(R);
   }
//   cout << "H matrix\n" << H;
//   cout << "R matrix\n" << R;
//   cout << "B matrix\n" << B; cout.flush();
//   cout << B*R*H << V0 << *this;
   twiss[1]=B[1][1]*B[1][1];
   twiss[0]=B[1][0]*B[1][1];
   twiss[3]=B[3][3]*B[3][3];
   twiss[2]=B[3][2]*B[3][3];
   twiss[5]=B[5][5]*B[5][5];
   twiss[4]=B[5][4]*B[5][5];

   twiss[6]=R[2][0];
   twiss[7]=R[2][1];
   twiss[8]=R[3][0];
   twiss[9]=R[3][1];

   twiss[10]=-H[0][5];
   twiss[11]=-H[1][5];
   twiss[12]=-H[2][5];
   twiss[13]=-H[3][5];
   twiss[14]=-H[0][4];
   twiss[15]=-H[1][4];
   twiss[16]=-H[2][4];
   twiss[17]=-H[3][4];

}

void PrintList(const lin_map& x)
{
  int i;
  cout << '{';
  PrintList((matrix) x);
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(15);
  cout << ",{{";
  for(i=0;i<N_cv*2;i++) {
    cout << x.eig[i];
    if(i!=N_cv*2-1) cout << ',';
  }
  cout << "},\n{";
  for(i=0;i<x.Ntwiss;i++) {
    cout << x.twiss[i]; 
    if(i!=x.Ntwiss-1) cout << ',';
  }
  cout << "}}};\n";
  cout.unsetf(ios::scientific);
  cout.unsetf(ios::floatfield);
}

void PrintTwiss(const lin_map& x)
{
  cout << "***      AX   ";
  cout.width(12);  cout << x.twiss[0];
  cout << "   ***      BX   ";
  cout.width(12);  cout << x.twiss[1] << '\n';
  cout << "***      AY   ";
  cout.width(12);  cout << x.twiss[2];
  cout << "   ***      BY   ";
  cout.width(12);  cout << x.twiss[3] << '\n';
  cout << "***      AZ   ";
  cout.width(12);  cout << x.twiss[4];
  cout << "   ***      BZ   ";
  cout.width(12);  cout << x.twiss[5] << '\n';

  cout << "***      R1   ";
  cout.width(12);  cout << x.twiss[6];
  cout << "   ***      R2   ";
  cout.width(12);  cout << x.twiss[7] << '\n';
  cout << "***      R3   ";
  cout.width(12);  cout << x.twiss[8];
  cout << "   ***      R4   ";
  cout.width(12);  cout << x.twiss[9] << '\n';

  cout << "***      EX   ";
  cout.width(12);  cout << x.twiss[10];
  cout << "   ***     EPX   ";
  cout.width(12);  cout << x.twiss[11] << '\n';
  cout << "***      EY   ";
  cout.width(12);  cout << x.twiss[12];
  cout << "   ***     EPY   ";
  cout.width(12);  cout << x.twiss[13] << '\n';
  cout << "***      ZX   ";
  cout.width(12);  cout << x.twiss[14];
  cout << "   ***     ZPX   ";
  cout.width(12);  cout << x.twiss[15] << '\n';
  cout << "***      ZY   ";
  cout.width(12);  cout << x.twiss[16];
  cout << "   ***     ZPY   ";
  cout.width(12);  cout << x.twiss[17] << '\n';

}
