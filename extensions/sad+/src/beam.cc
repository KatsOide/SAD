#include <iostream>
using std::ostream;
using std::cout;
#include <cstdio>
using std::sscanf;
#include <cstdlib>
using std::atof;
using std::exit;
#include <cmath>
using std::fabs;
using std::atan2;
using std::sin;
using std::cos;
using std::sqrt;

#include <cstring>
using std::strchr;
using std::strlen;
using std::strstr;

#include <element.h>

#define Ncanv 6

double tran(void);
double ran1(void);
double rgauss(void);
double tgauss(void);

void Normalize(matrix&,double*);

Beam::Beam(const Beam& BEAM) 
{
   energy=BEAM.energy;
   emx=BEAM.emx; emy=BEAM.emy; emz=BEAM.emz; 
   Nparticle=BEAM.Nparticle;
  Beam_envelope=new matrix(6,6);
}

ostream& operator<<(ostream& s,Beam& BEAM)
{
   if(BEAM.energy==0.) return s;
   s << "######### BEAM resource print out #############################\n";
   s << "Energy " << BEAM.energy << " GeV\n";
   s << " sigma_x(y) ";
   s.width(10); s << BEAM.emx;
   s.width(10); s << BEAM.emy;
   s.width(10); s << BEAM.emz;
   s << " \n Number of particles " << BEAM.Nparticle <<'\n';
   return s;
}


void Beam::print(void) 
{
   if(energy==0.) return;
   cout << "######### BEAM resource print out #############################\n";
   cout << "Energy " << energy << " GeV\n";
   cout << " sigma_x(y) ";
   cout.width(10); cout << emx;
   cout.width(10); cout << emy;
   cout.width(10); cout << emz;
   cout << " \n Number of particles " << Nparticle <<'\n'; 
}

void Beam::set_beam_envelope(matrix* NtoP)
{
   matrix Emit(6,6);
   double emits[6];
   emits[0]=emx;emits[1]=emx;
   emits[2]=emy;emits[3]=emy;
   emits[4]=emz;emits[5]=emz;
   Emit.DiagonalMatrix(emits);
   *Beam_envelope=((matrix)(*NtoP))*(matrix) Emit*Transpose((matrix)*NtoP);
}

Beam::Beam(const double* K)
{
   energy=K[0];
   emx=K[1];
   emy=K[2];
   emz=K[3];
   Nparticle=K[4];
  Beam_envelope=new matrix(6,6);
}

void m_zslice(matrix& B,matrix& A,double* v6)
{
  int i,j;
  for(j=0;j<4;j++) {
    for(i=0;i<4;i++) A[i][j]=B[i][j];
    for(i=0;i<4;i++) {
      v6[i]=B[4][i];
      A[4][i]=B[5][i];
      A[i][4]=B[i][5];
    }
    v6[5]=B[4][4];
    v6[4]=B[5][4];
    A[4][4]=B[5][5];
  }
}

double benv_norm_axis(matrix& B,double* v)
{
  int i,j;
  double t_angle,vec[5];
  if(B[0][2]==0.) return 0.;
  t_angle=atan2(2.*B[0][2],B[0][0]-B[2][2])*0.5;
  double cost=cos(t_angle);
  double sint=sin(t_angle);
  matrix U(5,5);
  U=cost;
  U[4][4]=1.;
  U[0][2]=sint;
  U[1][3]=sint;
  U[2][0]=-sint;
  U[3][1]=-sint;

  B=(matrix)U*B*Transpose(U);
  for(i=0;i<5;i++) {
    vec[i]=0.;
    for(j=0;j<5;j++) vec[i]+=U[i][j]*v[j];
  }
  for(i=0;i<5;i++) v[i]=vec[i];

  return t_angle;
}

//extern int N_ord,N_var,L_vec,N_cv;



EMIT::EMIT(const EMIT& A)
{
  Linear_map=new matrix(6,6);
  PtoN=new matrix(6,6);
  Beam_envelope=new matrix(6,6);
  Damping_matrix=new matrix(6,6);
  Diffusion_matrix=new matrix(6,6);
  
  *Linear_map=*A.Linear_map;
  *PtoN=*A.PtoN;
  *Beam_envelope=*A.Beam_envelope;
  *Damping_matrix=*A.Damping_matrix;
  *Diffusion_matrix=*A.Diffusion_matrix;
  for(int i=0;i<6;i++) cod[i]=A.cod[i];
  emx=A.emx; emy=A.emy; emz=A.emz;
}

void EMIT::set_beam_envelope(matrix* NtoP)
{
   matrix Emit(6,6);
   double emits[6];
   emits[0]=emx;emits[1]=emx;
   emits[2]=emy;emits[3]=emy;
   emits[4]=emz;emits[5]=emz;
   Emit.DiagonalMatrix(emits);
   *Beam_envelope=((matrix)*NtoP)*Emit*Transpose(*NtoP);
}

void EMIT::set_EMIT(void)
{
   matrix tmp(6,6);
   double eig[12];
   teigen(*Diffusion_matrix,*DifMeig,eig);
//   tmp=Transpose(*Diffusion_matrix);
//   teigen(tmp.t,eig,6,6);
//   *DifMeig=Transpose(tmp);
   Normalize(*DifMeig,eig);
//   cout << Det(*DifMeig);
//   cout << Transpose(*DifMeig)*(*Diffusion_matrix)*(*DifMeig);
//   cout.flush();
   for(int i=0;i<6;i++) {
     if(eig[2*i]<=0.) {
       cout << "Negative eigen values (set_beig)\n";
       exit(1);}
     if(fabs(eig[2*i+1])>1.e-30) {
       cout << "Imaginary eigen values (set_beig)\n";
       exit(1);}
     beig[i]=sqrt(eig[2*i]); 
   }
   (*Damping_matrix)=((matrix)*Damping_matrix)*SInverse(*Linear_map);
   *DifMeig_i_D=Transpose(*DifMeig)*(*Damping_matrix);
}

void EMIT::set_EMIT_N(void)
{
   matrix NtoP(6,6),I6(6,6);
   NtoP=SInverse(*PtoN);
   I6.I();
   //   cout << "comecome set_EMIT_N\n"; cout.flush();
   *Damping_matrix=I6-(*Damping_matrix);
   *DifMeig=NtoP;
   for(int i=0;i<6;i++) beig[i]=sqrt((*Diffusion_matrix)[i][i]);
   *DifMeig_i_D=((matrix)*Damping_matrix)*(*PtoN);
//   *DifMeig_i_D=*PtoN;
//   cout << "DifMeig_i_D \n" << *DifMeig_i_D;
}

void EMIT::DampRateReset(double drx,double dry,double drz,
			 double bx,double by,double bz)
{
  double k[6];
  matrix I6(6,6);
  I6.I();
  k[0]=drx; k[1]=drx; k[2]=dry; k[3]=dry; k[4]=drz; k[5]=drz;
  Damping_matrix->DiagonalMatrix(k);
  k[0]=bx; k[1]=bx; k[2]=by; k[3]=by; k[4]=bz; k[5]=bz;
  Diffusion_matrix->DiagonalMatrix(k);
  *Damping_matrix=I6-(*Damping_matrix);
  for(int i=0;i<6;i++) beig[i]=sqrt((*Diffusion_matrix)[i][i]);
  *DifMeig_i_D=((matrix)*Damping_matrix)*(*PtoN);
//   *DifMeig_i_D=*PtoN;
//   cout << "DifMeig_i_D \n" << *DifMeig_i_D;
}

ostream& operator<<(ostream& s, EMIT& A)
{
  int i;
  s << "\n##########  EMIT information ###############";
  s << "\n Linear Map \n" << *A.Linear_map;
  s << "\n Phys2Norm  \n" << *A.PtoN;
  s << "\n Damping matrix\n" << *A.Damping_matrix;
  s << "\n Diffusion matrix \n" << *A.Diffusion_matrix;
//  s << "\n Eigen matrix diagonalizing Diffusion matrix\n" << *A.DifMeig;
  s << "\n Norm2Phys\n" << *A.DifMeig;
  s << "\n Eigen values of the Diffusion matrix (sqrt)\n";
  for(i=0;i<6;i++) {s.width(12);  s << A.beig[i] << ' ';} s << '\n';
  s << "\n Beam envelope matrix \n" << *A.Beam_envelope;
  s << "\n Closed Orbit\n";
  for(i=0;i<6;i++) {s.width(12);  s << A.cod[i] << ' ';}
  s<< '\n';
  s.flush();
  return s;
}


void EMIT::SynchrotronRadiation(map_double& x)
{
  double x0,x1,x2,x3,x4,x5;

  x0=(*DifMeig_i_D[0][0])*x[0]+(*DifMeig_i_D[0][1])*x[3]+
     (*DifMeig_i_D[0][2])*x[1]+(*DifMeig_i_D[0][3])*x[4]+
     (*DifMeig_i_D[0][4])*x[2]+(*DifMeig_i_D[0][5])*x[5];
  x3=(*DifMeig_i_D[1][0])*x[0]+(*DifMeig_i_D[1][1])*x[3]+
     (*DifMeig_i_D[1][2])*x[1]+(*DifMeig_i_D[1][3])*x[4]+
     (*DifMeig_i_D[1][4])*x[2]+(*DifMeig_i_D[1][5])*x[5];
  x1=(*DifMeig_i_D[2][0])*x[0]+(*DifMeig_i_D[2][1])*x[3]+
     (*DifMeig_i_D[2][2])*x[1]+(*DifMeig_i_D[2][3])*x[4]+
     (*DifMeig_i_D[2][4])*x[2]+(*DifMeig_i_D[2][5])*x[5];
  x4=(*DifMeig_i_D[3][0])*x[0]+(*DifMeig_i_D[3][1])*x[3]+
     (*DifMeig_i_D[3][2])*x[1]+(*DifMeig_i_D[3][3])*x[4]+
     (*DifMeig_i_D[3][4])*x[2]+(*DifMeig_i_D[3][5])*x[5];
  x2=(*DifMeig_i_D[4][0])*x[0]+(*DifMeig_i_D[4][1])*x[3]+
     (*DifMeig_i_D[4][2])*x[1]+(*DifMeig_i_D[4][3])*x[4]+
     (*DifMeig_i_D[4][4])*x[2]+(*DifMeig_i_D[4][5])*x[5];
  x5=(*DifMeig_i_D[5][0])*x[0]+(*DifMeig_i_D[5][1])*x[3]+
     (*DifMeig_i_D[5][2])*x[1]+(*DifMeig_i_D[5][3])*x[4]+
     (*DifMeig_i_D[5][4])*x[2]+(*DifMeig_i_D[5][5])*x[5];

  x0+=beig[0]*tgauss();
  x3+=beig[1]*tgauss();
  x1+=beig[2]*tgauss();
  x4+=beig[3]*tgauss();
  x2+=beig[4]*tgauss();
  x5+=beig[5]*tgauss();

  x[0]=(*DifMeig[0][0])*x0+(*DifMeig[0][1])*x3+
       (*DifMeig[0][2])*x1+(*DifMeig[0][3])*x4+
       (*DifMeig[0][4])*x2+(*DifMeig[0][5])*x5;
  x[3]=(*DifMeig[1][0])*x0+(*DifMeig[1][1])*x3+
       (*DifMeig[1][2])*x1+(*DifMeig[1][3])*x4+
       (*DifMeig[1][4])*x2+(*DifMeig[1][5])*x5;
  x[1]=(*DifMeig[2][0])*x0+(*DifMeig[2][1])*x3+
       (*DifMeig[2][2])*x1+(*DifMeig[2][3])*x4+
       (*DifMeig[2][4])*x2+(*DifMeig[2][5])*x5;
  x[4]=(*DifMeig[3][0])*x0+(*DifMeig[3][1])*x3+
       (*DifMeig[3][2])*x1+(*DifMeig[3][3])*x4+
       (*DifMeig[3][4])*x2+(*DifMeig[3][5])*x5;
  x[2]=(*DifMeig[4][0])*x0+(*DifMeig[4][1])*x3+
       (*DifMeig[4][2])*x1+(*DifMeig[4][3])*x4+
       (*DifMeig[4][4])*x2+(*DifMeig[4][5])*x5;
  x[5]=(*DifMeig[5][0])*x0+(*DifMeig[5][1])*x3+
       (*DifMeig[5][2])*x1+(*DifMeig[5][3])*x4+
       (*DifMeig[5][4])*x2+(*DifMeig[5][5])*x5;
}

void EMIT::SynchrotronRadiation(pBeam& x)
{
  double x0,x1,x2,x3,x4,x5;

  for(int i=0;i<x.N_particle;i++) {
    x0=(*DifMeig_i_D)[0][0]*x.x[i]+(*DifMeig_i_D)[0][1]*x.px[i]+
       (*DifMeig_i_D)[0][2]*x.y[i]+(*DifMeig_i_D)[0][3]*x.py[i]+
       (*DifMeig_i_D)[0][4]*x.z[i]+(*DifMeig_i_D)[0][5]*x.pz[i];
    x3=(*DifMeig_i_D)[1][0]*x.x[i]+(*DifMeig_i_D)[1][1]*x.px[i]+
       (*DifMeig_i_D)[1][2]*x.y[i]+(*DifMeig_i_D)[1][3]*x.py[i]+
       (*DifMeig_i_D)[1][4]*x.z[i]+(*DifMeig_i_D)[1][5]*x.pz[i];
    x1=(*DifMeig_i_D)[2][0]*x.x[i]+(*DifMeig_i_D)[2][1]*x.px[i]+
       (*DifMeig_i_D)[2][2]*x.y[i]+(*DifMeig_i_D)[2][3]*x.py[i]+
       (*DifMeig_i_D)[2][4]*x.z[i]+(*DifMeig_i_D)[2][5]*x.pz[i];
    x4=(*DifMeig_i_D)[3][0]*x.x[i]+(*DifMeig_i_D)[3][1]*x.px[i]+
       (*DifMeig_i_D)[3][2]*x.y[i]+(*DifMeig_i_D)[3][3]*x.py[i]+
       (*DifMeig_i_D)[3][4]*x.z[i]+(*DifMeig_i_D)[3][5]*x.pz[i];
    x2=(*DifMeig_i_D)[4][0]*x.x[i]+(*DifMeig_i_D)[4][1]*x.px[i]+
       (*DifMeig_i_D)[4][2]*x.y[i]+(*DifMeig_i_D)[4][3]*x.py[i]+
       (*DifMeig_i_D)[4][4]*x.z[i]+(*DifMeig_i_D)[4][5]*x.pz[i];
    x5=(*DifMeig_i_D)[5][0]*x.x[i]+(*DifMeig_i_D)[5][1]*x.px[i]+
       (*DifMeig_i_D)[5][2]*x.y[i]+(*DifMeig_i_D)[5][3]*x.py[i]+
       (*DifMeig_i_D)[5][4]*x.z[i]+(*DifMeig_i_D)[5][5]*x.pz[i];

    x0+=beig[0]*tgauss();
    x3+=beig[1]*tgauss();
    x1+=beig[2]*tgauss();
    x4+=beig[3]*tgauss();
    x2+=beig[4]*tgauss();
    x5+=beig[5]*tgauss();
//    cout << x1 << " " << x4 << "\n";

    x.x[i]= (*DifMeig)[0][0]*x0+(*DifMeig)[0][1]*x3+
            (*DifMeig)[0][2]*x1+(*DifMeig)[0][3]*x4+
            (*DifMeig)[0][4]*x2+(*DifMeig)[0][5]*x5;
    x.px[i]=(*DifMeig)[1][0]*x0+(*DifMeig)[1][1]*x3+
            (*DifMeig)[1][2]*x1+(*DifMeig)[1][3]*x4+
            (*DifMeig)[1][4]*x2+(*DifMeig)[1][5]*x5;
    x.y[i]= (*DifMeig)[2][0]*x0+(*DifMeig)[2][1]*x3+
            (*DifMeig)[2][2]*x1+(*DifMeig)[2][3]*x4+
            (*DifMeig)[2][4]*x2+(*DifMeig)[2][5]*x5;
    x.py[i]=(*DifMeig)[3][0]*x0+(*DifMeig)[3][1]*x3+
            (*DifMeig)[3][2]*x1+(*DifMeig)[3][3]*x4+
            (*DifMeig)[3][4]*x2+(*DifMeig)[3][5]*x5;
    x.z[i]= (*DifMeig)[4][0]*x0+(*DifMeig)[4][1]*x3+
            (*DifMeig)[4][2]*x1+(*DifMeig)[4][3]*x4+
            (*DifMeig)[4][4]*x2+(*DifMeig)[4][5]*x5;
    x.pz[i]=(*DifMeig)[5][0]*x0+(*DifMeig)[5][1]*x3+
            (*DifMeig)[5][2]*x1+(*DifMeig)[5][3]*x4+
            (*DifMeig)[5][4]*x2+(*DifMeig)[5][5]*x5;
  }
}

void ReadBeam(char* bp,Beam& BEAM)
{
  double energy,emx,emy,emz,Nparticle;
  energy=get_parm(bp,"energy");
  emx=get_parm(bp,"emx");
  emy=get_parm(bp,"emy");
  emz=get_parm(bp,"emz");
  Nparticle=get_parm(bp,"N_particle");
  BEAM.SetEnergy(energy);
  BEAM.SetEmittance(emx,emy,emz);
  BEAM.SetNparticle(Nparticle);
}

void ReadMathList(char* s,matrix& M)
{
/* elist is stored PtoN */
  size_t i;
  char *p;
  double *elistR;

  p=strstr(s,"{{")+2;
  for(i=0;i<strlen(p);i++) {if(*(p+i)==',') *(p+i)=' ';}
  elistR=M.t;
  for(i=0;i<6;i++) {
    sscanf(p,"%lf %lf %lf %lf %lf %lf",elistR+i*6,elistR+1+i*6,elistR+2+i*6,
	   elistR+3+i*6,elistR+4+i*6,elistR+5+i*6);
    p=strchr(p,'{')+1;
  }
}


int ReadMathList(char* s,double* K)
{
/* elist is stored PtoN */
  int i;
  char *p,*p1;

  p=strchr(s,'{')+1;
  i=0;
  while((p1=strchr(p,','))!=NULL) {
    K[i]=atof(p);
    p=p1+1;
    i++;
  }
  K[i]=atof(p);
  return i+1;
}

